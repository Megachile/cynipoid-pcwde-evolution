library(ggtree)
library(ggplot2)
library(ape)

interface_file <- Sys.getenv("S3_INTERFACE")
tree_file <- Sys.getenv("S3_TREE_FILE")
output_dir <- Sys.getenv("S3_OUTPUT_DIR")
if (interface_file == "" || tree_file == "" || output_dir == "") {
  stop("S3_INTERFACE, S3_TREE_FILE, and S3_OUTPUT_DIR are required")
}
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

ann <- read.table(interface_file, sep="\t", header=TRUE, quote="",
                  comment.char="", stringsAsFactors=FALSE, check.names=FALSE)
ann$display_label <- gsub("Ceroptres ex Druon ignotum", "Ceroptres sp. ex Druon ignotum", ann$display_label, fixed=TRUE)
expected_tips <- as.integer(Sys.getenv("S3_EXPECTED_TIPS", "68"))
required <- c("tip_id", "biological_unit_id", "species", "display_label",
              "full_sequence_length_aa", "architecture_n_to_c",
              "copy_index_n_to_c", "synteny_style", "tip_class")
stopifnot(all(required %in% names(ann)), nrow(ann) == expected_tips,
          length(unique(ann$tip_id)) == expected_tips)

tr <- ladderize(read.tree(tree_file), right=FALSE)
stopifnot(setequal(tr$tip.label, ann$tip_id))

domain_colors <- c(B="#D98C00", B_other="#F4E3BB", GH30="#7570B3",
                   GH59="#5A4FCF", `PI-PLC-X`="#E41A1C")
block_colors <- c(chr6a="#D98C00", chr3a="#2166AC", chr3b="#2166AC",
                  sg_chr1="#A6761D", sg_chr3_novel="#2166AC",
                  sg_chr4="#009E73", outgroup="#000000", `-`="#777777")
block_display <- c(chr6a="chr6", chr3a="chr3", chr3b="chr3",
                   sg_chr1="chr1", sg_chr3_novel="chr3", sg_chr4="chr4",
                   outgroup="outgroup", `-`="synteny unresolved")

# Expand the full attached architecture once for each domain-tree terminal.
blocks <- list()
for (i in seq_len(nrow(ann))) {
  tokens <- strsplit(ann$architecture_n_to_c[i], ";", fixed=TRUE)[[1]]
  beta_index <- 0L
  for (token in tokens) {
    fields <- regmatches(token, regexec("^([^:]+):(\\d+)-(\\d+)$", token))[[1]]
    stopifnot(length(fields) == 4L)
    cls <- fields[2]
    if (cls == "B") beta_index <- beta_index + 1L
    focal <- cls == "B" && beta_index == as.integer(ann$copy_index_n_to_c[i])
    blocks[[length(blocks) + 1L]] <- data.frame(
      tip_id=ann$tip_id[i], class=cls, start=as.integer(fields[3]),
      end=as.integer(fields[4]), focal=focal,
      protein_length=as.integer(ann$full_sequence_length_aa[i]),
      stringsAsFactors=FALSE)
  }
}
dom <- do.call(rbind, blocks)
stopifnot(sum(dom$focal) == expected_tips)

# Every row is anchored on the center of the B-lectin domain represented by
# that tree terminal. Sibling lectin copies remain attached but are de-emphasized.
focal_centers <- setNames((dom$start[dom$focal] + dom$end[dom$focal]) / 2,
                          dom$tip_id[dom$focal])
stopifnot(setequal(names(focal_centers), ann$tip_id))

p <- ggtree(tr, branch.length="none", size=0.28) + scale_y_reverse()
tip_data <- p$data[p$data$isTip, ]
max_tip_x <- max(tip_data$x)
display_map <- setNames(ann$display_label, ann$tip_id)
block_map <- setNames(ann$synteny_style, ann$tip_id)
block_map[ann$tip_class == "bacterial"] <- "outgroup"
block_map[is.na(block_map) | block_map == ""] <- "-"

# Layout mirrors Figure 8: tree | terminal labels | complete protein architecture.
label_left_pad <- 1.2
label_x_start <- max_tip_x + label_left_pad
max_label_chars <- max(nchar(ann$display_label))
label_width <- max_label_chars * 0.36
domain_x_start <- label_x_start + label_width + 4.5
aa_per_unit <- 50
max_left_aa <- max(focal_centers[ann$tip_id] - 1)
anchor_x <- domain_x_start + max_left_aa / aa_per_unit
domain_height <- 0.34

tip_y <- setNames(tip_data$y, tip_data$label)
for (i in seq_len(nrow(ann))) {
  lab <- ann$tip_id[i]
  center <- focal_centers[lab]
  x0 <- anchor_x + (1 - center) / aa_per_unit
  x1 <- anchor_x + (as.integer(ann$full_sequence_length_aa[i]) - center) / aa_per_unit
  p <- p + annotate("segment", x=x0, xend=x1, y=tip_y[lab], yend=tip_y[lab],
                    color="grey70", linewidth=0.32)
}
for (i in seq_len(nrow(dom))) {
  lab <- dom$tip_id[i]
  center <- focal_centers[lab]
  xmin <- anchor_x + (dom$start[i] - center) / aa_per_unit
  xmax <- anchor_x + (dom$end[i] - center) / aa_per_unit
  fill <- if (dom$class[i] == "B") {
    if (dom$focal[i]) domain_colors["B"] else domain_colors["B_other"]
  } else if (dom$class[i] %in% names(domain_colors)) {
    domain_colors[dom$class[i]]
  } else "grey55"
  border <- if (dom$focal[i]) "#5C3A00" else "grey45"
  line_width <- if (dom$focal[i]) 0.42 else 0.14
  p <- p + annotate("rect", xmin=xmin, xmax=xmax,
                    ymin=tip_y[lab]-domain_height, ymax=tip_y[lab]+domain_height,
                    fill=fill, color=border, linewidth=line_width)
}

p$data$display_label <- ifelse(p$data$isTip, display_map[p$data$label], NA_character_)
p$data$block <- ifelse(p$data$isTip, block_map[p$data$label], NA_character_)
p <- p + geom_tiplab(aes(label=display_label, color=block), size=3.15,
                     offset=label_x_start-max_tip_x, fontface="bold", hjust=0) +
  scale_color_manual(values=block_colors, na.value="grey50")

# Show only supported nodes, as in the prior complete-domain supplemental tree.
numeric_support <- suppressWarnings(as.numeric(p$data$label))
p$data$support_label <- ifelse(!p$data$isTip & is.finite(numeric_support) &
                                numeric_support >= 70,
                                sprintf("%.0f", numeric_support), "")
p <- p + geom_text2(aes(subset=!isTip & support_label != "", label=support_label),
                    size=2.2, hjust=1.1, vjust=-0.35, color="grey25")

# Match the headed chromosome/domain key format of Figure 5.
n_tips <- nrow(ann)
legend_y <- n_tips + 3
used_raw <- intersect(names(block_display), unique(block_map))
seen <- character(); used <- character()
for (b in used_raw) {
  label <- block_display[b]
  if (!(label %in% seen)) { used <- c(used, b); seen <- c(seen, label) }
}
# Derive key membership from the exact marks expanded and drawn above.
drawn_class <- ifelse(dom$class == "B" & !dom$focal, "B_other", dom$class)
domain_labels <- c(B="tree-tip B-lectin domain", B_other="other B-lectin domain",
                   GH30="GH30", GH59="GH59", `PI-PLC-X`="PI-PLC-X")
used_domains <- intersect(names(domain_labels), unique(drawn_class))
stopifnot(setequal(used_domains, unique(drawn_class)))
domain_key <- data.frame(label=unname(domain_labels[used_domains]),
                         color=unname(domain_colors[used_domains]))
# Measure both groups at the saved type size so columns follow their labels.
right_edge <- max(vapply(seq_len(nrow(ann)), function(i) {
  anchor_x + (as.integer(ann$full_sequence_length_aa[i]) - focal_centers[ann$tip_id[i]]) / aa_per_unit
}, numeric(1)))
key_size <- 3.5
header_size <- 3.8
swatch_width <- 2.6
label_gap <- 0.65
column_gap <- 1.3
group_gap <- 3.0
# Resolve the same panel and coordinate expansion used for the saved figure.
measure_plot <- p + expand_limits(y=n_tips + 8.1) + xlim_tree(right_edge + 3) +
  coord_cartesian(xlim=c(-1, right_edge + 3), clip="off") +
  theme(legend.position="none", plot.margin=margin(5, 5, 5, 5))
grDevices::pdf(NULL, width=14, height=15.5)
grid::grid.newpage()
grid::grid.draw(ggplot2::ggplotGrob(measure_plot))
grid::grid.force()
viewports <- grid::grid.ls(viewports=TRUE, grobs=FALSE, print=FALSE)$name
panel_name <- viewports[grepl("^panel\\.", viewports)][1]
stopifnot(!is.na(panel_name))
grid::seekViewport(panel_name)
panel_inches <- grid::convertWidth(grid::unit(1, "npc"), "inches", valueOnly=TRUE)
panel_range <- ggplot2::ggplot_build(measure_plot)$layout$panel_params[[1]]$x.range
units_per_inch <- diff(panel_range) / panel_inches
measure_key <- function(labels, bold=FALSE, size=key_size) {
  vapply(labels, function(label) grid::convertWidth(
    grid::grobWidth(grid::textGrob(label, gp=grid::gpar(
      fontsize=size * ggplot2::.pt, fontface=if (bold) "bold" else "plain"))),
    "inches", valueOnly=TRUE), numeric(1)) *
    units_per_inch
}
block_widths <- measure_key(unname(block_display[used]), TRUE)
domain_widths <- measure_key(domain_key$label)
headers <- measure_key(c("Reference chromosome", "Protein domains"), TRUE, header_size)
grDevices::dev.off()
column_step <- function(widths) swatch_width + label_gap +
  max(widths[seq_along(widths) %% 2 == 1]) + column_gap
group_width <- function(widths) column_step(widths) + swatch_width +
  label_gap + max(widths[seq_along(widths) %% 2 == 0])
block_start <- 1
# Centre the protein-domain key under the drawn protein architectures
# (domain_x_start to right_edge), never closer than group_gap to the
# chromosome key.
domain_block_width <- max(group_width(domain_widths), headers[2])
domain_start <- max(
  block_start + max(group_width(block_widths), headers[1]) + group_gap,
  (domain_x_start + right_edge) / 2 - domain_block_width / 2)
stopifnot(domain_start + max(group_width(domain_widths), headers[2]) < right_edge + 3)
p <- p +
  annotate("text", x=block_start, y=legend_y-1.7, label="Reference chromosome",
           hjust=0, size=header_size, fontface="bold") +
  annotate("text", x=domain_start, y=legend_y-1.7, label="Protein domains",
           hjust=0, size=header_size, fontface="bold")
for (i in seq_along(used)) {
  b <- used[i]
  x <- block_start + ((i-1) %% 2) * column_step(block_widths)
  y <- legend_y + ((i-1) %/% 2) * 2.05
  p <- p + annotate("rect", xmin=x, xmax=x+swatch_width,
                    ymin=y-0.95, ymax=y+0.95,
                    fill=block_colors[b], color="grey30", linewidth=0.2) +
    annotate("text", x=x+swatch_width+label_gap, y=y, label=block_display[b],
             hjust=0, size=key_size, color=block_colors[b], fontface="bold")
}
for (i in seq_len(nrow(domain_key))) {
  x <- domain_start + ((i-1) %% 2) * column_step(domain_widths)
  y <- legend_y + ((i-1) %/% 2) * 2.05
  p <- p + annotate("rect", xmin=x, xmax=x+swatch_width,
                    ymin=y-0.95, ymax=y+0.95,
                    fill=domain_key$color[i], color="grey30", linewidth=0.2) +
    annotate("text", x=x+swatch_width+label_gap, y=y, label=domain_key$label[i],
             hjust=0, size=key_size)
}

# Preserve the original tree panel's vertical scale when an undrawn key
# category is removed: the old third key row ended at n_tips + 8.1.
p <- p + expand_limits(y=n_tips + 8.1) + xlim_tree(right_edge + 3) +
  coord_cartesian(xlim=c(-1, right_edge + 3), clip="off") +
  theme(legend.position="none", plot.margin=margin(5, 5, 5, 5))

base <- file.path(output_dir, "Figure_S4")
ggsave(paste0(base, ".svg"), p, width=14, height=15.5, limitsize=FALSE)

visual <- tip_data$label[order(tip_data$y)]
order_table <- ann[match(visual, ann$tip_id), ]
order_table$display_order <- seq_len(nrow(order_table))
write.table(order_table, paste0(base, "_terminal_order.tsv"), sep="\t",
            quote=FALSE, row.names=FALSE)
qc <- data.frame(
  metric=c("tree_tips", "interface_rows", "unique_tip_ids", "highlighted_tip_domains",
           "architecture_units", "display_geometry", "branch_lengths_displayed",
           "support_labels_threshold", "rooting", "figure_width_inches", "figure_height_inches"),
  value=c(length(tr$tip.label), nrow(ann), length(unique(ann$tip_id)), sum(dom$focal),
          length(unique(ann$biological_unit_id)), "full_architecture",
          "no", ">=70 UFBoot", "display-only bacterial-domain MRCA", 14, 15.5),
  status="PASS", stringsAsFactors=FALSE)
write.table(qc, paste0(base, "_render_qc.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("Wrote Figure-8-style Figure S6 with", nrow(ann), "domain tips\n")

write.table(dom, file.path(output_dir, "Figure_S4_rendered_domains.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(p$data, file.path(output_dir, "Figure_S4_plot_data.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Figure 7 B-lectin locus tree with domain
# architecture strip + tip labels colored by reference (BK) chromosome from the built copy table.
#
# Alignment mode: the inferred ancestral D2 lectin is the shared anchor.
# D1-like singlets are shifted to the D1 lane; D2-like singlets to D2.

library(ggtree)
library(treeio)
library(ggplot2)
library(ape)

# ── Domain colors/names ──
domain_colors <- c(
  A="#1B7837", B="#E6AB02", C="#A6D854", D="#E7298A", E="#984EA3",
  F="#FF7F00", G="#8DD3C7", H="#252525", I="#BEBADA", J="#80B1D3",
  K="#542788", L="#2166AC", M="#D98C00", N="#D6604D", O="#C2185B",
  P="#B2182B", Q="#E31A1C", R="#F781BF", S="#FFD92F", T="#A6761D",
  U="#7570B3", V="#BEAED4", W="#66C2A5", X="#FDB462", Y="#7FC97F",
  Z="#BDBDBD",
  u="#5A4FCF", v="#9897E0", w="#40A28C", y="#B8D88C", c="#006D2C",
  s="#525252", p="#E41A1C",
  m="#8B5A00", q="#7A0D16"
)

domain_names <- c(
  A="GH5_2", B="CBM2", C="CBM3", D="CBM5", E="CBM6", F="CBM10",
  G="CBM12", H="CBM17", I="CBM28", J="CBM73", K="CBM81", L="PL1_8",
  M="B-lectin domain", N="PL4", O="PL4_2", P="PL4_4", Q="PL4_5", R="CBM22",
  S="CBM32", T="CBM64", U="GH30", V="GH30_C", W="CBM4_9", X="fn3",
  Y="Ricin_B", Z="Big_5",
  u="GH59", v="GH30_2", w="CBM6", y="RicinB_2", c="CBM35_2",
  s="SLH", p="PI-PLC-X",
  m="B-lectin domain (weak)", q="Region A"
)

# Tip colour = the copy's reference (BK) chromosome in the built copy table
# (scripts/members/build_copy_table.py), the same placement Figure 1 and
# Tables S2/S3 use. chr3/chr6 keep their earlier colours.
block_colors <- c(
  chr1       = "#A6761D", chr2 = "#E69F00", chr3 = "#56B4E9", chr4 = "#CC79A7",
  chr5       = "#0072B2", chr6 = "#009E73", chr7 = "#D55E00", chr8 = "#F0E442",
  chr9       = "#882255", chr10 = "#44AA99",
  unplaced   = "#8C8C8C",
  transcript = "#BDBDBD",
  outgroup   = "#000000"
)

block_legend_order <- c(paste0("chr", 1:10), "unplaced", "transcript", "outgroup")

# ── Parse domain FAA ──
parse_domains <- function(faa_file) {
  lines <- readLines(faa_file)
  result <- list()
  cur_name <- NULL; cur_seq <- ""
  process_seq <- function(name, seq) {
    if (is.null(name) || nchar(seq) == 0) return(NULL)
    chars <- strsplit(seq, "")[[1]]
    blocks <- list(); cur_ch <- NULL; cur_start <- NULL
    for (i in seq_along(chars)) {
      ch <- chars[i]
      if (ch != "-") {
        if (is.null(cur_ch) || ch != cur_ch) {
          if (!is.null(cur_ch)) {
            blocks[[length(blocks)+1]] <- data.frame(
              label=name, start=cur_start, end=i-1, domain=cur_ch,
              prot_len=length(chars), stringsAsFactors=FALSE)
          }
          cur_ch <- ch; cur_start <- i
        }
      } else if (!is.null(cur_ch)) {
        blocks[[length(blocks)+1]] <- data.frame(
          label=name, start=cur_start, end=i-1, domain=cur_ch,
          prot_len=length(chars), stringsAsFactors=FALSE)
        cur_ch <- NULL
      }
    }
    if (!is.null(cur_ch)) {
      blocks[[length(blocks)+1]] <- data.frame(
        label=name, start=cur_start, end=length(chars), domain=cur_ch,
        prot_len=length(chars), stringsAsFactors=FALSE)
    }
    if (length(blocks) > 0) do.call(rbind, blocks) else NULL
  }
  for (line in lines) {
    line <- trimws(line)
    if (startsWith(line, ">")) {
      df <- process_seq(cur_name, cur_seq)
      if (!is.null(df)) result[[length(result)+1]] <- df
      cur_name <- sub("^>", "", line); cur_seq <- ""
    } else if (nchar(line) > 0) {
      cur_seq <- paste0(cur_seq, line)
    }
  }
  df <- process_seq(cur_name, cur_seq)
  if (!is.null(df)) result[[length(result)+1]] <- df
  do.call(rbind, result)
}

# ── Uniform-beta input adapter ──
# The historical visual implementation below remains authoritative. This
# adapter replaces only its retired roster/tree/domain input layer.
interface_file <- Sys.getenv("PCWDE_B_LECTIN_TREE_INTERFACE")
tree_file <- Sys.getenv("PCWDE_B_LECTIN_TREE_TREE_FILE")
alignment_file <- Sys.getenv("PCWDE_B_LECTIN_TREE_ALIGNMENT_FILE")
output_dir <- Sys.getenv("PCWDE_B_LECTIN_TREE_OUTPUT_DIR", "")
figure_width <- as.numeric(Sys.getenv("PCWDE_B_LECTIN_TREE_WIDTH", "6.89"))
figure_height <- as.numeric(Sys.getenv("PCWDE_B_LECTIN_TREE_HEIGHT", "8.6"))
expected_tips <- as.integer(Sys.getenv("PCWDE_B_LECTIN_TREE_EXPECTED_TIPS", "48"))
expected_d1_singletons <- as.integer(Sys.getenv("PCWDE_B_LECTIN_TREE_EXPECTED_D1_SINGLETONS", "3"))
expected_d2_singletons <- as.integer(Sys.getenv("PCWDE_B_LECTIN_TREE_EXPECTED_D2_SINGLETONS", "28"))
if (interface_file == "" || tree_file == "" || alignment_file == "") {
  stop("PCWDE_B_LECTIN_TREE_INTERFACE, PCWDE_B_LECTIN_TREE_TREE_FILE, and PCWDE_B_LECTIN_TREE_ALIGNMENT_FILE are required")
}
if (!nzchar(output_dir)) stop("PCWDE_B_LECTIN_TREE_OUTPUT_DIR is required")
if (dir.exists(output_dir) && length(list.files(output_dir, all.files=TRUE, no..=TRUE))) {
  stop("refusing non-empty output directory: ", output_dir)
}
ann <- read.table(interface_file, sep="\t", header=TRUE, quote="",
                  comment.char="", stringsAsFactors=FALSE,
                  check.names=FALSE)
required_columns <- c(
  "biological_unit_id", "species", "full_sequence_length_aa",
  "one_best_tip_id", "architecture_n_to_c", "synteny_style", "tip_role"
)
stopifnot(all(required_columns %in% names(ann)))
stopifnot(nrow(ann) == expected_tips, length(unique(ann$one_best_tip_id)) == expected_tips)

tr_input <- read.tree(tree_file)
stopifnot(setequal(tr_input$tip.label, ann$one_best_tip_id))
out_tips <- ann$one_best_tip_id[ann$tip_role == "outgroup"]
stopifnot(length(out_tips) == 1L)
tr_rooted <- root(tr_input, outgroup=out_tips, resolve.root=TRUE)
tr <- ladderize(tr_rooted, right=FALSE)

# Preserve every input support label in a diagnostic, and draw only strong
# joint SH-aLRT/UFBoot support. The rooted display may add an unlabeled node.
support_label <- tr$node.label
support_match <- regexec("^([0-9]+[.]?[0-9]*)/([0-9]+[.]?[0-9]*)$", support_label)
support_parts <- regmatches(support_label, support_match)
support_a <- vapply(support_parts, function(x) if (length(x) == 3L) as.numeric(x[2]) else NA_real_, numeric(1))
support_b <- vapply(support_parts, function(x) if (length(x) == 3L) as.numeric(x[3]) else NA_real_, numeric(1))
support_manifest <- data.frame(node=seq_along(support_label) + length(tr$tip.label),
  original_label=support_label, SH_aLRT=support_a, UFBoot=support_b,
  starred=!is.na(support_a) & !is.na(support_b) & support_a >= 80 & support_b >= 95)
if (any(is.na(support_a) & !is.na(support_label) & nzchar(support_label) & support_label != "Root")) {
  warning("missing or malformed B-lectin support labels: ",
          paste(support_manifest$node[is.na(support_a) & !is.na(support_label) & nzchar(support_label) & support_label != "Root"], collapse=", "))
}

domain_codes <- c(B="M", GH30="U", GH59="u", `PI-PLC-X`="p")
uniform_blocks <- list()
for (i in seq_len(nrow(ann))) {
  tokens <- strsplit(ann$architecture_n_to_c[i], ";", fixed=TRUE)[[1]]
  for (token in tokens) {
    match <- regexec("^([^:]+):(\\d+)-(\\d+)$", token)
    fields <- regmatches(token, match)[[1]]
    stopifnot(length(fields) == 4L, fields[2] %in% names(domain_codes))
    uniform_blocks[[length(uniform_blocks)+1]] <- data.frame(
      label=ann$one_best_tip_id[i], start=as.integer(fields[3]),
      end=as.integer(fields[4]), domain=unname(domain_codes[fields[2]]),
      prot_len=as.integer(ann$full_sequence_length_aa[i]),
      stringsAsFactors=FALSE
    )
  }
}
dom_df <- do.call(rbind, uniform_blocks)
stopifnot(setequal(unique(dom_df$label), ann$one_best_tip_id))
cat("Domain blocks:", nrow(dom_df), "\n")

# Visible labels are parsed scientific species names from the accepted
# interface. Repeated species names are intentional; unique identifiers stay
# in the terminal-order sidecar rather than being appended to the figure.
ann$display_label <- ann$species
ann$display_label <- gsub("Eschatocerus acacia", "Eschatocerus acaciae", ann$display_label, fixed=TRUE)
ann$display_label <- gsub("Isocolus centaureae", "Eubothrus centaureae", ann$display_label, fixed=TRUE)
ann$display_label <- gsub("Ceroptres masudai", "Ceroptres kovalevi", ann$display_label, fixed=TRUE)
ann$display_label <- gsub("Bradysia coprophila", "Bradysia tilicola", ann$display_label, fixed=TRUE)
plotmath_quote <- function(x) encodeString(x, quote='"')
scientific_plotmath_bold <- function(x) {
  if (is.na(x) || !nzchar(x)) return('""')
  if (grepl(" ex ", x, fixed=TRUE)) {
    parts <- strsplit(x, " ex ", fixed=TRUE)[[1]]
    return(paste0("bolditalic(", plotmath_quote(parts[1]),
                  ")~bold(' ex ')~bolditalic(", plotmath_quote(parts[2]), ")"))
  }
  tokens <- strsplit(x, " +")[[1]]
  if (length(tokens) >= 2 && tokens[2] %in% c("sp", "sp.")) {
    suffix <- if (length(tokens) > 2) paste(tokens[-c(1,2)], collapse=" ") else ""
    answer <- paste0("bolditalic(", plotmath_quote(tokens[1]), ")~bold(' sp.')")
    if (nzchar(suffix)) answer <- paste0(answer, "~bold(", plotmath_quote(paste0(" ", suffix)), ")")
    return(answer)
  }
  if (length(tokens) >= 2) {
    suffix <- if (length(tokens) > 2) paste(tokens[-c(1,2)], collapse=" ") else ""
    answer <- paste0("bolditalic(", plotmath_quote(paste(tokens[1:2], collapse=" ")), ")")
    if (nzchar(suffix)) answer <- paste0(answer, "~bold(", plotmath_quote(paste0(" ", suffix)), ")")
    return(answer)
  }
  paste0("bolditalic(", plotmath_quote(x), ")")
}
ann$display_label_math <- vapply(
  ann$display_label, scientific_plotmath_bold, character(1)
)
stopifnot(all(!is.na(ann$display_label) & nzchar(trimws(ann$display_label))))
stopifnot(length(unique(ann$display_label)) < nrow(ann))
stopifnot(!any(grepl("\\[[0-9]+/[0-9]+\\]", ann$display_label)))

# Lookups: accepted tip -> visible species label, accepted tip -> reference chromosome.
display_map <- setNames(ann$display_label, ann$one_best_tip_id)
display_math_map <- setNames(ann$display_label_math, ann$one_best_tip_id)
species_map <- setNames(ann$species, ann$one_best_tip_id)
copy_table_file <- Sys.getenv("PCWDE_COPY_TABLE")
if (!nzchar(copy_table_file)) stop("PCWDE_COPY_TABLE (built copy table) is required")
copies <- read.table(copy_table_file, sep="\t", header=TRUE, quote="",
                     comment.char="", stringsAsFactors=FALSE, check.names=FALSE)
copies <- copies[copies$family == "B-lectin", , drop=FALSE]
copy_chromosome <- character()
for (i in seq_len(nrow(copies))) {
  ids <- trimws(c(copies$locus_protein_id[i], strsplit(copies$protein_ids[i], ",")[[1]]))
  copy_chromosome[ids[nzchar(ids)]] <- copies$bk_chromosome[i]
}
placement_of <- function(row) {
  if (row[["tip_role"]] == "outgroup") return("outgroup")
  if (startsWith(row[["biological_unit_id"]], "transcript|")) return("transcript")
  for (key in c(row[["one_best_tip_id"]], row[["biological_unit_id"]])) {
    if (!is.na(copy_chromosome[key])) return(unname(copy_chromosome[key]))
  }
  stop("B-lectin genome tip is not a copy-table copy: ", row[["biological_unit_id"]])
}
block_map <- setNames(vapply(seq_len(nrow(ann)), function(i) placement_of(ann[i, ]), ""),
                      ann$one_best_tip_id)
block_map[out_tips] <- "outgroup"
stopifnot(all(block_map %in% names(block_colors)))

# ── D1/D2-aware architecture anchor ──
alignment <- read.table(alignment_file, sep="\t", header=TRUE, quote="",
                        comment.char="", stringsAsFactors=FALSE,
                        check.names=FALSE)
stopifnot(nrow(alignment) == expected_tips,
          setequal(alignment$locus_id, ann$one_best_tip_id),
          all(alignment$alignment_lane %in% c("D1", "D2")))
alignment <- alignment[match(ann$one_best_tip_id, alignment$locus_id), ]

lectin_only <- dom_df[dom_df$domain %in% c("M", "m"), ]
lectin_blocks <- split(lectin_only, lectin_only$label)
lectin_centers <- lapply(lectin_blocks, function(g) (g$start + g$end) / 2)
anchors <- setNames(numeric(nrow(alignment)), alignment$locus_id)
for (i in seq_len(nrow(alignment))) {
  lab <- alignment$locus_id[i]
  idx <- as.integer(alignment$anchor_copy_index_n_to_c[i])
  stopifnot(length(lectin_centers[[lab]]) >= idx)
  anchors[lab] <- lectin_centers[[lab]][idx]
}
paired_gaps <- numeric()
for (lab in names(lectin_centers)) {
  if (length(lectin_centers[[lab]]) == 2L) {
    paired_gaps <- c(paired_gaps, diff(lectin_centers[[lab]]))
  }
}
lane_gap_aa <- median(paired_gaps)
stopifnot(is.finite(lane_gap_aa), lane_gap_aa > 0)
target_offsets <- setNames(ifelse(alignment$alignment_lane == "D1", -lane_gap_aa, 0),
                           alignment$locus_id)
prot_lens <- sapply(split(dom_df, dom_df$label), function(g) g$prot_len[1])

# ── Build ggtree ──
p <- ggtree(tr, branch.length="none", size=0.3)
p <- p + scale_y_reverse()

tip_data <- p$data[p$data$isTip, ]
max_tip_x <- max(tip_data$x)

# Units note: ggtree x is in "branch-count units" (1 = one edge step).
aa_per_unit     <- 110     # 1 x-unit of domain strip == 110 aa of protein
domain_height   <- 0.38
units_per_char  <- 0.55    # heuristic x-units per character at tip_size 4.5 (tune)
label_left_pad  <- 1
label_right_pad <- 9
# Retain the historical physical gutters on the explicit 170-mm canvas. Text
# size remains unchanged, so coordinate-space gutters scale inversely by width.
canvas_scale <- 1
label_left_pad <- label_left_pad * canvas_scale
label_right_pad <- label_right_pad * canvas_scale
units_per_char <- units_per_char * canvas_scale

# ── Layout: tree | LABELS | DOMAINS ──
# Need display labels NOW so we can compute the label column width.
tip_labels_for_width <- tip_data$label
tip_disp_for_width   <- unname(display_map[tip_labels_for_width])
stopifnot(!any(is.na(tip_disp_for_width)))
max_label_chars  <- max(nchar(tip_disp_for_width), na.rm = TRUE)
label_x_start    <- max_tip_x + label_left_pad
label_width      <- max_label_chars * units_per_char
domain_x_start   <- label_x_start + label_width + label_right_pad
max_left_aa      <- max(anchors - 1 - target_offsets[names(anchors)], na.rm=TRUE)
anchor_x         <- domain_x_start + max_left_aa / aa_per_unit

# Protein backbones and domains use the selected ancestral-domain center.
# D2 centers land at anchor_x; D1-only singlets land one median paired-domain
# spacing to the left. The grey backbone therefore passes through each domain
# box at its vertical center and preserves the original protein geometry.
for (i in seq_len(nrow(tip_data))) {
  lab <- tip_data$label[i]
  a <- anchors[lab]
  target <- target_offsets[lab]
  plen <- prot_lens[lab]
  x0 <- anchor_x + target / aa_per_unit + (1 - a) / aa_per_unit
  x1 <- anchor_x + target / aa_per_unit + (plen - a) / aa_per_unit
  p <- p + annotate("segment",
    x=x0, xend=x1,
    y=tip_data$y[i], yend=tip_data$y[i],
    color="grey70", linewidth=0.35)
}

for (i in seq_len(nrow(dom_df))) {
  lab <- dom_df$label[i]
  tip_row <- tip_data[tip_data$label == lab, ]
  if (nrow(tip_row) == 0) next
  tip_y <- tip_row$y[1]
  a <- anchors[lab]
  target <- target_offsets[lab]
  xmin <- anchor_x + target / aa_per_unit + (dom_df$start[i] - a) / aa_per_unit
  xmax <- anchor_x + target / aa_per_unit + (dom_df$end[i] - a) / aa_per_unit
  col <- domain_colors[dom_df$domain[i]]
  if (is.na(col)) col <- "grey50"
  p <- p + annotate("rect",
    xmin=xmin, xmax=xmax,
    ymin=tip_y - domain_height, ymax=tip_y + domain_height,
    fill=col, color="grey40", linewidth=0.15)
}

# ── Display labels: parsed species-only names from accepted interface ──
get_display <- function(label) {
  if (is.na(label)) return(NA_character_)
  unname(display_map[label])
}

p$data$display_label <- sapply(p$data$label, get_display)
p$data$display_label_math <- unname(display_math_map[p$data$label])

p$data$block <- ifelse(p$data$isTip, block_map[p$data$label], NA_character_)
p$data$block[p$data$label %in% out_tips] <- "outgroup"
stopifnot(!any(is.na(p$data$block[p$data$isTip])))

# Layout swap: labels go BEFORE the domain strip.
# `x=X` in geom_tiplab is ignored by this ggtree version; use offset=.
label_x    <- label_x_start
tip_offset <- label_x_start - max_tip_x

p <- p + geom_tiplab(aes(label=display_label_math, color=block),
                     size=3.5, offset=tip_offset, parse=TRUE,
                     hjust=0) +
  scale_color_manual(values=block_colors, na.value="grey50",
                     breaks=block_legend_order,
                     name="Reference chromosome")

# ── Joint-support asterisks ──
p$data$support_star <- ifelse(p$data$node %in% support_manifest$node[support_manifest$starred], "*", "")
p <- p + geom_nodelab(aes(label=support_star),
                      size=3.0, hjust=1.9, vjust=-0.4,
                      color="grey30")

# ── Legends below the figure: domain key on the left, block key on the right ──
# We abandon ggplot's built-in color legend (so both keys sit together)
# and draw both as manual annotations in the y-extent BELOW the tree
# (negative y relative to tip rows — because scale_y_reverse() is active,
# "below" in the final plot = y greater than n_tips).
tree_domains <- sort(unique(dom_df$domain))
n_tips <- sum(p$data$isTip)
legend_y_base <- n_tips + 3

# Keys side-by-side under the tree: block key on the LEFT, dom key RIGHT.
# Side-by-side keeps the tree from getting vertically squeezed (stacked
# keys took ~20 extra rows below the tree and shrunk tip-label row height
# past font size). Per-entry spacing is wide enough for the longest
# labels ("outgroup", "B-lectin domain (weak)").
blk_x_start <- 1
blk_row_h   <- 1.9
dom_row_h   <- 1.9

# Domain key: 2 cols, right half of legend area
dom_cols    <- 2
dom_spacing <- 28
n_dom       <- length(tree_domains)
n_dom_rows  <- ceiling(n_dom / dom_cols)
dom_x_start <- 43    # right of the block key

p <- p +
  annotate("text", x=blk_x_start, y=legend_y_base - 1.5,
           label="Reference chromosome", hjust=0, size=3.0, fontface="bold") +
  annotate("text", x=dom_x_start, y=legend_y_base - 1.5,
           label="Protein domains", hjust=0, size=3.0, fontface="bold")

for (idx in seq_along(tree_domains)) {
  letter <- tree_domains[idx]
  col_i <- (idx - 1) %% dom_cols
  row_i <- (idx - 1) %/% dom_cols
  x <- dom_x_start + col_i * dom_spacing
  y <- legend_y_base + row_i * dom_row_h
  dname <- domain_names[letter]; if (is.na(dname)) dname <- letter
  dcol  <- domain_colors[letter]; if (is.na(dcol))  dcol  <- "grey50"
  p <- p + annotate("rect",
    xmin=x, xmax=x + 2.0, ymin=y - 0.85, ymax=y + 0.85, fill=dcol,
    color="grey30", linewidth=0.2) +
    annotate("text", x=x + 2.8, y=y, label=dname,
             hjust=0, size=2.8)
}

block_display <- setNames(names(block_colors), names(block_colors))

# Reference-chromosome key (2 columns).
# De-duplicate by display label: all raw chromosome-3 sub-block styles use
# one display shade and one truthful "chr3" key entry.  The accepted interface
# and terminal-order sidecars retain the raw chr3a/chr3b/sg_chr3_novel values.
used_blocks_raw <- intersect(block_legend_order,
                             unique(p$data$block[p$data$isTip]))
seen_labels <- character()
used_blocks <- character()
for (b in used_blocks_raw) {
  lbl <- block_display[b]
  if (!is.na(lbl) && !(lbl %in% seen_labels)) {
    used_blocks <- c(used_blocks, b)
    seen_labels <- c(seen_labels, lbl)
  }
}
# Block key sits LEFT of the domain key (side-by-side).
blk_cols    <- 2
blk_spacing <- 20
blk_y_start <- legend_y_base

for (idx in seq_along(used_blocks)) {
  blk <- used_blocks[idx]
  col_i <- (idx - 1) %% blk_cols
  row_i <- (idx - 1) %/% blk_cols
  x <- blk_x_start + col_i * blk_spacing
  y <- blk_y_start + row_i * blk_row_h
  lbl <- block_display[blk]; if (is.na(lbl)) lbl <- blk
  col <- block_colors[blk]; if (is.na(col)) col <- "grey50"
  p <- p + annotate("rect",
    xmin=x, xmax=x + 2.0, ymin=y - 0.85, ymax=y + 0.85, fill=col,
    color="grey30", linewidth=0.2) +
    annotate("text", x=x + 2.8, y=y, label=lbl,
             hjust=0, size=2.8, color=col, fontface="bold")
}

# ── x-limits — labels on the LEFT of domains, so right edge is domain_right_edge ──
domain_right_edge <- max(vapply(seq_len(nrow(ann)), function(i) {
  lab <- ann$one_best_tip_id[i]
  anchor_x + target_offsets[lab] / aa_per_unit + (prot_lens[lab] - anchors[lab]) / aa_per_unit
}, numeric(1)))
right_pad <- 1
left_limit <- -1.5 * canvas_scale
p <- p + xlim_tree(domain_right_edge + right_pad) +
  ggplot2::coord_cartesian(xlim=c(left_limit, domain_right_edge + right_pad), clip="off") +
  theme(legend.position="none",
        plot.margin = margin(5, 5, 5, 5))

# Measure every label-to-domain gutter using the actual bold font metrics and
# the rendered panel viewport. This is a QC sidecar and does not alter geometry.
measure_device <- tempfile(fileext=".png")
grDevices::png(measure_device, width=figure_width, height=figure_height,
               units="in", res=300, bg="white")
grid::grid.newpage()
grid::grid.draw(ggplot2::ggplotGrob(p))
grid::grid.force()
viewport_names <- grid::grid.ls(viewports=TRUE, grobs=FALSE,
                                print=FALSE)$name
panel_names <- viewport_names[grepl("^panel\\.", viewport_names)]
stopifnot(length(panel_names) >= 1L)
grid::seekViewport(panel_names[1])
panel_viewport <- grid::current.viewport()
panel_width_inches <- grid::convertWidth(grid::unit(1, "npc"), "inches",
                                         valueOnly=TRUE)
panel_x_units_per_inch <- diff(panel_viewport$xscale) / panel_width_inches
label_font_size_points <- 3.5 * ggplot2::.pt
label_width_inches <- vapply(ann$display_label, function(label) {
  label_grob <- grid::textGrob(
    label, gp=grid::gpar(fontsize=label_font_size_points, fontface="bold")
  )
  grid::convertWidth(grid::grobWidth(label_grob), "inches", valueOnly=TRUE)
}, numeric(1))
grDevices::dev.off()
unlink(measure_device)

tip_x_map <- setNames(tip_data$x, tip_data$label)
domain_left_map <- vapply(split(dom_df, dom_df$label), function(g) {
  lab <- g$label[1]
  a <- anchors[lab]
  min(anchor_x + target_offsets[lab] / aa_per_unit + (g$start - a) / aa_per_unit)
}, numeric(1))
label_x_by_tip <- unname(tip_x_map[ann$one_best_tip_id]) + tip_offset
domain_left_by_tip <- unname(domain_left_map[ann$one_best_tip_id])
label_right_by_tip <- label_x_by_tip + label_width_inches * panel_x_units_per_inch
clearance_table <- data.frame(
  tip_id=ann$one_best_tip_id,
  stable_id=ann$biological_unit_id,
  species=ann$species,
  label_x=label_x_by_tip,
  label_width_inches=label_width_inches,
  label_right_x=label_right_by_tip,
  domain_left_x=domain_left_by_tip,
  clearance_x=domain_left_by_tip - label_right_by_tip,
  panel_width_inches=panel_width_inches,
  stringsAsFactors=FALSE
)
stopifnot(nrow(clearance_table) == expected_tips,
          all(is.finite(clearance_table$clearance_x)),
          all(clearance_table$clearance_x > 0))

# ── Save ── (historical 12×9.5 in canvas; 170-mm width is explicit)
output_basename <- "b_lectin_tree"
render_qc <- file.path(output_dir, paste0(output_basename, "_render_qc.tsv"))
dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
out_base <- file.path(output_dir, output_basename)
write.table(support_manifest, paste0(out_base, "_support_manifest.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE, na="")
cat("Starred internal nodes:", sum(support_manifest$starred), "/", nrow(support_manifest), "\n")
ggsave(paste0(out_base, ".pdf"), p, width=figure_width, height=figure_height, limitsize=FALSE)
ggsave(paste0(out_base, ".svg"), p, width=figure_width, height=figure_height, limitsize=FALSE)
ggsave(paste0(out_base, ".png"), p, width=figure_width, height=figure_height,
       dpi=300, bg="white", limitsize=FALSE)

# ggtree y increases downward after scale_y_reverse(), so ascending raw y is
# the effective visual top-to-bottom order.
visual_tip_data <- tip_data[order(tip_data$y), ]
visual_order <- visual_tip_data$label
stopifnot(length(visual_order) == expected_tips, visual_order[1] == out_tips)
order_table <- ann[match(visual_order, ann$one_best_tip_id), ]
order_table$row <- seq_len(nrow(order_table))
write.table(
  order_table[, c("row", "one_best_tip_id", "biological_unit_id", "species",
                  "display_label", "architecture_n_to_c")],
  paste0(out_base, "_terminal_order.tsv"), sep="\t", quote=FALSE,
  row.names=FALSE, col.names=c("row", "tip_id", "stable_id", "species",
                              "display_label", "architecture_n_to_c")
)
write.table(
  clearance_table[match(visual_order, clearance_table$tip_id), ],
  paste0(out_base, "_label_clearance.tsv"), sep="\t", quote=FALSE,
  row.names=FALSE
)
write.tree(tr_rooted, paste0(out_base, "_rooted_pre_ladderize.treefile"))
write.tree(tr, paste0(out_base, "_rooted_ladderized.treefile"))

qc <- data.frame(
  metric=c("tree_tips", "manifest_tips", "outgroup_first",
           "terminal_order_rows", "unique_tip_ids", "unique_stable_ids",
           "unique_display_labels", "minimum_label_clearance_x", "ordinal_labels",
           "independent_row_sorting", "root_operation", "ladderize_right",
           "visual_y_direction", "figure_width_inches", "figure_height_inches",
           "displayed_domains", "displayed_synteny"),
  value=c(length(visual_order), nrow(ann), visual_order[1] == out_tips,
          length(visual_order), length(unique(ann$one_best_tip_id)),
          length(unique(ann$biological_unit_id)), length(unique(ann$display_label)),
          min(clearance_table$clearance_x),
          sum(grepl("\\[[0-9]+/[0-9]+\\]", ann$display_label)), "no",
          "ape::root(resolve.root=TRUE)", "FALSE", "scale_y_reverse",
          figure_width, figure_height, paste(sort(unique(dom_df$domain)), collapse=","),
          paste(sort(unique(block_map)), collapse=",")),
  stringsAsFactors=FALSE
)
write.table(qc, render_qc, sep="\t", quote=FALSE, row.names=FALSE)
# Drawn tip colours, checked against the copy table by scripts/release/check_copy_counts.py.
dir.create(file.path(output_dir, "evidence"), showWarnings=FALSE)
write.table(data.frame(one_best_tip_id=ann$one_best_tip_id,
                       biological_unit_id=ann$biological_unit_id,
                       tip_role=ann$tip_role,
                       drawn_placement=unname(block_map[ann$one_best_tip_id])),
            file.path(output_dir, "evidence", "drawn_placements.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
alignment$rendered_anchor_center_x <- anchor_x + target_offsets[alignment$locus_id] / aa_per_unit
alignment$lane_gap_aa <- lane_gap_aa
write.table(alignment, paste0(out_base, "_domain_lane_alignment.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
stopifnot(visual_order[1] == out_tips,
          sum(alignment$singleton_classification == "D1-like") == expected_d1_singletons,
          sum(alignment$singleton_classification == "D2-like") == expected_d2_singletons)
cat("Wrote", out_base, ".{pdf,svg,png} and exact ladderized order sidecars\n")

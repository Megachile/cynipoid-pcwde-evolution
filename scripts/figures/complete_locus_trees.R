library(ggtree)
library(treeio)
library(ggplot2)
library(ape)
library(patchwork)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])), "locus_tip_marks.R"))

required_env <- function(name) {
  value <- Sys.getenv(name, unset="")
  if (!nzchar(value)) stop(name, " is required")
  value
}
TREE_DIR <- required_env("PCWDE_COMPLETE_LOCUS_TREES_INPUT_DIR")
FIG_OUT <- required_env("PCWDE_COMPLETE_LOCUS_TREES_OUTPUT_DIR")
QC_OUT <- file.path(FIG_OUT, "evidence")
if (!dir.exists(TREE_DIR)) stop("Figure S3 input directory does not exist: ", TREE_DIR)
if (dir.exists(FIG_OUT) && length(list.files(FIG_OUT, all.files=TRUE, no..=TRUE))) {
  stop("refusing non-empty output directory: ", FIG_OUT)
}
dir.create(FIG_OUT, showWarnings=FALSE, recursive=TRUE)
dir.create(QC_OUT, showWarnings=FALSE, recursive=TRUE)

# Use the exact accepted Figure 6 per-tip expression codes, including the
# tissue/life-stage prefix and the dot/+/++ magnitude tier.
RNA_OVERLAY_FILE <- file.path(TREE_DIR, "main_figure6_rna_evidence_overlay.tsv")
if (!file.exists(RNA_OVERLAY_FILE)) stop("missing accepted Figure 6 RNA overlay")
rna_overlay <- read.delim(RNA_OVERLAY_FILE, stringsAsFactors=FALSE, check.names=FALSE)
if (!all(c("family", "tip", "displayed_code") %in% names(rna_overlay)) ||
    anyDuplicated(paste(rna_overlay$family, rna_overlay$tip))) {
  stop("invalid accepted Figure 6 RNA overlay")
}

# Accepted display palette
palette_file <- file.path(TREE_DIR, "species_colors.R")
if (!file.exists(palette_file)) stop("missing accepted species palette")
source(palette_file)


# Chr display: strip chr prefix, 1-9 as-is, 10, blank for >10/unknown
chr_display <- function(x) {
  if (is.na(x) || x %in% c("", "outgroup", "unknown", "transcript", "fragment")) return("")
  x <- sub("^chr", "", x)
  n <- suppressWarnings(as.integer(x))
  if (is.na(n)) return("")
  if (n >= 1 && n <= 10) return(as.character(n))
  return("")  # chr11+ = unplaced RagTag scaffolds
}

# Tissue display with life stage prefix
tissue_display <- function(x) {
  if (is.na(x) || x %in% c("-", "")) return("")
  mapping <- c(
    "Larval"="L", "Venom"="AV", "Salivary"="LS", "Gall"="LG",
    "Gut"="LG", "Malpighian"="LM", "Head"="LH", "Pelt"="LP",
    "Egg"="", "WholeBody"="LW",
    "Other"="", "Mixed"="LX", "Carcass"="AC", "Accessory"="AA",
    # Legacy single-char codes
    "G"="LG", "S"="LS", "M"="LX", "V"="AV",
    "L"="L", "E"="", "x"="L", "C"="AC", "CR"="AC"
  )
  result <- mapping[x]
  if (is.na(result)) return("")
  return(result)
}

# ── Display label: clean tip names ──────────────────────────────
make_display_label <- function(label, pruned = FALSE) {
  if (is.na(label)) return(label)
  if (!pruned) return(label)
  # Clean TRINITY labels: BKS→BK, drop everything after TRINITY
  if (grepl("TRINITY", label)) {
    label <- sub("^BKS_", "BK_", label)
    label <- sub("_TRINITY.*", "_TRINITY", label)
    return(label)
  }
  # Strip cynipoid| or bacterium| prefix
  if (grepl("\\|", label)) {
    parts <- strsplit(label, "\\|")[[1]]
    if (parts[1] == "bacterium") {
      # bacterium|accession|species -> species
      label <- parts[length(parts)]
    } else {
      # cynipoid|Species_name_N -> Species_name_N
      label <- paste(parts[-1], collapse = "|")
    }
  }
  # Truncate genomic labels to binomial + copy number
  label <- sub("_subsp\\..*", "", label)
  label <- sub("_pv\\..*", "", label)
  parts <- strsplit(label, "_")[[1]]
  parts <- parts[parts != ""]
  if (length(parts) >= 2) {
    label <- paste(parts[1:2], collapse = "_")
  } else if (length(parts) >= 2) {
    label <- paste(parts[1:2], collapse = "_")
  }
  # Fix known truncated names
  label <- sub("spectabi$", "spectabilis", label)
  label
}

plotmath_quote <- function(x) encodeString(x, quote='"')
scientific_plotmath <- function(x) {
  if (is.na(x) || !nzchar(x)) return('""')
  if (grepl(" ex ", x, fixed=TRUE)) {
    parts <- strsplit(x, " ex ", fixed=TRUE)[[1]]
    return(paste0("italic(", plotmath_quote(parts[1]),
                  ")~' ex '~italic(", plotmath_quote(parts[2]), ")"))
  }
  tokens <- strsplit(x, " +")[[1]]
  if (length(tokens) >= 2 && tokens[2] %in% c("sp", "sp.")) {
    suffix <- if (length(tokens) > 2) paste(tokens[-c(1,2)], collapse=" ") else ""
    answer <- paste0("italic(", plotmath_quote(tokens[1]), ")~' sp.'")
    if (nzchar(suffix)) answer <- paste0(answer, "~", plotmath_quote(paste0(" ", suffix)))
    return(answer)
  }
  if (length(tokens) >= 2) {
    suffix <- if (length(tokens) > 2) paste(tokens[-c(1,2)], collapse=" ") else ""
    answer <- paste0("italic(", plotmath_quote(paste(tokens[1:2], collapse=" ")), ")")
    if (nzchar(suffix)) answer <- paste0(answer, "~", plotmath_quote(paste0(" ", suffix)))
    return(answer)
  }
  paste0("italic(", plotmath_quote(x), ")")
}
support_at_least_70 <- function(x) {
  numeric_value <- suppressWarnings(as.numeric(as.character(x)))
  ifelse(!is.na(numeric_value) & numeric_value >= 70, sprintf("%.0f", numeric_value), "")
}

# ── Remove redundant TRINITY tips ───────────────────────────────
# Keep one representative per unique tissue type per gene cluster
find_redundant_tips <- function(meta_df) {
  trinity <- meta_df[grepl("TRINITY", meta_df$label), ]
  if (nrow(trinity) == 0) return(character(0))
  # Cluster = everything before _iN.p1
  trinity$cluster <- sub("_i[0-9]+\\.p1$", "", trinity$label)
  trinity$tis_disp <- sapply(trinity$Tissue, tissue_display)

  to_drop <- character(0)
  for (clust in unique(trinity$cluster)) {
    grp <- trinity[trinity$cluster == clust, ]
    if (nrow(grp) <= 1) next
    keep <- character(0)
    for (tis in unique(grp$tis_disp)) {
      matches <- grp[grp$tis_disp == tis, ]
      keep <- c(keep, matches$label[1])
    }
    to_drop <- c(to_drop, setdiff(grp$label, keep))
  }
  to_drop
}


# Safely root on a requested outgroup, falling back gracefully:
#   1. Filter requested tips to those actually in the tree.
#   2. Try full multi-tip root with resolve.root=TRUE.
#   3. If that fails (paraphyletic, etc.), try each single tip in order
#      until one works.
#   4. If all fail, keep tree as-is and warn.
safe_root <- function(tr, outgroup, label="tree") {
  og <- outgroup[outgroup %in% tr$tip.label]
  if (length(og) == 0) {
    cat("  [safe_root]", label, "- no requested outgroup tips present; leaving as-is\n")
    return(tr)
  }
  # Unroot first to avoid degenerate re-rooting of already-rooted trees
  # (re-rooting a rooted NHX tree creates a structure that hangs ggtree)
  if (is.rooted(tr)) {
    original_tips <- tr$tip.label
    unrooted <- try(unroot(tr), silent = TRUE)
    if (inherits(unrooted, "phylo")) {
      tr <- unrooted
      # unroot() on NHX trees can convert a labeled root node into a phantom tip
      phantom <- setdiff(tr$tip.label, original_tips)
      if (length(phantom) > 0) {
        cat("  [safe_root]", label, "- dropping", length(phantom), "phantom tip(s) from unroot\n")
        tr <- drop.tip(tr, phantom)
      }
    } else {
      cat("  [safe_root]", label, "- unroot failed; trying root from original tree\n")
    }
  }
  if (length(og) > 1) {
    rooted <- try(root(tr, outgroup = og, resolve.root = TRUE), silent = TRUE)
    if (inherits(rooted, "phylo")) return(rooted)
    cat("  [safe_root]", label, "- multi-tip outgroup not monophyletic; falling back to single tip\n")
  }
  for (tip in og) {
    rooted <- try(root(tr, outgroup = tip, resolve.root = TRUE), silent = TRUE)
    if (inherits(rooted, "phylo")) {
      cat("  [safe_root]", label, "- rooted on", tip, "\n")
      return(rooted)
    }
  }
  cat("  [safe_root]", label, "- ALL rooting attempts failed; leaving as-is\n")
  return(tr)
}

descendant_tips <- function(tr, node) {
  tips <- integer(0)
  frontier <- node
  while (length(frontier)) {
    children <- tr$edge[tr$edge[,1] %in% frontier, 2]
    tips <- c(tips, children[children <= length(tr$tip.label)])
    frontier <- children[children > length(tr$tip.label)]
  }
  unique(tips)
}

edge_signature <- function(tr) {
  lengths <- if (is.null(tr$edge.length)) rep(NA_real_, nrow(tr$edge)) else tr$edge.length
  sort(sprintf("%d>%d:%s", tr$edge[,1], tr$edge[,2],
               ifelse(is.na(lengths), "NA", sprintf("%.17g", lengths))))
}

orient_phylo_tip_first <- function(tr, target) {
  if (!inherits(tr, "phylo")) stop("orientation input is not phylo")
  target_index <- match(target, tr$tip.label)
  if (is.na(target_index)) stop(paste("orientation target missing", target))
  before_edges <- edge_signature(tr)
  before_tips <- sort(tr$tip.label)
  before_labels <- tr$node.label
  root_node <- length(tr$tip.label) + 1L
  target_path <- integer(0)
  node <- tr$edge[match(target_index, tr$edge[,2]), 1]
  while (!is.na(node) && node != root_node) {
    target_path <- c(target_path, node)
    node <- tr$edge[match(node, tr$edge[,2]), 1]
  }
  edge_order <- integer(0)
  visit <- function(node) {
    rows <- which(tr$edge[,1] == node)
    if (node %in% target_path && length(rows) > 1L) {
      child_has_target <- vapply(tr$edge[rows,2], function(child) {
        if (child <= length(tr$tip.label)) child == target_index
        else target_index %in% descendant_tips(tr, child)
      }, logical(1))
      rows <- c(rows[child_has_target], rows[!child_has_target])
    }
    for (row in rows) {
      edge_order <<- c(edge_order, row)
      child <- tr$edge[row,2]
      if (child > length(tr$tip.label)) visit(child)
    }
  }
  visit(root_node)
  if (length(edge_order) != nrow(tr$edge) || length(unique(edge_order)) != nrow(tr$edge)) {
    stop("orientation did not produce a complete unique edge order")
  }
  tr$edge <- tr$edge[edge_order, , drop=FALSE]
  if (!is.null(tr$edge.length)) tr$edge.length <- tr$edge.length[edge_order]
  list(tree=tr,
       topology_invariant=identical(before_edges, edge_signature(tr)),
       tips_invariant=identical(before_tips, sort(tr$tip.label)),
       support_invariant=identical(before_labels, tr$node.label))
}

# Select the largest exclusively Diplolepididae clade. If several clades have
# the same number of tips, use the longest stem branch as a display-only
# tie-break; this is deterministic and does not change the inferred splits.
select_diplolepididae_clade <- function(tr) {
  candidates <- list()
  for (node in seq.int(length(tr$tip.label) + 1L,
                       length(tr$tip.label) + tr$Nnode)) {
    idx <- descendant_tips(tr, node)
    labs <- tr$tip.label[idx]
    if (length(labs) >= 2L && all(grepl("^cynipoid\\|Diplolepis_", labs))) {
      edge_row <- match(node, tr$edge[, 2])
      stem <- if (!is.na(edge_row) && !is.null(tr$edge.length)) tr$edge.length[edge_row] else 0
      support <- if (!is.null(tr$node.label)) tr$node.label[node - length(tr$tip.label)] else NA_character_
      candidates[[length(candidates) + 1L]] <- list(
        node=node, tips=sort(labs), size=length(labs), stem=stem, support=support
      )
    }
  }
  if (!length(candidates)) stop("no exclusively Diplolepididae clade with at least two tips")
  max_size <- max(vapply(candidates, function(x) x$size, integer(1)))
  largest <- candidates[vapply(candidates, function(x) x$size == max_size, logical(1))]
  stem_lengths <- vapply(largest, function(x) x$stem, numeric(1))
  winner <- largest[[which.max(stem_lengths)]]
  winner$size_ties <- length(largest)
  winner
}

# Accepted complete Figure S3 trees
full_trees <- list(
  list(prefix="complete_locus_trees_GH5_2", file="complete_locus_GH5_2.rooted.treefile",
       family="GH5_2", w=28, h=60,
       tip_size=6, node_size=5, branch_size=0.4,
       annot_size=5, png_dpi=150,
       annot_offset=12, sp_w=0.6, chr_w=0.9, tis_w=2.3,
       pruned=FALSE),
  list(prefix="complete_locus_trees_PL1_8", file="complete_locus_PL1_8.rooted.treefile",
       family="PL1_8", w=28, h=70,
       tip_size=6, node_size=5, branch_size=0.4,
       annot_size=5, png_dpi=150,
       annot_offset=8, sp_w=0.6, chr_w=0.9, tis_w=2.3,
       pruned=FALSE),
  list(prefix="complete_locus_trees_PL4_4", file="complete_locus_PL4_4.rooted.treefile",
       family="PL4_4", w=28, h=85,
       tip_size=6, node_size=5, branch_size=0.4,
       annot_size=5, png_dpi=150,
       annot_offset=12, sp_w=0.6, chr_w=0.9, tis_w=2.3,
       pruned=FALSE)
)

trees <- full_trees

# Cell half-height (full = 1.0, matches tip spacing)
cell_h <- 0.5

# Store pruned plots for combined figure
pruned_plots <- list()

for (t in trees) {
  cat("Rendering", t$prefix, "\n")

  # ── Read NHX, extract metadata, root ──
  nhx <- read.nhx(file.path(TREE_DIR, t$file))
  meta <- as_tibble(nhx)

  # Extract metadata for tips
  meta_df <- data.frame(
    label   = meta$label,
    Sp      = as.character(meta$Sp),
    Chr     = as.character(meta$Chr),
    Tissue  = as.character(meta$Tissue),
    stringsAsFactors = FALSE
  )
  meta_df <- meta_df[!is.na(meta_df$label), ]

  # Display-only transformation: remove the single bacterial tip, root on the
  # largest exclusively Diplolepididae clade, and ladderize. The accepted
  # inferred cynipoid topology, branch lengths, and support values are retained.
  tr <- as.phylo(nhx)
  original_tips <- tr$tip.label
  bacterial_tips <- grep("^bacterium\\|", tr$tip.label, value=TRUE)
  if (length(bacterial_tips) != 1L) {
    stop(paste("expected exactly one bacterial display root in", t$family,
               "but found", length(bacterial_tips)))
  }
  tr <- drop.tip(tr, bacterial_tips)
  meta_df <- meta_df[meta_df$label %in% tr$tip.label, , drop=FALSE]
  ov <- rna_overlay[rna_overlay$family == t$family, , drop=FALSE]
  overlay_index <- match(meta_df$label, ov$tip)
  if (anyNA(overlay_index)) {
    stop(paste("Figure 6 expression overlay does not cover all", t$family, "tips"))
  }
  meta_df$Tissue <- ov$displayed_code[overlay_index]
  cat("  Figure 6 high-resolution expression overlay:", nrow(meta_df), "tips matched\n")
  selected <- select_diplolepididae_clade(tr)
  tr <- root(unroot(tr), outgroup=selected$tips, resolve.root=TRUE, edgelabel=TRUE)

  # Drop redundant TRINITY tips (pruned trees only)
  is_pruned <- isTRUE(t$pruned)
  if (is_pruned) {
    redundant <- find_redundant_tips(meta_df)
    if (length(redundant) > 0) {
      cat("  Dropping redundant:", paste(redundant, collapse=", "), "\n")
      tr <- drop.tip(tr, redundant)
      meta_df <- meta_df[!meta_df$label %in% redundant, ]
    }
  }

  tr <- ladderize(tr, right = FALSE)

  # Build ggtree cladogram and join metadata
  p <- ggtree(tr, branch.length = "none", ladderize = FALSE, size = t$branch_size) %<+% meta_df

  tip_data <- p$data[p$data$isTip, ]
  max_tip_x <- max(tip_data$x)
  n_tips <- nrow(tip_data)
  metadata_tip_rows <- length(unique(meta_df$label[meta_df$label %in% tr$tip.label]))
  root_children <- tr$edge[tr$edge[,1] == length(tr$tip.label) + 1L, 2]
  root_descendants <- lapply(root_children, function(node) {
    if (node <= length(tr$tip.label)) tr$tip.label[node]
    else tr$tip.label[descendant_tips(tr, node)]
  })
  root_exact <- any(vapply(root_descendants, function(x) setequal(x, selected$tips), logical(1)))
  tips_exact <- setequal(tr$tip.label, setdiff(original_tips, bacterial_tips))
  no_bacteria <- !any(grepl("^bacterium\\|", tr$tip.label))
  gate_status <- if (root_exact && tips_exact && no_bacteria && n_tips == metadata_tip_rows) "PASS" else "FAIL"
  gate <- data.frame(
    family=t$family, input_tips=length(original_tips), displayed_tips=n_tips,
    removed_bacterial_tip=paste(bacterial_tips, collapse=";"),
    diplolepididae_clade_size=selected$size,
    largest_clade_ties=selected$size_ties,
    tie_break="longest_stem_branch",
    selected_stem_length=selected$stem,
    selected_support=selected$support,
    selected_tips=paste(selected$tips, collapse=";"),
    expression_nonblank=sum(!is.na(meta_df$Tissue) & nzchar(meta_df$Tissue)),
    expression_dot=sum(meta_df$Tissue == "·", na.rm=TRUE),
    expression_plus=sum(grepl("[+]$", meta_df$Tissue), na.rm=TRUE),
    expression_plusplus=sum(grepl("[+][+]$", meta_df$Tissue), na.rm=TRUE),
    root_exact=root_exact, tips_exact=tips_exact, no_bacteria=no_bacteria,
    metadata_tip_rows=metadata_tip_rows, status=gate_status,
    stringsAsFactors=FALSE
  )
  write.table(gate, file=file.path(QC_OUT, paste0(t$family, "_display_root_gate.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE)
  write.tree(tr, file=file.path(QC_OUT, paste0("complete_locus_", t$family,
                                               ".diplolepididae_rooted_ladderized.treefile")))
  if (gate_status != "PASS") stop(paste("display-root gate failed", t$family))

  # ── Column x-ranges (flush, no gaps) — per-tree sizing ──
  sp_xmin  <- max_tip_x + t$annot_offset
  sp_xmax  <- sp_xmin + t$sp_w
  chr_xmin <- sp_xmax
  chr_xmax <- chr_xmin + t$chr_w
  tis_xmin <- chr_xmax
  tis_xmax <- tis_xmin + t$tis_w

  # Column centers for text
  chr_cx <- (chr_xmin + chr_xmax) / 2
  tis_cx <- (tis_xmin + tis_xmax) / 2

  # ── Build annotation dataframes (vectorized — avoids O(n²) layer accumulation) ──
  # Species color rects
  sp_fills <- sapply(as.character(tip_data$Sp), function(sp) {
    if (!is.na(sp) && sp %in% names(species_colors)) species_colors[sp] else "white"
  })
  rect_sp <- data.frame(
    xmin = sp_xmin, xmax = sp_xmax,
    ymin = tip_data$y - cell_h, ymax = tip_data$y + cell_h,
    fill = sp_fills, stringsAsFactors = FALSE
  )

  # Chr rects + text
  chr_vals <- sapply(locus_tip_chromosomes(t$family, as.character(tip_data$label)), chr_display)
  marks <- locus_tip_marks(t$family, as.character(tip_data$label))
  chr_vals[nzchar(marks)] <- marks[nzchar(marks)]
  write_locus_chr_labels(QC_OUT, t$family, "complete", as.character(tip_data$label), chr_vals)
  rect_chr <- data.frame(
    xmin = chr_xmin, xmax = chr_xmax,
    ymin = tip_data$y - cell_h, ymax = tip_data$y + cell_h,
    stringsAsFactors = FALSE
  )
  chr_text <- data.frame(
    x = chr_cx, y = tip_data$y, label = chr_vals,
    colour = ifelse(chr_vals %in% c("T", "U"), "grey58", "black"),
    stringsAsFactors = FALSE
  )
  chr_text <- chr_text[nchar(chr_text$label) > 0, , drop = FALSE]

  # Tissue rects + text
  tis_vals <- as.character(tip_data$Tissue)
  tis_vals[is.na(tis_vals)] <- ""
  rect_tis <- data.frame(
    xmin = tis_xmin, xmax = tis_xmax,
    ymin = tip_data$y - cell_h, ymax = tip_data$y + cell_h,
    stringsAsFactors = FALSE
  )
  tis_text <- data.frame(
    x = tis_cx, y = tip_data$y, label = tis_vals,
    stringsAsFactors = FALSE
  )
  tis_text <- tis_text[nchar(tis_text$label) > 0, , drop = FALSE]

  # Add all annotation layers in 5 calls (not n_tips × 5)
  p <- p +
    geom_rect(data = rect_sp,  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill = rect_sp$fill, color = "grey70", linewidth = 0.15, inherit.aes = FALSE) +
    geom_rect(data = rect_chr, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill = "white", color = "grey70", linewidth = 0.15, inherit.aes = FALSE) +
    geom_rect(data = rect_tis, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill = "white", color = "grey70", linewidth = 0.15, inherit.aes = FALSE)
  if (nrow(chr_text) > 0) {
    p <- p + geom_text(data = chr_text, aes(x=x, y=y, label=label),
                       color=chr_text$colour, size = t$annot_size, fontface = "bold", inherit.aes = FALSE)
  }
  if (nrow(tis_text) > 0) {
    p <- p + geom_text(data = tis_text, aes(x=x, y=y, label=label),
                       size = t$annot_size, fontface = "plain", inherit.aes = FALSE)
  }

  # Create display labels from NHX Sp attribute (correct species names)
  # Display-only species name renames (DB-canonical; underlying NHX Sp
  # attribute and tree-tip labels keep the legacy underscore form).
  apply_display_renames <- function(x) {
    if (is.na(x)) return(x)
    x <- gsub("Isocolus centaureae", "Eubothrus centaureae", x, fixed = TRUE)
    x <- gsub("Ceroptres masudai", "Ceroptres kovalevi", x, fixed = TRUE)
    x <- gsub("Bradysia coprophila", "Bradysia tilicola", x, fixed = TRUE)
    x <- gsub("Eschatocerus acacia", "Eschatocerus acaciae", x, fixed = TRUE)
    if (grepl("Synergus.*D[.] quercuslanigerum.*inquiline", x) ||
        grepl("Synergus.*quercuslanigerum.*inquiline", x)) {
      return("Synergus ex Druon quercuslanigerum")
    }
    if (grepl("Synergus.*DI.*inquiline", x) ||
        grepl("Synergus.*D[.] ignotum.*inquiline", x) ||
        grepl("Synergus.*Druon ignotum.*inquiline", x)) {
      return("Synergus ex Druon ignotum")
    }
    x
  }
  p$data$display_label <- sapply(seq_len(nrow(p$data)), function(i) {
    row <- p$data[i, ]
    if (!is.na(row$Sp) && row$Sp != "") {
      apply_display_renames(gsub("_", " ", row$Sp))
    } else if (!is.na(row$label)) {
      apply_display_renames(gsub("_", " ", make_display_label(row$label, pruned = is_pruned)))
    } else {
      ""
    }
  })
  p$data$display_label_math <- vapply(
    p$data$display_label, scientific_plotmath, character(1)
  )
  header_text <- data.frame(
    x=c((sp_xmin + sp_xmax) / 2,
        (chr_xmin + chr_xmax) / 2,
        (tis_xmin + tis_xmax) / 2),
    y=rep(0, 3),
    label=c("Species", "Chr.", "Expression")
  )

  p <- p +
    scale_y_reverse() +
    geom_text(data=header_text, aes(x=x, y=y, label=label),
              size=3.2, fontface="bold", inherit.aes=FALSE) +
    geom_tiplab(aes(label = display_label_math), size = t$tip_size,
                offset = 0.3, parse=TRUE) +
    geom_nodelab(aes(label = support_at_least_70(label)),
                 size = t$node_size, hjust = 1.3, vjust = 1.5,
                 color = "grey30") +
    xlim_tree(tis_xmax + 2)

  # Store pruned plots for combined figure
  if (is_pruned) {
    pruned_plots[[t$family]] <- p
  }

  outpdf <- file.path(FIG_OUT, paste0(t$prefix, "_ggtree.pdf"))
  outsvg <- file.path(FIG_OUT, paste0(t$prefix, "_ggtree.svg"))
  outpng <- file.path(FIG_OUT, paste0(t$prefix, "_ggtree.png"))

  ggsave(outpdf, p, width = t$w, height = t$h, limitsize = FALSE)
  ggsave(outsvg, p, width = t$w, height = t$h, limitsize = FALSE)
  ggsave(outpng, p, width = t$w, height = t$h, dpi = t$png_dpi,
         bg = "white", limitsize = FALSE)
  cat("  ->", outpdf, outsvg, outpng, "\n")
}
cat("Done!\n")

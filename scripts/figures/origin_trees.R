library(ggtree)
library(treeio)
library(ggplot2)
library(ape)

required_env <- function(name) {
  value <- Sys.getenv(name, unset="")
  if (!nzchar(value)) stop(name, " is required")
  value
}

# Preserve the accepted Figure 3 bytes. R's PDF device emits the render time in
# two fixed-width metadata fields even when all graphical content is unchanged.
canonicalize_figure3_pdf_timestamp <- function(path) {
  stamp <- "20260825211502"
  bytes <- readBin(path, what="raw", n=file.info(path)$size)
  replace_one <- function(payload, label) {
    prefix <- charToRaw(paste0("/", label, " (D:"))
    starts <- which(vapply(
      seq_len(length(payload) - length(prefix) + 1L),
      function(start) identical(payload[start:(start + length(prefix) - 1L)], prefix),
      logical(1)
    ))
    if (length(starts) != 1L) stop("expected one PDF ", label, " field in ", path)
    date_start <- starts[[1]] + length(prefix)
    date_end <- date_start + 13L
    observed <- payload[date_start:date_end]
    if (!identical(payload[date_end + 1L], charToRaw(")")) ||
        any(!observed %in% charToRaw("0123456789"))) {
      stop("malformed PDF ", label, " field in ", path)
    }
    payload[date_start:date_end] <- charToRaw(stamp)
    payload
  }
  bytes <- replace_one(bytes, "CreationDate")
  bytes <- replace_one(bytes, "ModDate")
  writeBin(bytes, path)
}
INPUT_DIR <- required_env("PCWDE_ORIGIN_TREES_INPUT_DIR")
FIG_OUT <- required_env("PCWDE_ORIGIN_TREES_OUTPUT_DIR")
# Display trees and key nodes are build intermediates: true prunings of the
# complete origin trees (scripts/trees/prune_origin_display_trees.R).
TREE_DIR <- required_env("PCWDE_ORIGIN_DISPLAY_TREE_DIR")
KEY_NODES <- read.delim(file.path(TREE_DIR, "origin_key_nodes.tsv"), sep="\t", quote="",
                        colClasses="character", na.strings=character(0))
DOMAIN_DIR <- INPUT_DIR
EVIDENCE_DIR <- file.path(FIG_OUT, "evidence")
if (dir.exists(FIG_OUT) && length(list.files(FIG_OUT, all.files=TRUE, no..=TRUE))) {
  stop("refusing non-empty output directory: ", FIG_OUT)
}
dir.create(FIG_OUT, showWarnings=FALSE, recursive=TRUE)
dir.create(EVIDENCE_DIR, showWarnings=FALSE, recursive=TRUE)
RENDERER_DIR <- INPUT_DIR
source(file.path(RENDERER_DIR, "renderer_native_order.R"))
# Shared key layout (scripts/figures/key_layout.R), located beside this renderer.
.renderer_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
if (length(.renderer_file) != 1) stop("cannot locate the running renderer script")
source(file.path(dirname(normalizePath(.renderer_file)), "key_layout.R"))

# ── Eukaryote category detection ──────────────────────────────
euk_prefixes <- c("cynipoid", "midge", "bee", "beetle",
                  "nematode", "oomycete", "fungus", "plant", "protist")

get_category <- function(label) {
  prefix <- sub("\\|.*", "", label)
  if (prefix %in% euk_prefixes) return(prefix)
  return(NA_character_)
}

# Strip type prefix and optionally accession from tip labels
make_display_label <- function(label, remove_accession = FALSE) {
  if (is.na(label) || !grepl("\\|", label)) return(label)
  # Remove type prefix (everything before first |)
  cleaned <- sub("^[^|]*\\|", "", label)
  if (remove_accession) {
    if (grepl("\\|", cleaned)) {
      # Format: accession|species - take species part
      cleaned <- sub("^[^|]*\\|", "", cleaned)
    }
    # Inquiline-naming-convention shortcut (catch BEFORE subsp truncation,
    # since the keep-set labels still carry the messy underscore-padded form)
    if (grepl("Synergus__BK7_inquiline", cleaned)) return("Synergus ex Belonocnema kinseyi")
    if (grepl("Ceroptres.*DI.*inquiline", cleaned)) return("Ceroptres ex Druon ignotum")
    if (grepl("Synergus.*quercuslanigerum.*inquiline|Synergus.*DL.*inquiline|Synergus__DI_inquiline", cleaned)) return("Synergus ex Druon quercuslanigerum")
    if (grepl("Synergus.*DC.*inquiline|Synergus__DC", cleaned)) return("Synergus ex Druon")
    # Truncate to binomial: strip subsp/pv/var suffixes, then take Genus_species
    cleaned <- sub("_subsp\\..*", "", cleaned)
    cleaned <- sub("_pv\\..*", "", cleaned)
    cleaned <- sub("_var\\..*", "", cleaned)
    parts <- strsplit(cleaned, "_")[[1]]
    parts <- parts[parts != ""]
    if (length(parts) >= 2) {
      cleaned <- paste(parts[1:2], collapse = "_")
    }
    # Strip underscores from the final label
    cleaned <- gsub("_", " ", cleaned)
    # Display-only species name renames (DB-canonical; underlying tree-tip
    # labels keep the legacy underscore form for lookup).
    cleaned <- gsub("Eschatocerus acacia", "Eschatocerus acaciae", cleaned, fixed = TRUE)
    cleaned <- gsub("Isocolus centaureae", "Eubothrus centaureae", cleaned, fixed = TRUE)
    cleaned <- gsub("Ceroptres masudai", "Ceroptres kovalevi", cleaned, fixed = TRUE)
    cleaned <- gsub("Bradysia coprophila", "Bradysia tilicola", cleaned, fixed = TRUE)
  }
  cleaned
}

# Colors for clade shading
clade_colors <- c(
  cynipoid = "#2171b5",
  midge    = "#8c6d31",
  bee      = "#e7ba52",
  beetle   = "#607D8B",
  nematode = "#d62728",
  oomycete = "#637939",
  fungus   = "#9467bd",
  plant    = "#2ca02c",
  protist  = "#17becf"
)

# ── Domain info ───────────────────────────────────────────────
domain_colors <- c(
  A = "#1B7837",   # GH5_2 - Emerald (paper anchor)
  B = "#00A6D6",   # CBM2 - Cyan blue
  C = "#A6D854",   # CBM3 - Chartreuse
  D = "#E7298A",   # CBM5 - Magenta
  E = "#984EA3",   # CBM6 - Purple
  F = "#FF7F00",   # CBM10 - Orange
  G = "#8DD3C7",   # CBM12 - Pale Teal
  H = "#252525",   # CBM17 - Near Black
  I = "#BEBADA",   # CBM28 - Lavender
  J = "#80B1D3",   # CBM73 - Sky Blue
  K = "#542788",   # CBM81 - Deep Violet
  L = "#2166AC",   # PL1_8 - Sapphire (paper anchor)
  M = "#D98C00",   # CBM13 / B_lectin - Amber (4th HGT family anchor)
  N = "#D6604D",   # PL4 - Brick (near garnet)
  O = "#C2185B",   # PL4_2 - Raspberry (near garnet)
  P = "#B2182B",   # PL4_4 - Garnet (paper anchor)
  Q = "#E31A1C",   # PL4_5 - Scarlet (near garnet)
  R = "#F781BF",   # CBM22 - Pink
  S = "#FFD92F",   # CBM32 - Yellow
  T = "#A6761D",   # CBM64 - Sienna
  U = "#7570B3",   # Glyco_hydro_30 - Slate (GH30 catalytic)
  V = "#BEAED4",   # Glyco_hydro_30C - Pale Slate (C-terminal β-sandwich)
  W = "#66C2A5",   # CBM_4_9 - Teal
  X = "#FDB462",   # fn3 - Pale Orange
  Y = "#7FC97F",   # Ricin_B_lectin (plant) - Pale Green
  Z = "#BDBDBD",   # Big_5 - Grey
  u = "#5A4FCF",   # Glyco_hydro_59 - Deep Indigo (GH59 variant)
  v = "#9897E0",   # Glyco_hydr_30_2 - Periwinkle (GH30 variant 2)
  w = "#40A28C",   # CBM_6 - Dark Teal
  y = "#B8D88C",   # RicinB_lectin_2 - Pale Lime
  c = "#006D2C",   # CBM_35_2 - Forest
  s = "#525252",   # SLH - Dark Grey
  p = "#E41A1C",   # PI-PLC-X - Crimson (bacterial phospholipase C donor marker)
  # Dual-shade CBM13 classing (from build_domain_faa.py DUAL_SHADE_FAMILIES):
  m = "#8B5A00",   # B_lectin (degraded) — dark amber, PF01453 c-E 1e-3..1e-1
  q = "#7A0D16"    # Region A only — dark crimson, custom HMM donor-region residual
)

domain_names <- c(
  A = "GH5_2",  B = "CBM2",  C = "CBM3",  D = "CBM5",
  E = "CBM6",   F = "CBM10", G = "CBM12", H = "CBM17",
  I = "CBM28",  J = "CBM73", K = "CBM81", L = "PL1_8",
  M = "CBM13",  N = "PL4",   O = "PL4_2", P = "PL4_4",
  Q = "PL4_5",  R = "CBM22", S = "CBM32", T = "CBM64",
  U = "GH30",   V = "GH30_C",W = "CBM4_9",X = "fn3",
  Y = "Ricin_B",Z = "Big_5",
  u = "GH59",   v = "GH30_2",w = "CBM6",  y = "RicinB_2",
  c = "CBM35_2",s = "SLH",   p = "PI-PLC-X",
  m = "CBM13 (weak)", q = "Region A only"
)

# ── Parse domain FAA file ────────────────────────────────────
parse_domains <- function(faa_file) {
  lines <- readLines(faa_file)
  result <- list()
  cur_name <- NULL
  cur_seq <- ""

  process_seq <- function(name, seq) {
    if (is.null(name) || nchar(seq) == 0) return(NULL)
    chars <- strsplit(seq, "")[[1]]
    blocks <- list()
    cur_ch <- NULL
    cur_start <- NULL

    plen <- length(chars)
    for (i in seq_along(chars)) {
      ch <- chars[i]
      if (ch != "-") {
        if (is.null(cur_ch) || ch != cur_ch) {
          if (!is.null(cur_ch)) {
            blocks[[length(blocks) + 1]] <- data.frame(
              label = name, start = cur_start, end = i - 1,
              domain = cur_ch, prot_len = plen, stringsAsFactors = FALSE
            )
          }
          cur_ch <- ch
          cur_start <- i
        }
      } else {
        if (!is.null(cur_ch)) {
          blocks[[length(blocks) + 1]] <- data.frame(
            label = name, start = cur_start, end = i - 1,
            domain = cur_ch, prot_len = plen, stringsAsFactors = FALSE
          )
          cur_ch <- NULL
        }
      }
    }
    if (!is.null(cur_ch)) {
      blocks[[length(blocks) + 1]] <- data.frame(
        label = name, start = cur_start, end = plen,
        domain = cur_ch, prot_len = plen, stringsAsFactors = FALSE
      )
    }
    if (length(blocks) > 0) do.call(rbind, blocks) else NULL
  }

  for (line in lines) {
    line <- trimws(line)
    if (startsWith(line, ">")) {
      df <- process_seq(cur_name, cur_seq)
      if (!is.null(df)) result[[length(result) + 1]] <- df
      cur_name <- sub("^>", "", line)
      cur_seq <- ""
    } else if (nchar(line) > 0) {
      cur_seq <- paste0(cur_seq, line)
    }
  }
  df <- process_seq(cur_name, cur_seq)
  if (!is.null(df)) result[[length(result) + 1]] <- df

  do.call(rbind, result)
}

# ── Build shading rectangles from tip y-positions ─────────────
build_shade_rects <- function(p, strip_width=3) {
  d <- p$data
  tip_data <- d[d$isTip, ]
  max_tip_x <- max(tip_data$x)
  xmin_rect <- max_tip_x
  xmax_rect <- max_tip_x + strip_width
  tip_data <- tip_data[order(tip_data$y), ]
  tip_data$cat <- sapply(tip_data$label, get_category)

  pad <- 0.49

  rects <- list()
  cur_cat <- NA
  cur_ymin <- NA
  cur_ymax <- NA

  for (i in seq_len(nrow(tip_data))) {
    cat_i <- tip_data$cat[i]
    yi <- tip_data$y[i]

    if (!is.na(cat_i) && identical(cat_i, cur_cat)) {
      cur_ymax <- yi + pad
    } else {
      if (!is.na(cur_cat)) {
        rects[[length(rects) + 1]] <- data.frame(
          xmin = xmin_rect, xmax = xmax_rect,
          ymin = cur_ymin, ymax = cur_ymax,
          cat = cur_cat, stringsAsFactors = FALSE
        )
      }
      if (!is.na(cat_i)) {
        cur_cat <- cat_i
        cur_ymin <- yi - pad
        cur_ymax <- yi + pad
      } else {
        cur_cat <- NA
      }
    }
  }
  if (!is.na(cur_cat)) {
    rects[[length(rects) + 1]] <- data.frame(
      xmin = xmin_rect, xmax = xmax_rect,
      ymin = cur_ymin, ymax = cur_ymax,
      cat = cur_cat, stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rects)
}

# ── Tree configs ──────────────────────────────────────────────
# Pruned trees (subset)
pruned_trees <- list(
  list(prefix="origin_trees_GH5_2", file="origin_display_GH5_2_fig2.treefile",
       family="GH5_2", scope="figure3", key_view="fig2",
       domains="GH5_2_subset_domains.faa", w=20, h=14,
       tip_size=2.3, node_size=2.3, legend_size=2.3, branch_size=0.8,
       dom_height=0.35, png_dpi=600, pruned=TRUE, dominant="A"),
  list(prefix="origin_trees_PL1_8", file="origin_display_PL1_8_fig2.treefile",
       family="PL1_8", scope="figure3", key_view="fig2",
       domains="PL1_8_subset_domains.faa", w=20, h=14,
       tip_size=2.3, node_size=2.3, legend_size=2.3, branch_size=0.8,
       dom_height=0.35, png_dpi=600, pruned=TRUE, dominant="L"),
  list(prefix="origin_trees_PL4_4", file="origin_display_PL4_4_fig2.treefile",
       family="PL4_4", scope="figure3", key_view="fig2",
       domains="PL4_4_subset_domains.faa", w=20, h=14,
       tip_size=2.3, node_size=2.3, legend_size=2.3, branch_size=0.8,
       dom_height=0.35, png_dpi=600, pruned=TRUE, dominant="P")
)

# Fixed Figure S1 trees: 4/4/3 focal cynipoids plus all 166/63/79 accepted outgroups.
full_trees <- list(
  list(prefix="origin_trees_full_GH5_2", file="origin_display_GH5_2_s1.treefile",
       family="GH5_2", scope="s1", key_view="s1",
       domains="CynOGv3_000234_domains.faa", w=28, h=44,
       tip_size=5, node_size=4, legend_size=6, branch_size=0.35,
       dom_height=0.25, png_dpi=200, pruned=TRUE, dominant="A"),
  list(prefix="origin_trees_full_PL1_8", file="origin_display_PL1_8_s1.treefile",
       family="PL1_8", scope="s1", key_view="s1",
       domains="CynOGv3_000148_domains.faa", w=28, h=26,
       tip_size=5, node_size=4, legend_size=6, branch_size=0.35,
       dom_height=0.25, png_dpi=200, pruned=TRUE, dominant="L"),
  list(prefix="origin_trees_full_PL4_4", file="origin_display_PL4_4_s1.treefile",
       family="PL4_4", scope="s1", key_view="s1",
       domains="CynOGv3_000168_domains.faa", w=28, h=30,
       tip_size=5, node_size=4, legend_size=6, branch_size=0.35,
       dom_height=0.25, png_dpi=200, pruned=TRUE, dominant="P")
)

trees <- c(pruned_trees, full_trees)

pruned_plots <- list()   # store pruned plots for combined figure
rendered_order_rows <- list()
ladderization_rows <- list()
s1_sister_rows <- list()
drawn_key_rows <- list()
drawn_label_rows <- list()

strip_width     <- 5      # narrow clade-color band between tips and domain strip
domain_offset   <- 18    # legacy: floor for unpruned trees only
domain_width    <- 12    # legacy: only used for full (unpruned) trees
aa_per_unit     <- 50    # 1 x-unit of domain strip == 50 aa of protein (pruned trees)
label_gap       <- 4     # gap between label area and domain strip
units_per_char  <- 0.52  # accepted label-width estimate at tip_size 8
label_left_pad  <- 1     # compact gap between tree and tip labels
label_right_pad <- 1.5   # accepted label-to-domain clearance

# Figure 2 support labels: the cynipoid clade and the node joining it to its
# sister lineage are always labelled, whatever the value; every other internal
# node is labelled at >=70 except two-tip cherries. The two key-node values must
# equal BUILD/trees/origin_key_nodes.tsv; the drawn labels are written as evidence.
figure2_support_labels <- function(p, tr, t) {
  n <- Ntip(tr)
  cyn <- tr$tip.label[startsWith(tr$tip.label, "cynipoid|")]
  cnode <- ape::getMRCA(tr, cyn)
  jnode <- tr$edge[tr$edge[, 2] == cnode, 1]
  if (length(jnode) != 1) stop(t$family, ": cynipoid clade has no joining node")
  kids <- split(tr$edge[, 2], tr$edge[, 1])
  d <- p$data[!p$data$isTip, c("node", "x", "y")]
  d$support <- vapply(d$node, function(k) {
    v <- tr$node.label[k - n]; if (is.na(v)) "" else v }, character(1))
  d$role <- ifelse(d$node == cnode, "cynipoid_clade", ifelse(d$node == jnode, "join", ""))
  cherry <- vapply(d$node, function(k) all(kids[[as.character(k)]] <= n), logical(1))
  d$lab <- ifelse(nzchar(d$role), d$support,
                  ifelse(cherry, "", support_at_least_70(d$support)))
  for (role in c("cynipoid_clade", "join")) {
    want <- KEY_NODES$support[KEY_NODES$family == t$family & KEY_NODES$view == t$key_view &
                              KEY_NODES$node == role]
    got <- d$lab[d$role == role]
    if (length(want) != 1 || length(got) != 1 || !nzchar(got) || !identical(want, got)) {
      stop(t$family, " ", role, ": drawn support '", paste(got, collapse=","),
           "' differs from origin_key_nodes.tsv '", paste(want, collapse=","), "'")
    }
  }
  d
}

support_at_least_70 <- function(x) {
  numeric_value <- suppressWarnings(as.numeric(as.character(x)))
  ifelse(!is.na(numeric_value) & numeric_value >= 70,
         as.character(x), "")
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

for (t in trees) {
  cat("Rendering", t$prefix, "\n")

  domain_height <- t$dom_height

  # Ladderization changes only child display order. The rooted input treefile
  # remains byte-frozen and is never rewritten.
  tr <- ape::ladderize(read.tree(file.path(TREE_DIR, t$file)), right=TRUE)

  ladder_plot <- ggtree(tr, branch.length="none", size=t$branch_size)
  family <- t$family
  scope <- t$scope
  centered <- center_focal_rendered(ladder_plot, tr, family, scope)
  p <- centered$plot
  rendered_order_rows[[length(rendered_order_rows) + 1]] <- centered$evidence
  ladderization_rows[[length(ladderization_rows) + 1]] <- ladderization_qc(
    ladder_plot, p, tr, family, scope, centered$exempt_identity_keys
  )
  if (scope == "s1") {
    s1_sister_rows[[length(s1_sister_rows) + 1]] <- s1_direct_sister_evidence(tr, family)
  }

  # ── Pre-compute label dimensions so the clade shading can wrap them ──
  tip_data_pre <- p$data[p$data$isTip, ]
  max_tip_x_pre <- max(tip_data_pre$x)
  is_pruned_pre <- isTRUE(t$pruned)
  if (is_pruned_pre) {
    tip_disp_pre <- sapply(tip_data_pre$label, function(l)
                           make_display_label(l, remove_accession = TRUE))
    max_label_chars_pre <- max(nchar(tip_disp_pre), na.rm = TRUE)
    label_x_start_pre   <- max_tip_x_pre + label_left_pad
    label_width_pre     <- max_label_chars_pre * units_per_char
    shade_xmin_override <- label_x_start_pre - 0.5
    shade_xmax_override <- label_x_start_pre + label_width_pre + 0.5
  } else {
    shade_xmin_override <- NA_real_
    shade_xmax_override <- NA_real_
  }

  # ── Clade shading (wraps tip labels for pruned trees) ──
  shade <- build_shade_rects(p, strip_width)
  if (is_pruned_pre && !is.na(shade_xmin_override)) {
    shade$xmin <- shade_xmin_override
    shade$xmax <- shade_xmax_override
  }
  cat("  Shading blocks:", nrow(shade), "\n")

  for (i in seq_len(nrow(shade))) {
    p <- p + annotate("rect",
      xmin = shade$xmin[i], xmax = shade$xmax[i],
      ymin = shade$ymin[i], ymax = shade$ymax[i],
      fill = clade_colors[shade$cat[i]],
      alpha = 0.30
    )
  }

  # ── Domain architecture ──
  dom_df <- parse_domains(file.path(DOMAIN_DIR, t$domains))
  cat("  Domain blocks:", nrow(dom_df), "\n")

  tip_data <- p$data[p$data$isTip, ]
  max_tip_x <- max(tip_data$x)
  is_pruned <- isTRUE(t$pruned)

  if (is_pruned) {
    # ── Per-protein anchor: leftmost occurrence of dominant catalytic letter ──
    # Falls back to anchor=1 (no shift) when the dominant letter is absent.
    dom_letter <- t$dominant
    anchors <- sapply(split(dom_df, dom_df$label), function(g) {
      hit <- g[g$domain == dom_letter, ]
      if (nrow(hit) == 0) return(1L)
      min(hit$start)
    })
    prot_lens <- sapply(split(dom_df, dom_df$label), function(g) g$prot_len[1])

    max_prot          <- max(prot_lens, na.rm = TRUE)
    max_anchor        <- max(anchors,   na.rm = TRUE)

    # ── Layout: tree | strip | LABELS | DOMAINS ──
    # Compute label width from the longest display label among tips.
    tip_labels       <- tip_data$label
    tip_disp         <- sapply(tip_labels, function(l) make_display_label(l, remove_accession = TRUE))
    max_label_chars  <- max(nchar(tip_disp), na.rm = TRUE)
    label_x_start    <- max_tip_x + label_left_pad  # match shading override
    label_width      <- max_label_chars * units_per_char
    domain_x_start   <- label_x_start + label_width + label_right_pad
    anchor_x         <- domain_x_start + (max_anchor - 1) / aa_per_unit

    # Protein backbones at REAL aa length, anchored on dominant-domain start.
    for (i in seq_len(nrow(tip_data))) {
      lab <- tip_data$label[i]
      a <- anchors[lab]; if (is.na(a)) a <- 1
      plen <- prot_lens[lab]; if (is.na(plen)) plen <- 0
      x0 <- anchor_x - (a - 1) / aa_per_unit
      x1 <- anchor_x + (plen - a + 1) / aa_per_unit
      p <- p + annotate("segment",
        x = x0, xend = x1,
        y = tip_data$y[i], yend = tip_data$y[i],
        color = "grey70", linewidth = 0.3)
    }

    # Domain blocks shifted so each protein's leftmost dominant-letter == anchor_x
    for (i in seq_len(nrow(dom_df))) {
      lab <- dom_df$label[i]
      tip_row <- tip_data[tip_data$label == lab, ]
      if (nrow(tip_row) == 0) next
      tip_y <- tip_row$y[1]
      a <- anchors[lab]; if (is.na(a)) a <- 1
      xmin <- anchor_x + (dom_df$start[i] - a) / aa_per_unit
      xmax <- anchor_x + (dom_df$end[i]   - a + 1) / aa_per_unit
      col <- domain_colors[dom_df$domain[i]]
      if (is.na(col)) col <- "grey50"
      p <- p + annotate("rect",
        xmin = xmin, xmax = xmax,
        ymin = tip_y - domain_height, ymax = tip_y + domain_height,
        fill = col, color = "grey40", linewidth = 0.15)
    }

    # Bound the panel by the farthest backbone actually rendered. Using the
    # independent maxima of protein length and anchor position adds false space.
    domain_right_edge <- max(anchor_x + (prot_lens - anchors + 1) / aa_per_unit, na.rm=TRUE)
    label_x           <- label_x_start  # labels start here, not after domains
  } else {
    # Legacy normalized layout for full (unpruned) trees — too many tips
    # for per-protein backbones to be useful.
    seq_len_total  <- max(dom_df$end)
    domain_x_start <- max_tip_x + domain_offset
    for (i in seq_len(nrow(tip_data))) {
      p <- p + annotate("segment",
        x = domain_x_start, xend = domain_x_start + domain_width,
        y = tip_data$y[i], yend = tip_data$y[i],
        color = "grey70", linewidth = 0.3)
    }
    for (i in seq_len(nrow(dom_df))) {
      tip_row <- tip_data[tip_data$label == dom_df$label[i], ]
      if (nrow(tip_row) == 0) next
      tip_y <- tip_row$y[1]
      xmin <- domain_x_start + (dom_df$start[i] / seq_len_total) * domain_width
      xmax <- domain_x_start + (dom_df$end[i]   / seq_len_total) * domain_width
      col <- domain_colors[dom_df$domain[i]]
      if (is.na(col)) col <- "grey50"
      p <- p + annotate("rect",
        xmin = xmin, xmax = xmax,
        ymin = tip_y - domain_height, ymax = tip_y + domain_height,
        fill = col, color = "grey40", linewidth = 0.15)
    }
    domain_right_edge <- domain_x_start + domain_width
    label_x           <- NA_real_
  }

  # ── Colored domain header above domain column (pruned only) ──
  if (FALSE) {  # Include all domains in legend
    dom_letter <- t$dominant
    p <- p + annotate("text",
      x = domain_x_start + domain_width / 2,
      y = -0.5,
      label = domain_names[dom_letter],
      color = domain_colors[dom_letter],
      size = t$tip_size * 1.2, fontface = "bold", hjust = 0.5
    )
  }

  # ── Domain legend (skip dominant domain for pruned trees) ──
  tree_domains <- sort(unique(dom_df$domain))
  if (FALSE) {  # Include all domains in legend
    tree_domains <- tree_domains  # keep all domains
  }
  n_tips <- sum(p$data$isTip)
  legend_y_start <- n_tips + 2
  legend_cols <- 3
  col_spacing <- if (is_pruned) 13 else 15
  legend_row_h <- if (t$scope == "figure3") 2.4 else 1.9
  # Centre the full key, including the last swatch and its widest label, under the
  # architecture column (Figure 2 and Figure S3 alike).
  legend_labels <- vapply(tree_domains, function(letter) {
    name <- domain_names[letter]
    if (is.na(name)) letter else name
  }, character(1))
  # The saved plot is t$w inches wide; the panel occupies roughly its width
  # after the left labels and ggplot margins. The panel preflight verifies it.
  key_label_units_dom <- key_label_units(p, legend_labels, t$legend_size, t$w)
  # The last text extends beyond the swatch anchor; account for the panel's
  # right-side plotting clearance when centring the complete visible extent.
  right_clearance <- if (length(tree_domains) == 5) 4 else 4.4
  legend_x_start <- centred_key_x(domain_x_start, domain_right_edge, length(tree_domains),
                                  legend_cols, col_spacing, key_label_units_dom, right_clearance)

  for (idx in seq_along(tree_domains)) {
    letter <- tree_domains[idx]
    col_i <- (idx - 1) %% legend_cols
    row_i <- (idx - 1) %/% legend_cols
    x <- legend_x_start + col_i * col_spacing
    y <- legend_y_start + row_i * legend_row_h

    dname <- domain_names[letter]
    if (is.na(dname)) dname <- letter
    dcol <- domain_colors[letter]
    if (is.na(dcol)) dcol <- "grey50"

    p <- p + annotate("rect",
      xmin = x, xmax = x + 3.0,
      ymin = y - 0.5, ymax = y + 0.5,
      fill = dcol
    ) + annotate("text",
      x = x + 3.5, y = y,
      label = dname,
      hjust = 0, size = t$legend_size
    )
  }

  # ── Clade taxonomy key (categories actually present in this tree) ──
  clade_display <- c(
    cynipoid = "Cynipoid wasps", midge = "Midge (Diptera)",
    bee      = "Bee",            beetle = "Beetle",
    nematode = "Nematode",       oomycete = "Oomycete",
    fungus   = "Fungus",         plant   = "Plant",
    protist  = "Protist"
  )
  clade_legend_order <- c("cynipoid","bee","beetle","midge",
                          "nematode","fungus","plant",
                          "oomycete","protist")
  present_clades <- intersect(clade_legend_order, unique(shade$cat))
  if (length(present_clades) > 0 && !grepl('^origin_trees_(GH5_2|PL1_8|PL4_4)$', t$prefix)) {
    domain_key_rows <- ceiling(length(tree_domains) / legend_cols)
    # Full origin panels reserve a separate lower row for the framed taxon key.
    taxon_y_start <- legend_y_start + (domain_key_rows + 1) * legend_row_h
    p <- add_framed_taxon_key(
      p, present_clades, unname(clade_display[present_clades]), clade_colors,
      domain_x_start, domain_right_edge,
      y_start = taxon_y_start,
      cols = legend_cols, col_spacing = col_spacing, row_h = legend_row_h,
      legend_size = t$legend_size, width_in = t$w, right_clearance = right_clearance,
      left_x = 0, swatch_width = 5, swatch_half_height = 0.7)
  }

  # Create display labels: strip type prefix; also strip accession on pruned trees
  p$data$display_label <- sapply(p$data$label, function(l) {
    make_display_label(l, remove_accession = is_pruned)
  })
  p$data$display_label_math <- vapply(
    p$data$display_label, scientific_plotmath, character(1)
  )

  if (is_pruned) {
    # `geom_tiplab(x=X, hjust=0)` is ignored by this ggtree version —
    # labels draw at the default tip position and crash into the
    # domain strip.  Use `offset=` which ggtree honors reliably.
    tip_offset <- label_x - max_tip_x
    # Figure 2 panels use the key-node support rule; Figure S1 keeps plain >=70 labels.
    support_layer <- if (t$scope == "figure3") {
      fig2_labels <- figure2_support_labels(p, tr, t)
      drawn_key_rows[[length(drawn_key_rows) + 1]] <- data.frame(
        family=t$family, node=fig2_labels$role[nzchar(fig2_labels$role)],
        drawn_support=fig2_labels$lab[nzchar(fig2_labels$role)], stringsAsFactors=FALSE)
      drawn_label_rows[[length(drawn_label_rows) + 1]] <- data.frame(
        family=t$family, node=fig2_labels$node, role=fig2_labels$role,
        support=fig2_labels$support, drawn=fig2_labels$lab, stringsAsFactors=FALSE)
      geom_text(data=fig2_labels[nzchar(fig2_labels$lab), ], aes(x=x, y=y, label=lab),
                inherit.aes=FALSE, size=t$node_size, hjust=1.3, vjust=1.5,
                color="grey30")
    } else {
      geom_nodelab(aes(label=support_at_least_70(label)),
                   size=t$node_size, hjust=1.3, vjust=1.5, color="grey30")
    }
    p <- p +
      scale_y_reverse() +
      geom_tiplab(aes(label=display_label_math),
                  size=t$tip_size, offset=tip_offset, hjust=0, parse=TRUE) +
      support_layer +
      xlim_tree(domain_right_edge)
  } else {
    p <- p +
      scale_y_reverse() +
      geom_tiplab(aes(label=display_label_math), size=t$tip_size, offset=0.3,
                  parse=TRUE) +
      geom_nodelab(aes(label=support_at_least_70(label)),
                   size=t$node_size, hjust=1.3, vjust=1.5,
                   color="grey30") +
      xlim_tree(domain_x_start + domain_width + 2)
  }

  outpdf <- file.path(FIG_OUT, paste0(t$prefix, "_ggtree.pdf"))
  outsvg <- file.path(FIG_OUT, paste0(t$prefix, "_ggtree.svg"))
  outpng <- file.path(FIG_OUT, paste0(t$prefix, "_ggtree.png"))

  ggsave(outpdf, p, width=t$w, height=t$h, limitsize=FALSE)
  ggsave(outsvg, p, width=t$w, height=t$h, limitsize=FALSE)
  ggsave(outpng, p, width=t$w, height=t$h, dpi=t$png_dpi, bg="white", limitsize=FALSE)
  cat("  ->", outpdf, outsvg, outpng, "\n")

  if (isTRUE(t$pruned)) pruned_plots[[t$prefix]] <- p
}

# ── Combined subset figure: PL1_8, PL4_4, GH5_2 top to bottom ──
# CBM13/B_lectin moved to its own dedicated Fig 8 (tree + domain + chr locus
# in one panel, monophyletic cynipoid clade, donor-architecture bacteria only).
library(patchwork)
compact_combo_panel <- function(p) {
  # Text has physical width beyond its data anchor. Preserve panel geometry
  # while allowing the domain-key labels to use the remaining page margin.
  p$coordinates$clip <- "off"
  p + theme(plot.margin = margin(1, -6, 1, -6, unit = "pt"))
}
clade_key <- data.frame(
  cat = factor(names(clade_colors), levels=names(clade_colors)),
  label = c('Cynipoid wasps', 'Midge (Diptera)', 'Bee', 'Beetle',
            'Nematode', 'Oomycete', 'Fungus', 'Plant', 'Protist'),
  x = rep(1:3, each=3), y = rep(3:1, 3)
)
key <- ggplot(clade_key, aes(x=x, y=y, fill=cat)) +
  annotate('rect', xmin=.99, xmax=4.08, ymin=.52, ymax=3.48,
           fill=NA, color='grey40', linewidth=.25) +
  geom_rect(aes(xmin=x, xmax=x+.20, ymin=y-.30, ymax=y+.30), alpha=.45) +
  geom_text(aes(x=x+.25, label=label), hjust=0, size=2.4) +
  scale_fill_manual(values=clade_colors) +
  coord_cartesian(xlim=c(.97,4.1), ylim=c(.5,3.5), expand=FALSE) +
  theme_void() + theme(legend.position='none', plot.margin=margin(0,0,0,0))
combo <- (compact_combo_panel(pruned_plots[["origin_trees_PL1_8"]]) /
  compact_combo_panel(pruned_plots[["origin_trees_PL4_4"]]) /
  compact_combo_panel(pruned_plots[["origin_trees_GH5_2"]]) /
  ((key | plot_spacer()) + plot_layout(widths=c(.75,.25)))) +
  plot_layout(heights=c(1,1,1,.18))
ggsave(file.path(FIG_OUT, "origin_trees.pdf"), combo, width=6.89, height=9.25, limitsize=FALSE)
canonicalize_figure3_pdf_timestamp(file.path(FIG_OUT, "origin_trees.pdf"))
ggsave(file.path(FIG_OUT, "origin_trees.svg"), combo, width=6.89, height=9.25, limitsize=FALSE)
ggsave(file.path(FIG_OUT, "origin_trees.png"), combo, width=6.89, height=9.25, dpi=300, bg="white", limitsize=FALSE)
cat("  -> ", FIG_OUT, "/origin_trees.{pdf,svg,png}\n", sep="")

rendered_order_qc <- do.call(rbind, rendered_order_rows)
write.table(rendered_order_qc,
            file=file.path(EVIDENCE_DIR, "rendered_tip_order_qc.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
if (nrow(rendered_order_qc) != 6 || any(rendered_order_qc$outer_15pct_gate != "PASS")) {
  stop("renderer-native y-coordinate QC did not pass all six panels")
}

ladderization_qc_table <- do.call(rbind, ladderization_rows)
write.table(ladderization_qc_table,
            file=file.path(EVIDENCE_DIR, "ladderization_qc.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE, na="")
if (any(ladderization_qc_table$ladderization_gate == "FAIL")) {
  stop("ladderization QC contains a failed internal node")
}

s1_sister_qc <- do.call(rbind, s1_sister_rows)
write.table(s1_sister_qc,
            file=file.path(EVIDENCE_DIR, "s1_direct_sister_species.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE, na="not defined")
if (nrow(s1_sister_qc) != 3 || any(s1_sister_qc$direct_sister_tip_count < 1)) {
  stop("S1 direct-sister evidence is incomplete")
}

drawn_keys <- do.call(rbind, drawn_key_rows)
if (is.null(drawn_keys) || nrow(drawn_keys) != 6) stop("Figure 2 key-node labels missing")
write.table(drawn_keys, file=file.path(EVIDENCE_DIR, "figure2_key_node_labels.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(do.call(rbind, drawn_label_rows),
            file=file.path(EVIDENCE_DIR, "figure2_support_labels.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

cat("Done!\n")

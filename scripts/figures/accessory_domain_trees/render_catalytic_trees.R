library(ggtree)
library(treeio)
library(ggplot2)
library(ape)
options(warn=1)

INPUT_DIR <- Sys.getenv("PCWDE_ACCESSORY_DOMAIN_TREES_INPUT_DIR", unset="")
FIG_OUT <- Sys.getenv("PCWDE_ACCESSORY_DOMAIN_TREES_OUTPUT_DIR", unset="")
if (!nzchar(INPUT_DIR)) stop("PCWDE_ACCESSORY_DOMAIN_TREES_INPUT_DIR is required")
if (!nzchar(FIG_OUT)) stop("PCWDE_ACCESSORY_DOMAIN_TREES_OUTPUT_DIR is required")
if (dir.exists(FIG_OUT) && length(list.files(FIG_OUT, all.files=TRUE, no..=TRUE))) {
  stop("refusing non-empty output directory: ", FIG_OUT)
}
dir.create(FIG_OUT, showWarnings=FALSE, recursive=TRUE)
# Shared key layout (scripts/figures/key_layout.R), one directory above this renderer.
.renderer_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
if (length(.renderer_file) != 1) stop("cannot locate the running renderer script")
source(file.path(dirname(dirname(normalizePath(.renderer_file))), "key_layout.R"))

# Display-only S4 crosswalks from the validated current review authority.
S4_INPUT <- file.path(INPUT_DIR, "label_inputs")
METADATA_TSV <- file.path(INPUT_DIR, "metadata.tsv")
s4_labels <- c()
metadata_labels <- c()
if (nchar(METADATA_TSV) > 0) {
  metadata <- read.delim(METADATA_TSV, stringsAsFactors=FALSE,
                         check.names=FALSE, quote="")
  metadata_labels <- setNames(metadata$species, metadata$protein_id)
}
if (nchar(S4_INPUT) > 0) {
  # `tip_crosswalk.tsv` is the complete 246-tip identifier authority, while
  # `display_labels.tsv` is its row-aligned public-label authority.  The old
  # adapter treated every crosswalk row as focal and read its `species`
  # column; that field is intentionally blank for all 223 reference tips.
  # Appending reference_panel_manifest.tsv then created duplicate names, and
  # named-vector lookup returned the first (blank) value.  Join the two
  # complete authorities directly and fail closed if their alignment drifts.
  label_contract <- file.path(INPUT_DIR, "s4_label_contract.R")
  source(label_contract, local=TRUE, encoding="UTF-8")
  crosswalk <- read_s4_tsv_utf8(file.path(S4_INPUT, "tip_crosswalk.tsv"))
  accepted_labels <- read_s4_tsv_utf8(
    file.path(S4_INPUT, "display_labels.tsv")
  )
  label_map <- build_s4_label_map(crosswalk, accepted_labels)
  s4_labels <- setNames(label_map$display_label, label_map$tip_id)
}

# ── Eukaryote category detection ──────────────────────────────
euk_prefixes <- c("cynipoid", "midge", "bee", "beetle",
                  "nematode", "oomycete", "fungus", "plant", "protist")

get_category <- function(label) {
  if (grepl("^[^|]*_ENA\\|CDX", label)) return("cynipoid")
  if (grepl("^[A-Z0-9_]+_F[0-9]+__", label)) return("cynipoid")
  if (grepl("__bacterium_", label)) return(NA_character_)
  if (grepl("__fungus_", label)) return("fungus")
  if (grepl("__plant_", label)) return("plant")
  prefix <- sub("\\|.*", "", label)
  if (prefix %in% euk_prefixes) return(prefix)
  return(NA_character_)
}

# Strip type prefix and optionally accession from tip labels
make_display_label <- function(label, remove_accession = FALSE) {
  if (label %in% names(s4_labels)) return(unname(s4_labels[label]))
  if (label %in% names(metadata_labels)) {
    cleaned <- unname(metadata_labels[label])
    cleaned <- gsub("Isocolus centaureae", "Eubothrus centaureae", cleaned, fixed=TRUE)
    cleaned <- gsub("Ceroptres masudai", "Ceroptres kovalevi", cleaned, fixed=TRUE)
    cleaned <- gsub("Bradysia coprophila", "Bradysia tilicola", cleaned, fixed=TRUE)
    return(cleaned)
  }
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
    if (grepl("Synergus__BK7_inquiline", cleaned)) return("Synergus sp. ex Belonocnema kinseyi")
    if (grepl("Ceroptres.*DI.*inquiline", cleaned)) return("Ceroptres sp. ex Druon ignotum")
    if (grepl("Synergus.*quercuslanigerum.*inquiline|Synergus.*DL.*inquiline", cleaned)) return("Synergus sp. ex Druon quercuslanigerum")
    if (grepl("Synergus__DI_inquiline|Synergus.*DC.*inquiline|Synergus__DC", cleaned)) return("Synergus sp. ex Druon ignotum")
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

# Preserve the accepted unobtrusive paralog suffix policy.  Existing
# cynipoid labels carry their stable paralog number; any remaining collision
# receives a deterministic bracketed rank based on raw tip ID.
make_display_labels <- function(labels, remove_accession = FALSE) {
  displayed <- vapply(labels, make_display_label, character(1),
                      remove_accession=remove_accession)
  explicit <- grepl("^cynipoid\\|.*_[0-9]+$", labels)
  suffix <- sub(".*_([0-9]+)$", "\\1", labels)
  displayed[explicit] <- paste0(displayed[explicit], " [", suffix[explicit], "]")
  for (value in unique(displayed[duplicated(displayed)])) {
    idx <- which(displayed == value)
    order_idx <- idx[order(labels[idx])]
    displayed[order_idx] <- paste0(value, " [", seq_along(order_idx), "]")
  }
  unname(displayed)
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
  B = "#E6AB02",   # CBM2 - Amber
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
  q = "#7A0D16",   # Region A only — dark crimson, custom HMM donor-region residual
  o = "#A0A0A0"    # Other Pfam domain
)

domain_names <- c(
  A = "GH5_2",  B = "CBM2",  C = "CBM3",  D = "CBM5",
  E = "CBM6",   F = "CBM10", G = "CBM12", H = "CBM17",
  I = "CBM28",  J = "CBM73", K = "CBM81", L = "PL1_8",
  M = "B-lectin (PF01453)", N = "PL4", O = "PL4_2", P = "PL4_4",
  Q = "PL4_5",  R = "CBM22", S = "CBM32", T = "CBM64",
  U = "GH30",   V = "GH30_C",W = "CBM4_9",X = "fn3",
  Y = "Ricin_B",Z = "Big_5",
  u = "GH59",   v = "GH30_2",w = "CBM6",  y = "RicinB_2",
  c = "CBM35_2",s = "SLH",   p = "PI-PLC-X",
  m = "CBM13 (weak)", q = "Region A only", o = "Other Pfam"
)

# The S4 architecture table is generated once from a uniform full-protein
# Pfam scan. Every tip has a protein length and at least one mapped call.
parse_s4_architecture <- function(path) {
  x <- read.delim(path, stringsAsFactors=FALSE, check.names=FALSE)
  required <- c("tip_id", "protein_length_aa", "domain", "start", "end")
  stopifnot(all(required %in% names(x)), nrow(x) > 0,
            all(x$start >= 1), all(x$end >= x$start),
            all(x$end <= x$protein_length_aa))
  code <- c(B_lectin="M", Glyco_hydro_30="U", Glyco_hydro_30C="V",
            `PI-PLC-X`="p", Glyco_hydro_59="u", CBM_4_9="W",
            CBM_6="w", Ricin_B_lectin="Y", RicinB_lectin_2="y",
            fn3="X", SLH="s")
  x$code <- unname(code[x$domain]); x$code[is.na(x$code)] <- "o"
  data.frame(label=x$tip_id, start=x$start, end=x$end,
             domain=x$code, prot_len=x$protein_length_aa,
             stringsAsFactors=FALSE)
}

plotmath_quote <- function(x) encodeString(x, quote='"')
italicize_s4_label <- function(x) {
  # "unclassified Genus[ - accession]": qualifier roman, genus italic.
  uc <- regmatches(x, regexec("^(unclassified )([A-Z][a-z]+)(( - .*)?)$", x))[[1]]
  if (length(uc) == 5L)
    return(paste0("plain(", plotmath_quote(uc[2]), ")*italic(",
                  plotmath_quote(uc[3]), ")*plain(", plotmath_quote(uc[4]), ")"))
  # "Genus sp. ex Host species[ - ...]": both taxa italic, "sp. ex" roman.
  ex <- regmatches(x, regexec(
    "^([A-Z][a-z]+)( sp\\. ex )([A-Z][a-z]+ [a-z][a-z-]+)(.*)$", x))[[1]]
  if (length(ex) == 5L)
    return(paste0("italic(", plotmath_quote(ex[2]), ")*plain(",
                  plotmath_quote(ex[3]), ")*italic(", plotmath_quote(ex[4]),
                  ")*plain(", plotmath_quote(ex[5]), ")"))
  sp <- regexec("^([A-Z][A-Za-z-]+)( sp\\..*)$", x)
  parts <- regmatches(x, sp)[[1]]
  if (length(parts) == 3L)
    return(paste0("italic(", plotmath_quote(parts[2]), ")*plain(",
                  plotmath_quote(parts[3]), ")"))
  match <- regexec("^([A-Z][A-Za-z-]+) ([a-z][a-z-]+)(.*)$", x)
  parts <- regmatches(x, match)[[1]]
  if (length(parts) == 4L) {
    return(paste0("italic(", plotmath_quote(paste(parts[2], parts[3])),
                  ")*plain(", plotmath_quote(parts[4]), ")"))
  }
  if (grepl("^[A-Z][a-z]+$", x))
    return(paste0("italic(", plotmath_quote(x), ")"))
  paste0("plain(", plotmath_quote(x), ")")
}

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
  # Category values must be unnamed: sapply()/vapply() name results by tip
  # label, and identical() compares those names, so adjacent same-category
  # tips with different labels never coalesced into one block.
  tip_data$cat <- vapply(tip_data$label, get_category, character(1),
                         USE.NAMES=FALSE)

  pad <- 0.49

  rects <- list()
  cur_cat <- NA
  cur_ymin <- NA
  cur_ymax <- NA

  for (i in seq_len(nrow(tip_data))) {
    cat_i <- tip_data$cat[i]
    yi <- tip_data$y[i]

    if (!is.na(cat_i) && identical(unname(cat_i), unname(cur_cat))) {
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

# Contiguous same-category tip runs in y order; unshaded (NA) tips break runs.
# Tip IDs are kept separately from the unnamed category values.
shade_runs <- function(tip_data) {
  ordered <- tip_data[order(tip_data$y), ]
  cats <- vapply(ordered$label, get_category, character(1), USE.NAMES=FALSE)
  n <- length(cats)
  new_run <- c(TRUE, if (n > 1L) vapply(seq.int(2L, n), function(j)
    is.na(cats[j]) || is.na(cats[j - 1L]) || cats[j] != cats[j - 1L],
    logical(1)) else logical(0))
  runs <- split(seq_len(n), cumsum(new_run))
  runs <- Filter(function(idx) !is.na(cats[idx[1]]), runs)
  lapply(unname(runs), function(idx) list(
    cat = cats[idx[1]], tips = ordered$label[idx], y = ordered$y[idx]))
}

# Compare every rectangle with the tip run it should cover.  Fails when the
# block count differs, a category differs, or either edge is more than
# `bound` rows from the run's tip extent padded by `pad`.
check_shade_geometry <- function(shade, tip_data, prefix, pad=0.49, bound=0.5) {
  runs <- shade_runs(tip_data)
  cat(sprintf("  SHADE_COUNT %s producer_blocks=%d independent_runs=%d\n",
              prefix, nrow(shade), length(runs)))
  if (length(runs) != nrow(shade))
    stop("S4 shade block count ", nrow(shade), " != tip runs ", length(runs),
         " in ", prefix)
  out <- data.frame(block=seq_along(runs), class=NA_character_, tips=0L,
                    ymin=shade$ymin, ymax=shade$ymax,
                    tip_ymin=NA_real_, tip_ymax=NA_real_,
                    mismatch_rows=NA_real_, stringsAsFactors=FALSE)
  for (i in seq_along(runs)) {
    r <- runs[[i]]
    expected <- range(r$y) + c(-pad, pad)
    mismatch <- max(abs(c(shade$ymin[i], shade$ymax[i]) - expected))
    out$class[i] <- r$cat
    out$tips[i] <- length(r$tips)
    out$tip_ymin[i] <- min(r$y)
    out$tip_ymax[i] <- max(r$y)
    out$mismatch_rows[i] <- mismatch
    cat(sprintf("  SHADE_GEOMETRY %s block=%d class=%s tips=%d mismatch_rows=%.3f\n",
                prefix, i, r$cat, length(r$tips), mismatch))
    if (!identical(unname(shade$cat[i]), unname(r$cat)) || mismatch > bound)
      stop("S4 shade extent/category mismatch in ", prefix, " block ", i)
  }
  invisible(out)
}

# ── Tree configs ──────────────────────────────────────────────
# Two accepted catalytic-origin panels.
pruned_trees <- list(
  list(prefix="accessory_domain_trees_GH30_3", file="GH30_3.unrooted.treefile",
       domains="architecture.tsv", w=21, h=42, tip_size=4.2, node_size=3.1,
       legend_size=6, branch_size=0.35, dom_height=0.25, png_dpi=160,
       pruned=FALSE, dominant="U"),
  list(prefix="accessory_domain_trees_PI_PLC_X", file="PI_PLC_X.unrooted.treefile",
       domains="architecture.tsv", w=21, h=19, tip_size=8, node_size=5.2,
       legend_size=7, branch_size=0.5, dom_height=0.3, png_dpi=200,
       pruned=FALSE, dominant="p")
)

trees <- pruned_trees

pruned_plots <- list()   # store pruned plots for combined figure

strip_width     <- 5      # narrow clade-color band between tips and domain strip
domain_offset   <- 18    # legacy: floor for unpruned trees only
domain_width    <- 12    # legacy: only used for full (unpruned) trees
aa_per_unit     <- 50    # 1 x-unit of domain strip == 50 aa of protein (pruned trees)
label_gap       <- 4     # gap between label area and domain strip
units_per_char  <- 0.62  # conservative text-width estimate at tip_size 8
label_left_pad  <- 2     # gap between strip end and label start
label_right_pad <- 5     # explicit clearance between longest label and domain strip

for (t in trees) {
  cat("Rendering", t$prefix, "\n")

  domain_height <- t$dom_height

  tr <- read.tree(file.path(INPUT_DIR, t$file))
  tr <- ladderize(tr, right=FALSE)

  p <- ggtree(tr, branch.length="none", size=t$branch_size)

  # ── Pre-compute label dimensions so the clade shading can wrap them ──
  tip_data_pre <- p$data[p$data$isTip, ]
  max_tip_x_pre <- max(tip_data_pre$x)
  is_pruned_pre <- isTRUE(t$pruned)
  if (is_pruned_pre) {
    tip_disp_pre <- make_display_labels(tip_data_pre$label, remove_accession=TRUE)
    max_label_chars_pre <- max(nchar(tip_disp_pre), na.rm = TRUE)
    label_x_start_pre   <- max_tip_x_pre + label_left_pad
    label_width_pre     <- max_label_chars_pre * units_per_char
    shade_xmin_override <- label_x_start_pre - 0.5
    shade_xmax_override <- label_x_start_pre + label_width_pre + 0.5
  } else {
    shade_xmin_override <- NA_real_
    shade_xmax_override <- NA_real_
  }

  # ── Clade shading ──
  shade <- build_shade_rects(p, strip_width)
  if (!is_pruned_pre) {
    # Full S4 labels occupy the column between the tree and architecture.
    # The old five-unit strip ended near the start of each name.
    shade$xmin <- max_tip_x_pre
    shade$xmax <- max_tip_x_pre + label_left_pad +
      max(nchar(make_display_labels(tip_data_pre$label,
                                    remove_accession=FALSE)), na.rm=TRUE) *
      units_per_char * (t$tip_size / 8) + 0.5
  }
  if (is_pruned_pre && !is.na(shade_xmin_override)) {
    shade$xmin <- shade_xmin_override
    shade$xmax <- shade_xmax_override
  }
  cat("  Shading blocks:", nrow(shade), "\n")
  check_shade_geometry(shade, tip_data_pre, t$prefix)

  for (i in seq_len(nrow(shade))) {
    p <- p + annotate("rect",
      xmin = shade$xmin[i], xmax = shade$xmax[i],
      ymin = shade$ymin[i], ymax = shade$ymax[i],
      fill = clade_colors[shade$cat[i]],
      alpha = 0.30
    )
  }

  # ── Domain architecture ──
  has_domains <- !is.na(t$domains)
  dom_df <- if (has_domains) parse_s4_architecture(file.path(INPUT_DIR, t$domains)) else data.frame(
    label=character(), start=integer(), end=integer(), domain=character(),
    prot_len=integer(), stringsAsFactors=FALSE)
  cat("  Domain blocks:", nrow(dom_df), "\n")

  tip_data <- p$data[p$data$isTip, ]
  max_tip_x <- max(tip_data$x)
  is_pruned <- isTRUE(t$pruned)

  if (!has_domains) {
    # S4 has no protein-architecture track in its validated authority.
    domain_x_start <- max_tip_x + domain_offset
    domain_right_edge <- domain_x_start
    label_x <- NA_real_
  } else if (is_pruned) {
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
    tip_disp         <- make_display_labels(tip_labels, remove_accession=TRUE)
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

    domain_right_edge <- anchor_x + (max_prot - 1) / aa_per_unit
    label_x           <- label_x_start  # labels start here, not after domains
  } else {
    # Full S4 trees: tree | names | source proteins at a common aa scale, each protein
    # aligned on the start of its catalytic domain (GH30 or PI-PLC-X), as Figure S3 anchors
    # on the dominant catalytic domain. Whole proteins are drawn; nothing is clipped.
    tip_disp <- make_display_labels(tip_data$label, remove_accession=FALSE)
    label_x <- max_tip_x + label_left_pad
    label_width <- max(nchar(tip_disp), na.rm=TRUE) * units_per_char *
                   (t$tip_size / 8)
    domain_x_start <- label_x + label_width + label_right_pad
    aa_per_unit <- 50
    anchors <- sapply(split(dom_df, dom_df$label), function(g) {
      hit <- g[g$domain == t$dominant, ]
      if (nrow(hit) == 0) return(NA_integer_)
      min(hit$start)
    })
    missing_anchor <- setdiff(tip_data$label, names(anchors)[!is.na(anchors)])
    if (length(missing_anchor)) {
      stop("S4 tips without their catalytic domain: ", paste(head(missing_anchor), collapse=", "))
    }
    prot_lens <- sapply(split(dom_df, dom_df$label), function(g) {
      stopifnot(length(unique(g$prot_len)) == 1L)
      g$prot_len[1]
    })
    anchor_x <- domain_x_start + (max(anchors[tip_data$label]) - 1) / aa_per_unit
    for (i in seq_len(nrow(tip_data))) {
      lab <- tip_data$label[i]
      a <- anchors[[lab]]; plen <- prot_lens[[lab]]
      p <- p + annotate("segment",
        x = anchor_x - (a - 1) / aa_per_unit, xend = anchor_x + (plen - a + 1) / aa_per_unit,
        y = tip_data$y[i], yend = tip_data$y[i],
        color = "grey70", linewidth = 0.3)
    }
    for (i in seq_len(nrow(dom_df))) {
      tip_row <- tip_data[tip_data$label == dom_df$label[i], ]
      if (nrow(tip_row) == 0) next
      tip_y <- tip_row$y[1]
      a <- anchors[[dom_df$label[i]]]
      xmin <- anchor_x + (dom_df$start[i] - a) / aa_per_unit
      xmax <- anchor_x + (dom_df$end[i] - a + 1) / aa_per_unit
      col <- domain_colors[dom_df$domain[i]]
      if (is.na(col)) col <- "grey50"
      p <- p + annotate("rect",
        xmin = xmin, xmax = xmax,
        ymin = tip_y - domain_height, ymax = tip_y + domain_height,
        fill = col, color = "grey40", linewidth = 0.15)
    }
    domain_right_edge <- max(anchor_x + (prot_lens[tip_data$label] -
                                         anchors[tip_data$label] + 1) / aa_per_unit)
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
  col_spacing <- 15
  legend_row_h <- 1.9
  # Domain key centred under the domain bars (shared Figure 2 key layout).
  legend_labels <- vapply(tree_domains, function(letter) {
    name <- domain_names[letter]
    if (is.na(name)) letter else name
  }, character(1))
  right_clearance <- 4.4
  # Measure label widths against the final x range: the xlim_tree(domain_right_edge + 4)
  # added before saving widens the panel, so measuring on p alone underestimates them.
  key_measure_p <- p + xlim_tree(domain_right_edge + 4)
  legend_label_units <- key_label_units(key_measure_p, legend_labels, t$legend_size, t$w)
  # Columns never closer than the widest label plus its swatch (same rule as the taxon key),
  # so a long label such as "B-lectin (PF01453)" cannot run into the next swatch.
  col_spacing <- max(col_spacing, 3.5 + legend_label_units + 1.5)
  legend_x_start <- centred_key_x(domain_x_start, domain_right_edge, length(tree_domains),
                                  legend_cols, col_spacing, legend_label_units,
                                  right_clearance)

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
  if (length(present_clades) > 0) {
    clade_n_rows  <- ceiling(length(tree_domains) / legend_cols)
    p <- add_framed_taxon_key(
      p, present_clades, unname(clade_display[present_clades]), clade_colors,
      domain_x_start, domain_right_edge,
      y_start = legend_y_start + clade_n_rows * legend_row_h + 1.2,
      cols = legend_cols, col_spacing = col_spacing, row_h = legend_row_h,
      legend_size = t$legend_size, width_in = t$w, right_clearance = right_clearance,
      measure_p = key_measure_p)
  }

  # Create display labels: strip type prefix; also strip accession on pruned trees
  p$data$display_label <- make_display_labels(p$data$label,
                                               remove_accession=is_pruned)
  p$data$display_label_math <- vapply(p$data$display_label,
                                      italicize_s4_label, character(1))

  if (is_pruned) {
    # `geom_tiplab(x=X, hjust=0)` is ignored by this ggtree version —
    # labels draw at the default tip position and crash into the
    # domain strip.  Use `offset=` which ggtree honors reliably.
    tip_offset <- label_x - max_tip_x
    p <- p +
      scale_y_reverse() +
      geom_tiplab(aes(label=display_label),
                  size=t$tip_size, offset=tip_offset, hjust=0) +
      geom_nodelab(aes(label=label),
                   size=t$node_size, hjust=1.3, vjust=1.5,
                   color="grey30") +
      xlim_tree(domain_right_edge + 4)
  } else {
    tip_offset <- label_x - max_tip_x
    p <- p +
      scale_y_reverse() +
      geom_tiplab(aes(label=display_label_math), size=t$tip_size,
                  offset=tip_offset, hjust=0, parse=TRUE) +
      geom_nodelab(aes(label=label),
                   size=t$node_size, hjust=1.3, vjust=1.5,
                   color="grey30") +
      xlim_tree(domain_right_edge + 4)
  }

  outpdf <- file.path(FIG_OUT, paste0(t$prefix, "_ggtree.pdf"))
  outsvg <- file.path(FIG_OUT, paste0(t$prefix, "_ggtree.svg"))
  outpng <- file.path(FIG_OUT, paste0(t$prefix, "_ggtree.png"))

  ggsave(outpdf, p, width=t$w, height=t$h, limitsize=FALSE)
  ggsave(outsvg, p, width=t$w, height=t$h, limitsize=FALSE)
  ggsave(outpng, p, width=t$w, height=t$h, dpi=t$png_dpi, bg="white", limitsize=FALSE)
  cat("  ->", outpdf, outsvg, outpng, "\n")

  pruned_plots[[t$prefix]] <- p
}

cat("Done!\n")

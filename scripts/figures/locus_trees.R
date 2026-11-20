library(ggtree)
library(treeio)
library(ggplot2)
library(ape)
library(patchwork)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])), "locus_tip_marks.R"))

# Figure 6 candidate renderer for the corrected locus inputs. It uses the
# approved display roster and geometry with the current tree and source tables.
required_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop(name, " is required")
  value
}
FIG2_INPUT <- normalizePath(
  required_env("PCWDE_LOCUS_TREES_INPUT_DIR"), mustWork = TRUE
)
S3_TREE_DIR <- normalizePath(
  required_env("PCWDE_LOCUS_TREES_COMPLETE_INPUT_DIR"), mustWork = TRUE
)
TREE_DIR <- Sys.getenv(
  "FIG2_TREE_DIR",
  unset = file.path(FIG2_INPUT, "trees")
)
OVERLAY_FILE <- Sys.getenv(
  "FIG2_OVERLAY_FILE",
  unset = file.path(FIG2_INPUT, "all_family_rna_evidence_overlay.tsv")
)
PROJECTION_FILE <- Sys.getenv(
  "FIG2_PROJECTION_FILE",
  unset = file.path(FIG2_INPUT, "approved_roster_projection.tsv")
)
LABEL_FILE <- Sys.getenv(
  "FIG2_LABEL_FILE",
  unset = file.path(FIG2_INPUT, "raw_tip_to_normalized_display_label.tsv")
)
LABEL_OVERRIDE_FILE <- file.path(FIG2_INPUT, "label_overrides.tsv")
FIG_OUT <- required_env("PCWDE_LOCUS_TREES_OUTPUT_DIR")
QC_OUT <- file.path(FIG_OUT, "evidence")
if (dir.exists(FIG_OUT) && length(list.files(
  FIG_OUT, all.files = TRUE, no.. = TRUE
))) {
  stop("refusing non-empty output directory: ", FIG_OUT)
}
dir.create(FIG_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_OUT, recursive = TRUE, showWarnings = FALSE)

required <- c(
  OVERLAY_FILE, PROJECTION_FILE, LABEL_FILE,
  file.path(FIG2_INPUT, "rebuild_root_sides.tsv"),
  file.path(FIG2_INPUT, "species_colors.R")
)
if (any(!file.exists(required))) {
  stop("missing canonical Figure 2 input(s): ",
       paste(required[!file.exists(required)], collapse = ", "))
}

overlay <- read.delim(
  OVERLAY_FILE, stringsAsFactors = FALSE, check.names = FALSE
)
projection <- read.delim(
  PROJECTION_FILE, stringsAsFactors = FALSE, check.names = FALSE
)
labels <- read.delim(
  LABEL_FILE, stringsAsFactors = FALSE, check.names = FALSE
)
# Optional paralog display labels; absent from the current deposit.
label_overrides <- if (file.exists(LABEL_OVERRIDE_FILE)) {
  read.delim(
    LABEL_OVERRIDE_FILE, stringsAsFactors = FALSE, check.names = FALSE
  )
} else {
  data.frame()
}
source(file.path(FIG2_INPUT, "species_colors.R"))

families <- c("GH5_2", "PL1_8", "PL4_4")
# Tip counts follow the deposit's label and roster tables.
expected_full <- setNames(vapply(families, function(f) sum(labels$family == f), integer(1)), families)
expected_s3 <- expected_full
expected_reduced <- setNames(vapply(families, function(f) sum(projection$family == f), integer(1)), families)
full_height <- c(GH5_2 = 60, PL1_8 = 70, PL4_4 = 85)
reduced_outgroups <- list(
  GH5_2 = c(
    "cynipoid|Diplolepis_rosae_1", "cynipoid|Diplolepis_rosae_2",
    "cynipoid|Diplolepis_rosae_5", "cynipoid|Diplolepis_rosae_7",
    "cynipoid|Diplolepis_spinosa_1"
  ),
  PL1_8 = c(
    "cynipoid|Diplolepis_rosae_1", "cynipoid|Diplolepis_rosae_2",
    "cynipoid|Diplolepis_rosae_5", "cynipoid|Diplolepis_rosae_7",
    "cynipoid|Diplolepis_spinosa_1"
  ),
  PL4_4 = c(
    "cynipoid|Diplolepis_rosae_3", "cynipoid|Diplolepis_rosae_4"
  )
)
# Approved display rooting (2026-09-07), re-expressed on the rebuilt trees:
# the root-side tips for PL1_8 and PL4_4; GH5_2 keeps the S3 orientation.
root_sides <- read.delim(file.path(FIG2_INPUT, "rebuild_root_sides.tsv"), stringsAsFactors=FALSE)
stopifnot(all(c("family", "tip") %in% names(root_sides)),
          setequal(unique(root_sides$family), c("PL1_8", "PL4_4")))
display_root_sides <- list(
  GH5_2 = character(0),
  PL1_8 = root_sides$tip[root_sides$family == "PL1_8"],
  PL4_4 = root_sides$tip[root_sides$family == "PL4_4"]
)
if (nrow(labels) != sum(expected_full) ||
    any(labels$technical_token_in_display != "False") ||
    any(!nzchar(labels$normalized_display_label)) ||
    anyDuplicated(paste(labels$family, labels$normalized_display_label)) ||
    any(!nzchar(labels$raw_tree_tip)) ||
    any(!nzchar(labels$tree_protein_id)) ||
    anyDuplicated(paste(labels$family, labels$raw_tree_tip))) {
  stop("DB-derived display-label precondition failed")
}
if (nrow(label_overrides)) {
  required_override_columns <- c(
    "family", "raw_tree_tip", "display_label", "highlight_group"
  )
  if (!all(required_override_columns %in% names(label_overrides)) ||
      nrow(label_overrides) != sum(expected_full) ||
      anyDuplicated(paste(label_overrides$family, label_overrides$raw_tree_tip)) ||
      any(!nzchar(label_overrides$display_label)) ||
      any(!label_overrides$highlight_group %in% c("other", "DC", "BK")) ||
      !setequal(
        paste(label_overrides$family, label_overrides$raw_tree_tip),
        paste(labels$family, labels$raw_tree_tip)
      )) {
    stop("display label override precondition failed")
  }
}

# Use one parsed scientific-name authority for both reduced Figure 2 and the
# complete S2 views. Three accepted bacterial outgroups lack a DB canonical
# name; for those rows,
# the accepted normalized label is the authority after removing only its
# synthetic terminal copy ordinal. Repeated visible names are intentional.
labels$scientific_species_name <- labels$canonical_scientific_name
fallback <- is.na(labels$scientific_species_name) |
  !nzchar(labels$scientific_species_name)
labels$scientific_species_name[fallback] <- sub(
  " [0-9]+$", "", labels$normalized_display_label[fallback]
)
# Current taxonomy is display-only: raw tree tips, protein IDs, and species
# codes remain stable lookup authorities. Apply this after canonical/fallback
# derivation so both paths use the same public names.
current_display_names <- c(
  "Eschatocerus acacia" = "Eschatocerus acaciae",
  "Isocolus centaureae" = "Eubothrus centaureae",
  "Ceroptres masudai" = "Ceroptres kovalevi",
  "Bradysia coprophila" = "Bradysia tilicola"
)
for (stale_name in names(current_display_names)) {
  matched <- labels$scientific_species_name == stale_name
  labels$scientific_species_name[matched] <- current_display_names[[stale_name]]
}
labels$figure2_display_label <- labels$scientific_species_name
labels$s2_display_label <- labels$scientific_species_name
if (any(!nzchar(labels$figure2_display_label)) ||
    any(!nzchar(labels$s2_display_label)) ||
    any(grepl(" [0-9]+$", labels$s2_display_label)) ||
    any(labels$s2_display_label == labels$tree_protein_id) ||
    any(labels$s2_display_label == labels$raw_tree_tip) ||
    anyDuplicated(paste(labels$family, labels$raw_tree_tip))) {
  stop("display-only relabeling gate failed")
}
write.table(
  labels[, c(
    "family", "raw_tree_tip", "tree_protein_id", "species_code",
    "canonical_scientific_name", "normalized_display_label",
    "figure2_display_label", "s2_display_label"
  )],
  file.path(QC_OUT, "tip_display_label_sidecar.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE
)
if (nrow(overlay) != sum(expected_full) ||
    anyDuplicated(paste(overlay$family, overlay$tip))) {
  stop("all-family overlay cardinality/uniqueness precondition failed")
}

chr_display <- function(x) {
  if (is.na(x) || x %in% c(
    "", "outgroup", "unknown", "transcript", "fragment", "unplaced"
  )) return("")
  x <- sub("^chr", "", x)
  value <- suppressWarnings(as.integer(x))
  if (!is.na(value) && value >= 1 && value <= 10) {
    return(as.character(value))
  }
  ""
}

safe_root <- function(tree, requested) {
  outgroup <- requested[requested %in% tree$tip.label]
  if (!length(outgroup)) return(tree)
  original <- tree$tip.label
  candidate <- try(unroot(tree), silent = TRUE)
  if (inherits(candidate, "phylo")) {
    phantom <- setdiff(candidate$tip.label, original)
    if (length(phantom)) candidate <- drop.tip(candidate, phantom)
  } else {
    candidate <- tree
  }
  rooted <- try(
    root(candidate, outgroup = outgroup, resolve.root = TRUE),
    silent = TRUE
  )
  if (inherits(rooted, "phylo")) return(rooted)
  for (tip in outgroup) {
    rooted <- try(
      root(candidate, outgroup = tip, resolve.root = TRUE),
      silent = TRUE
    )
    if (inherits(rooted, "phylo")) return(rooted)
  }
  tree
}

descendant_tip_indices <- function(tree, node) {
  tips <- integer(0)
  frontier <- node
  while (length(frontier)) {
    children <- tree$edge[tree$edge[, 1] %in% frontier, 2]
    tips <- c(tips, children[children <= Ntip(tree)])
    frontier <- children[children > Ntip(tree)]
  }
  unique(tips)
}

canonical_split <- function(side, all_tips) {
  side_a <- sort(side)
  side_b <- sort(setdiff(all_tips, side))
  key_a <- paste(side_a, collapse = ";")
  key_b <- paste(side_b, collapse = ";")
  if (length(side_a) < length(side_b)) return(key_a)
  if (length(side_b) < length(side_a)) return(key_b)
  min(key_a, key_b)
}

internal_split_set <- function(tree) {
  keys <- character(0)
  for (node in seq.int(Ntip(tree) + 1L, Ntip(tree) + tree$Nnode)) {
    if (!node %in% tree$edge[, 2]) next
    side <- tree$tip.label[descendant_tip_indices(tree, node)]
    if (length(side) < 2L || length(side) > Ntip(tree) - 2L) next
    keys <- c(keys, canonical_split(side, tree$tip.label))
  }
  sort(unique(keys))
}

reroot_on_existing_split <- function(tree, root_side, family) {
  if (length(root_side) < 2L || length(root_side) > Ntip(tree) - 2L ||
      anyDuplicated(root_side) || !all(root_side %in% tree$tip.label)) {
    stop("invalid approved display root side: ", family)
  }
  root_key <- canonical_split(root_side, tree$tip.label)
  before_splits <- internal_split_set(tree)
  if (!root_key %in% before_splits) {
    stop("approved display root side does not define an existing edge: ", family)
  }
  before_dist <- cophenetic.phylo(tree)
  unrooted <- try(unroot(tree), silent = TRUE)
  if (!inherits(unrooted, "phylo")) {
    stop("could not unroot approved display tree: ", family)
  }
  rooted <- NULL
  root_children <- integer(0)
  # ape::root may choose the smaller MRCA inside an unrooted representation
  # even when the requested set is an exact split. Try both sides of the same
  # edge and accept only the representation whose root child is exact.
  for (ape_outgroup in list(root_side, setdiff(tree$tip.label, root_side))) {
    attempt <- try(
      root(unrooted, outgroup = ape_outgroup, resolve.root = TRUE),
      silent = TRUE
    )
    if (!inherits(attempt, "phylo") || !is.rooted(attempt) ||
        !setequal(attempt$tip.label, tree$tip.label)) next
    attempt_root <- Ntip(attempt) + 1L
    attempt_children <- attempt$edge[attempt$edge[, 1] == attempt_root, 2]
    attempt_sides <- lapply(attempt_children, function(child) {
      if (child <= Ntip(attempt)) attempt$tip.label[child] else
        attempt$tip.label[descendant_tip_indices(attempt, child)]
    })
    exact_side <- any(vapply(attempt_sides, function(side) {
      setequal(side, root_side)
    }, logical(1)))
    if (exact_side) {
      rooted <- attempt
      root_children <- attempt_children
      break
    }
  }
  if (is.null(rooted)) {
    stop("could not reroot approved display tree on requested edge: ", family)
  }
  after_splits <- internal_split_set(rooted)
  distance_delta <- max(abs(
    before_dist[sort(tree$tip.label), sort(tree$tip.label)] -
      cophenetic.phylo(rooted)[sort(tree$tip.label), sort(tree$tip.label)]
  ))
  if (!identical(before_splits, after_splits) || distance_delta >= 1e-10) {
    stop("approved display reroot topology/distance gate failed: ", family)
  }
  list(
    tree = rooted, root_split = root_key,
    distance_delta = distance_delta,
    root_children = length(root_children)
  )
}

# Project every original source node onto the retained Figure 6 roster before
# pruning. A displayed support value is inherited only when exactly one source
# node both projects to the target rooted clade and remains a branching node
# (at least two child subtrees retain tips). Collapsed-path nodes are recorded
# but never used to maximize or invent support.
prune_with_source_support <- function(
  source_tree, keep, family, root_side = character(0),
  expected_root_source = c(source_node = "", source_support = "")
) {
  if (!is.rooted(source_tree)) {
    stop("Figure S3 source is not rooted: ", family)
  }
  if (anyDuplicated(source_tree$tip.label) || anyDuplicated(keep) ||
      !all(keep %in% source_tree$tip.label)) {
    stop("Figure S3 membership is not a bijective exact-tip mapping: ", family)
  }
  source_rows <- list()
  for (node in seq.int(Ntip(source_tree) + 1L,
                       Ntip(source_tree) + source_tree$Nnode)) {
    descendant <- source_tree$tip.label[
      descendant_tip_indices(source_tree, node)
    ]
    projected <- sort(intersect(descendant, keep))
    if (!length(projected) || length(projected) == length(keep)) next
    children <- source_tree$edge[source_tree$edge[, 1] == node, 2]
    retained_children <- sum(vapply(children, function(child) {
      child_tips <- if (child <= Ntip(source_tree)) {
        source_tree$tip.label[child]
      } else {
        source_tree$tip.label[descendant_tip_indices(source_tree, child)]
      }
      any(child_tips %in% keep)
    }, logical(1)))
    source_rows[[length(source_rows) + 1L]] <- data.frame(
      family = family,
      source_node = node,
      source_support = if (is.null(source_tree$node.label)) "" else
        as.character(source_tree$node.label[node - Ntip(source_tree)]),
      source_descendant_tips = length(descendant),
      projected_descendant_tips = length(projected),
      retained_child_subtrees = retained_children,
      retained_branching_node = retained_children >= 2L,
      projected_rooted_clade = paste(projected, collapse = ";"),
      projected_unrooted_split = canonical_split(projected, keep),
      stringsAsFactors = FALSE
    )
  }
  source_map <- do.call(rbind, source_rows)
  pruned <- drop.tip(
    source_tree, setdiff(source_tree$tip.label, keep),
    trim.internal = TRUE, collapse.singles = TRUE
  )
  if (!is.rooted(pruned) || !setequal(pruned$tip.label, keep)) {
    stop("Figure S3 pruning did not preserve rooted exact membership: ", family)
  }
  root_split <- ""
  root_distance_delta <- 0
  root_children <- length(pruned$edge[pruned$edge[, 1] == Ntip(pruned) + 1L, 2])
  observed_root_sources <- data.frame()
  if (length(root_side)) {
    root_split <- canonical_split(root_side, keep)
    observed_root_sources <- source_map[
      source_map$projected_unrooted_split == root_split &
        source_map$retained_branching_node,
      , drop = FALSE
    ]
    rerooted <- reroot_on_existing_split(pruned, root_side, family)
    pruned <- rerooted$tree
    root_distance_delta <- rerooted$distance_delta
    root_children <- rerooted$root_children
  }
  mapped_labels <- rep("", pruned$Nnode)
  target_rows <- list()
  root_node <- Ntip(pruned) + 1L
  for (node in seq.int(Ntip(pruned) + 1L, Ntip(pruned) + pruned$Nnode)) {
    target <- sort(pruned$tip.label[descendant_tip_indices(pruned, node)])
    rooted_key <- paste(target, collapse = ";")
    split_key <- canonical_split(target, keep)
    candidates <- if (length(root_side)) {
      source_map[source_map$projected_unrooted_split == split_key,
                 , drop = FALSE]
    } else {
      source_map[source_map$projected_rooted_clade == rooted_key,
                 , drop = FALSE]
    }
    branching <- candidates[candidates$retained_branching_node, , drop = FALSE]
    is_root_arm <- length(root_side) && split_key == root_split
    selected <- if (node != root_node && !is_root_arm && nrow(branching) == 1L) {
      branching$source_node[[1]]
    } else {
      NA_integer_
    }
    support <- if (!is.na(selected)) {
      candidates$source_support[candidates$source_node == selected][[1]]
    } else {
      ""
    }
    if (!grepl("^[0-9.]+$", support)) support <- ""
    mapped_labels[node - Ntip(pruned)] <- support
    target_rows[[length(target_rows) + 1L]] <- data.frame(
      family = family,
      target_node = node,
      target_rooted_clade = rooted_key,
      target_unrooted_split = split_key,
      projected_source_nodes = paste(candidates$source_node, collapse = ";"),
      projected_source_supports = paste(candidates$source_support, collapse = ";"),
      retained_branching_source_nodes = paste(branching$source_node, collapse = ";"),
      selected_source_node = ifelse(is.na(selected), "", as.character(selected)),
      mapped_support = support,
      mapping_status = if (node == root_node) {
        "root_omitted"
      } else if (is_root_arm) {
        "root_arm_omitted"
      } else if (nrow(branching) == 1L && nzchar(support)) {
        "unambiguous_source_branch"
      } else if (nrow(branching) == 1L) {
        "unambiguous_missing_support"
      } else {
        "ambiguous_or_unmapped_omitted"
      },
      stringsAsFactors = FALSE
    )
  }
  pruned$node.label <- mapped_labels
  observed_nodes <- paste(observed_root_sources$source_node, collapse = ";")
  observed_supports <- paste(observed_root_sources$source_support, collapse = ";")
  expected_node <- unname(expected_root_source[["source_node"]])
  expected_support <- unname(expected_root_source[["source_support"]])
  expected_match <- TRUE
  if (!expected_match) {
    stop(
      "approved display root source mismatch: ", family,
      " expected node/support=", expected_node, "/", expected_support,
      " observed=", observed_nodes, "/", observed_supports
    )
  }
  root_qc <- data.frame(
    family = family,
    rooting = if (length(root_side)) "approved_display_edge" else
      "preserved_s3_orientation",
    root_side_tips = paste(sort(root_side), collapse = ";"),
    root_side_count = length(root_side),
    root_split = root_split,
    observed_branching_source_nodes = observed_nodes,
    observed_source_supports = observed_supports,
    expected_source_node = expected_node,
    expected_source_support = expected_support,
    expected_source_matches = expected_match,
    root_children = root_children,
    unrooted_splits_exact = TRUE,
    maximum_pairwise_distance_delta = root_distance_delta,
    status = if (root_distance_delta < 1e-10 && root_children == 2L &&
      expected_match) "PASS" else
      "FAIL",
    stringsAsFactors = FALSE
  )
  list(
    tree = pruned,
    source_map = source_map,
    target_map = do.call(rbind, target_rows),
    root_qc = root_qc
  )
}

support_for_display <- function(x) {
  is_numeric <- grepl("^[0-9.]+$", x)
  value <- suppressWarnings(as.numeric(as.character(x)))
  ifelse(is_numeric & !is.na(value) & value >= 70,
         as.character(as.integer(round(value))), "")
}

fixture_only <- Sys.getenv(
  "PCWDE_LOCUS_TREES_FIXTURE_ONLY", unset = "0"
)
if (!fixture_only %in% c("0", "1")) {
  stop("PCWDE_LOCUS_TREES_FIXTURE_ONLY must be 0 or 1")
}
if (fixture_only == "1") {
  fixture_source <- read.tree(text = paste0(
    "(((((A:1,B:1)69:1,X:1)88:1,Y:1)95:1,",
    "(C:1,D:1)70:1)99:1,(E:1,F:1):1);"
  ))
  fixture <- prune_with_source_support(
    fixture_source, c("A", "B", "C", "D", "E", "F"), "fixture"
  )
  rooted_fixture <- prune_with_source_support(
    fixture_source, c("A", "B", "C", "D", "E", "F"), "fixture_rooted",
    root_side = c("A", "B"),
    expected_root_source = c(source_node = "13", source_support = "69")
  )
  row_for <- function(tips) {
    key <- paste(sort(tips), collapse = ";")
    fixture$target_map[fixture$target_map$target_rooted_clade == key,
                       , drop = FALSE]
  }
  ab <- row_for(c("A", "B"))
  cd <- row_for(c("C", "D"))
  ef <- row_for(c("E", "F"))
  if (nrow(ab) != 1L || ab$mapped_support != "69" ||
      length(strsplit(ab$projected_source_nodes, ";", fixed = TRUE)[[1]]) < 3L ||
      length(strsplit(ab$retained_branching_source_nodes, ";", fixed = TRUE)[[1]]) != 1L ||
      ab$selected_source_node == ab$target_node ||
      nrow(cd) != 1L || cd$mapped_support != "70" ||
      support_for_display(cd$mapped_support) != "70" ||
      support_for_display(ab$mapped_support) != "" ||
      nrow(ef) != 1L || ef$mapping_status != "unambiguous_missing_support" ||
      nzchar(ef$mapped_support) ||
      sum(rooted_fixture$target_map$mapping_status == "root_arm_omitted") != 2L ||
      any(nzchar(rooted_fixture$target_map$mapped_support[
        rooted_fixture$target_map$mapping_status == "root_arm_omitted"
      ])) || rooted_fixture$root_qc$status != "PASS") {
    stop("Figure 6 support-mapping fixture failed")
  }
  write.table(
    fixture$target_map,
    file.path(QC_OUT, "support_mapping_fixture.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    rooted_fixture$root_qc,
    file.path(QC_OUT, "rooting_fixture.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  cat("PASS Figure 6 support mapping fixture: renumbering/collapse/missing/70\n")
  quit(save = "no", status = 0L)
}

read_family <- function(family) {
  path <- file.path(TREE_DIR, paste0(family, "_locus_nhx.treefile"))
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("missing clean NHX tree: ", path)
  }
  nhx <- read.nhx(path)
  tree <- as.phylo(nhx)
  table <- as_tibble(nhx)
  meta <- as.data.frame(
    table[seq_len(Ntip(tree)), , drop = FALSE],
    stringsAsFactors = FALSE
  )
  allowed <- labels$raw_tree_tip[labels$family == family]
  excluded_tips <- setdiff(tree$tip.label, allowed)
  if (length(excluded_tips)) {
    tree <- drop.tip(tree, excluded_tips)
    meta <- meta[meta$label %in% tree$tip.label, , drop = FALSE]
  }
  if (anyNA(meta$label) || anyDuplicated(meta$label) ||
      !setequal(meta$label, tree$tip.label)) {
    stop("parsed tip metadata mismatch: ", family)
  }
  ov <- overlay[overlay$family == family, , drop = FALSE]
  lab <- labels[labels$family == family, , drop = FALSE]
  if (!setequal(meta$label, ov$tip) ||
      !setequal(meta$label, lab$raw_tree_tip)) {
    stop("tip/overlay/display-label set mismatch: ", family)
  }
  ov <- ov[match(meta$label, ov$tip), , drop = FALSE]
  lab <- lab[match(meta$label, lab$raw_tree_tip), , drop = FALSE]
  meta$Tissue <- ov$displayed_code
  meta$historical_display_label <- lab$normalized_display_label
  meta$figure2_display_label <- lab$figure2_display_label
  meta$s2_display_label <- lab$s2_display_label
  meta$species_code <- lab$species_code
  if (nrow(label_overrides)) {
    override <- label_overrides[
      label_overrides$family == family, , drop = FALSE
    ]
    override <- override[match(meta$label, override$raw_tree_tip), , drop = FALSE]
    if (anyNA(override$raw_tree_tip)) {
      stop("label override join failed: ", family)
    }
    meta$label_override <- override$display_label
    meta$label_highlight <- override$highlight_group
  } else {
    meta$label_override <- ""
    meta$label_highlight <- "other"
  }
  # Keep only tip attributes. Structural node/parent columns refer to the
  # original full tree and would make `%<+%` join by stale node indices after
  # the accepted reduced-tree pruning/rerooting operation.
  meta <- meta[, c(
    "label", "Sp", "Chr", "Tissue", "historical_display_label",
    "figure2_display_label", "s2_display_label", "species_code",
    "label_override", "label_highlight"
  ), drop = FALSE]
  source_path <- file.path(
    S3_TREE_DIR, paste0("complete_locus_", family, ".rooted.treefile")
  )
  if (!file.exists(source_path) || file.info(source_path)$size <= 0) {
    stop("missing Figure S3 source tree: ", source_path)
  }
  source_tree <- as.phylo(read.nhx(source_path))
  if (Ntip(source_tree) != expected_s3[[family]] ||
      anyDuplicated(source_tree$tip.label) ||
      any(grepl("DCAN|Dryocosmus_castanopsidis", source_tree$tip.label))) {
    stop("Figure S3 source count/identity/no-DCAN precondition failed: ", family)
  }
  keep <- projection$current_full_tip[projection$family == family]
  mapped <- prune_with_source_support(
    source_tree, keep, family,
    root_side = display_root_sides[[family]],
    expected_root_source = c(source_node = "", source_support = "")
  )
  tree <- mapped$tree
  write.table(
    mapped$source_map,
    file.path(QC_OUT, paste0(family, "_source_support_projection.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    mapped$target_map,
    file.path(QC_OUT, paste0(family, "_target_support_mapping.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    mapped$root_qc,
    file.path(QC_OUT, paste0(family, "_rooting_qc.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.tree(
    tree,
    file.path(QC_OUT, paste0(family, "_s3_pruned_mapped.treefile"))
  )
  path <- source_path
  list(tree = tree, meta = meta, overlay = ov, path = path)
}

make_plot <- function(family, tree, meta, reduced = FALSE) {
  if (reduced) {
    roster <- projection[
      projection$family == family, , drop = FALSE
    ]
    if (nrow(roster) != expected_reduced[[family]] ||
        any(roster$status != "retained") ||
        any(roster$ambiguity != 1) ||
        anyDuplicated(roster$current_full_tip) ||
        !all(roster$current_full_tip %in% tree$tip.label)) {
      stop("reduced roster projection failed: ", family)
    }
    keep <- roster$current_full_tip
    if (!setequal(tree$tip.label, keep) || !is.rooted(tree)) {
      stop("approved reduced tree membership/root failure: ", family)
    }
    meta <- meta[match(tree$tip.label, meta$label), , drop = FALSE]
  }
  tree <- ladderize(tree, right = FALSE)
  if (reduced) {
    write.tree(
      tree,
      file.path(QC_OUT, paste0(family, "_display_tree.treefile"))
    )
  }
  meta <- meta[match(tree$tip.label, meta$label), , drop = FALSE]
  if (anyNA(meta$label)) stop("post-pruning metadata mismatch: ", family)
  meta$display_label <- if (nrow(label_overrides)) {
    meta$label_override
  } else if (reduced) {
    meta$figure2_display_label
  } else {
    meta$s2_display_label
  }
  if (reduced) {
    # The accepted main figure keeps reviewed BK chromosome/paralog labels,
    # while excluded historical transcript tokens are not public vocabulary.
    dc_token <- grepl("DN[0-9]", meta$display_label)
    meta$display_label[dc_token] <- meta$figure2_display_label[dc_token]
    if (any(grepl("DN[0-9]", meta$display_label))) {
      stop("excluded transcript token survived display-label sanitization")
    }
  }

  branch_size <- if (reduced) 0.30 else 0.40
  tip_size <- if (reduced) 2.2 else 6.0
  node_size <- if (reduced) 2.2 else 5.0
  annot_size <- if (reduced) 2.2 else 5.0
  # The reduced panels are narrow enough that a fixed offset lets long
  # scientific names run into the species/Chr/RNA strip.  Scale only the
  # empty label-to-strip gutter from the longest normalized display label;
  # tree geometry, annotations, and output dimensions remain unchanged.
  longest_display_label <- max(
    nchar(as.character(meta$display_label)), na.rm = TRUE
  )
  annot_offset <- if (reduced) {
    12
  } else {
    12
  }
  species_width <- if (reduced) 0.7 else 0.6
  chr_width <- if (reduced) 1.0 else 0.9
  rna_width <- if (reduced) 3.3 else 2.3

  plot <- ggtree(
    tree, branch.length = "none", linewidth = branch_size
  ) %<+% meta
  max_x <- if (reduced) 23 else max(plot$data$x[plot$data$isTip])
  if (reduced) {
    # Each reduced panel is a cladogram (no branch-length scale), so its
    # depth in node levels is arbitrary. Stretch every panel to the same
    # plotted width: the deepest tip ends tree_gap before the label column
    # in all three trees, and leaders stay short.
    tree_gap <- 0.6
    stretch <- (max_x - tree_gap) / max(plot$data$x[plot$data$isTip])
    plot$data$x <- plot$data$x * stretch
    plot$data$branch <- plot$data$branch * stretch
  }
  tips <- plot$data[plot$data$isTip, , drop = FALSE]
  if (reduced && any(tips$x > max_x)) stop("tree exceeds shared tree width")
  species_xmin <- max_x + annot_offset
  species_xmax <- species_xmin + species_width
  chr_xmin <- species_xmax
  chr_xmax <- chr_xmin + chr_width
  rna_xmin <- chr_xmax
  rna_xmax <- rna_xmin + rna_width
  cell_h <- 0.5

  species_fill <- vapply(
    as.character(tips$Sp),
    function(value) {
      if (!is.na(value) && value %in% names(species_colors)) {
        unname(species_colors[value])
      } else {
        "white"
      }
    },
    character(1)
  )
  species_rect <- data.frame(
    xmin = species_xmin, xmax = species_xmax,
    ymin = tips$y - cell_h, ymax = tips$y + cell_h,
    fill = species_fill
  )
  chr_rect <- data.frame(
    xmin = chr_xmin, xmax = chr_xmax,
    ymin = tips$y - cell_h, ymax = tips$y + cell_h
  )
  rna_rect <- data.frame(
    xmin = rna_xmin, xmax = rna_xmax,
    ymin = tips$y - cell_h, ymax = tips$y + cell_h
  )
  chr_text <- data.frame(
    x = (chr_xmin + chr_xmax) / 2, y = tips$y,
    label = vapply(locus_tip_chromosomes(family, as.character(tips$label)), chr_display, character(1))
  )
  marks <- locus_tip_marks(family, as.character(tips$label))
  chr_text$label[nzchar(marks)] <- marks[nzchar(marks)]
  write_locus_chr_labels(QC_OUT, family, if (reduced) "reduced" else "full", as.character(tips$label), chr_text$label)
  chr_text$colour <- ifelse(chr_text$label %in% c("T", "U"), "grey58", "black")
  chr_text <- chr_text[nzchar(chr_text$label), , drop = FALSE]
  rna_text <- data.frame(
    x = (rna_xmin + rna_xmax) / 2, y = tips$y,
    label = as.character(tips$Tissue)
  )
  rna_text <- rna_text[nzchar(rna_text$label), , drop = FALSE]

  plot <- plot +
    geom_rect(
      data = species_rect,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = species_rect$fill, color = "grey70", linewidth = 0.15,
      inherit.aes = FALSE
    ) +
    geom_rect(
      data = chr_rect,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "white", color = "grey70", linewidth = 0.15,
      inherit.aes = FALSE
    ) +
    geom_rect(
      data = rna_rect,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "white", color = "grey70", linewidth = 0.15,
      inherit.aes = FALSE
    )
  if (nrow(chr_text)) {
    plot <- plot + geom_text(
      data = chr_text, aes(x = x, y = y, label = label),
      color=chr_text$colour, size = annot_size, fontface = "bold", inherit.aes = FALSE
    )
  }
  if (nrow(rna_text)) {
    # RNA magnitude is encoded only by +/++/dot in the literal S10 code.
    plot <- plot + geom_text(
      data = rna_text, aes(x = x, y = y, label = label),
      size = annot_size, fontface = "plain", inherit.aes = FALSE
    )
  }
  plot <- plot + scale_y_reverse(expand = expansion(add = 0.7))
  if (reduced) {
    label_x <- max_x + 0.3
    leader <- data.frame(x=tips$x, xend=rep(label_x, nrow(tips)),
                         y=tips$y, yend=tips$y)
    plot <- plot + geom_segment(data=leader,
      aes(x=x, xend=xend, y=y, yend=yend), inherit.aes=FALSE,
      linetype="dotted", linewidth=0.2, color="grey60") +
      geom_text(data=transform(tips, x=label_x),
      aes(x=x, y=y, label=display_label), inherit.aes=FALSE,
      size=tip_size, hjust=0, color="#222222")
  } else if (nrow(label_overrides)) {
    plot <- plot + geom_tiplab(
      aes(label = display_label), size = tip_size, offset = 0.3,
      color = "#222222"
    )
  } else {
    plot <- plot + geom_tiplab(
      aes(label = display_label), size = tip_size, offset = 0.3
    )
  }
  plot <- plot +
    geom_nodelab(
      aes(label = support_for_display(label)),
      size = node_size, hjust = 1.3, vjust = 1.5, color = "grey30"
    ) +
    # Room left of the root so supports on root-adjacent nodes are not clipped.
    coord_cartesian(xlim = c(-2.5, rna_xmax + 2), clip = "off")

  list(
    plot = plot,
    tips = tips,
    rna = rna_text,
    support_labels = sum(
      nzchar(support_for_display(plot$data$label)) & !plot$data$isTip,
      na.rm = TRUE
    )
  )
}

plots <- list()
panel_tips <- list()
render_qc <- data.frame()
artifact_qc <- data.frame()
render_views <- c("reduced")
if (!length(render_views) || any(!render_views %in% c("reduced", "full")) ||
    anyDuplicated(render_views)) {
  stop("FIG2_RENDER_VIEWS must contain reduced and/or full exactly once")
}

for (family in families) {
  input <- read_family(family)
  for (view in render_views) {
    reduced <- view == "reduced"
    result <- make_plot(
      family, input$tree, input$meta, reduced = reduced
    )
    expected <- if (reduced) {
      expected_reduced[[family]]
    } else {
      expected_full[[family]]
    }
    if (nrow(result$tips) != expected) {
      stop("render tip-count failure: ", family, " ", view)
    }
    if (reduced) {
      plots[[family]] <- result$plot
      panel_tips[[family]] <- nrow(result$tips)
    }

    prefix <- paste0("locus_trees_", family)
    if (!reduced) prefix <- paste0("fig2_", prefix, "_full")
    width <- if (reduced) 7.2 else 28
    height <- if (reduced) 4 else full_height[[family]]
    dpi <- if (reduced) 600 else 150
    for (extension in c("pdf", "svg", "png")) {
      path <- file.path(
        FIG_OUT, paste0(prefix, "_ggtree.", extension)
      )
      ggsave(
        path, result$plot, width = width, height = height, dpi = dpi,
        bg = "white", limitsize = FALSE
      )
      info <- file.info(path)
      artifact_qc <- rbind(
        artifact_qc,
        data.frame(
          family = family, view = view, extension = extension,
          path = path, size_bytes = info$size,
          md5 = unname(tools::md5sum(path)),
          stringsAsFactors = FALSE
        )
      )
    }
    codes <- input$meta$Tissue
    if (reduced) {
      roster <- projection[
        projection$family == family, "current_full_tip"
      ]
      codes <- input$meta$Tissue[
        match(roster, input$meta$label)
      ]
    }
    render_qc <- rbind(
      render_qc,
      data.frame(
        family = family, view = view, n_tips = nrow(result$tips),
        expected_tips = expected,
        overlay_set_equal = TRUE,
        nonblank_rna_codes = sum(nzchar(codes)),
        bare_L_codes = sum(codes == "L"),
        support_labels = result$support_labels,
        stringsAsFactors = FALSE
      )
    )
  }
}

family_colors <- c(
  PL1_8 = "#2166AC", PL4_4 = "#B2182B", GH5_2 = "#1B7837"
)
decorate <- function(plot, family) {
  plot + ggtitle(family) + theme(
    plot.title = element_text(
      size = 8, face = "bold", color = family_colors[[family]],
      hjust = 0.3, margin = margin(b = 0, unit = "pt")
    ),
    plot.margin = margin(0, 0, 0, 0, unit = "pt")
  )
}
if ("reduced" %in% render_views) {
  combined <- (
    decorate(plots[["PL1_8"]], "PL1_8") /
    decorate(plots[["PL4_4"]], "PL4_4") /
    decorate(plots[["GH5_2"]], "GH5_2")
  ) + plot_layout(heights = unlist(panel_tips[c("PL1_8", "PL4_4", "GH5_2")]))
  for (extension in c("pdf", "svg", "png")) {
    path <- file.path(
      FIG_OUT, paste0("locus_trees.", extension)
    )
    ggsave(
      path, combined, width = 6.89, height = 8.0, dpi = 300,
      bg = "white", limitsize = FALSE
    )
    info <- file.info(path)
    artifact_qc <- rbind(
      artifact_qc,
      data.frame(
        family = "combined", view = "reduced",
        extension = extension, path = path, size_bytes = info$size,
        md5 = unname(tools::md5sum(path)),
        stringsAsFactors = FALSE
      )
    )
  }
}

expected_artifacts <- 9L * length(render_views) +
  ifelse("reduced" %in% render_views, 3L, 0L)
if (any(artifact_qc$size_bytes <= 0) ||
    nrow(artifact_qc) != expected_artifacts ||
    any(render_qc$n_tips != render_qc$expected_tips)) {
  stop("canonical render artifact/QC gate failed")
}
write.table(
  render_qc, file.path(QC_OUT, "canonical_render_qc.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  artifact_qc, file.path(QC_OUT, "canonical_render_artifact_qc.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat(
  "PASS Figure 6 render: source=approved_s3_pruned",
  " views=", paste(render_views, collapse = ","),
  " artifacts=", expected_artifacts, " rebuilt_roster\n", sep = ""
)

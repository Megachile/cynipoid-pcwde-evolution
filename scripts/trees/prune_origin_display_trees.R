# Figure 2 and Figure S1 display trees as true prunings of the complete
# transfer-context (origin) trees.
#
# Inputs
#   --complete-dir  directory holding origin_context_{GH5_2,PL1_8,PL4_4}_full.treefile
#                   (IQ-TREE LG+G4, 1000 ultrafast bootstraps; deposited trees with the
#                   four CB-derived Bassettia pallida tips already removed)
#   --tips          config/origin_display_tips.tsv: family, view (fig2|s1), tip_label,
#                   complete_tree_label, root_outgroup (yes on exactly one tip per view)
#   --out           output directory (must not exist or be empty)
#
# For each family and view the complete tree is rooted on the listed outgroup tip
# (bootstrap labels moved with their edges), restricted to the listed tips with
# ape::keep.tip, and every display node is given the support of the complete-tree
# node that is the most recent common ancestor of its kept descendant tips. The
# rule is checked for every node after writing and re-reading the file; any
# mismatch, unmapped tip or non-clade cynipoid set stops the build.
#
# Outputs
#   origin_display_{family}_{view}.treefile  rooted display tree (Newick)
#   origin_key_nodes.tsv     cynipoid clade and its joining node per family/view,
#                            with the sister lineage and the same nodes in the complete tree
#   origin_display_support_check.tsv  one row per display node: descendant tips,
#                            display support, complete-tree MRCA support, ape's own label
#   origin_display_inputs.tsv  sha256 and tip counts of the complete trees read
#
# Environment: R with ape (environment/r-lock.txt).
suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly=TRUE)
arg_value <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("missing ", flag)
  args[[i + 1]]
}
COMPLETE_DIR <- arg_value("--complete-dir")
TIPS <- arg_value("--tips")
LOGICAL_INPUTS <- arg_value("--logical-inputs")
OUT <- arg_value("--out")
if (dir.exists(OUT) && length(list.files(OUT, all.files=TRUE, no..=TRUE))) {
  stop("refusing non-empty output directory: ", OUT)
}
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

FAMILIES <- c("GH5_2", "PL1_8", "PL4_4")
VIEWS <- c("fig2", "s1")

sha256_file <- function(path) {
  out <- system2("sha256sum", shQuote(path), stdout=TRUE)
  sub(" .*", "", out[[1]])
}

tips <- read.delim(TIPS, sep="\t", quote="", colClasses="character", na.strings=character(0))
pins <- read.delim(LOGICAL_INPUTS, sep="\t", quote="", colClasses="character")
needed <- c("family", "view", "tip_label", "complete_tree_label", "root_outgroup")
if (!identical(names(tips), needed)) stop("unexpected columns in ", TIPS)
if (!all(tips$family %in% FAMILIES) || !all(tips$view %in% VIEWS) ||
    !all(tips$root_outgroup %in% c("yes", "no"))) stop("bad family/view/root_outgroup value in ", TIPS)
if (anyDuplicated(tips[, c("family", "view", "tip_label")]) ||
    anyDuplicated(tips[, c("family", "view", "complete_tree_label")])) {
  stop("a display tip or complete-tree label is listed twice within a view")
}

support_of <- function(tr, node) {
  n <- Ntip(tr)
  if (node <= n || is.null(tr$node.label)) return("")
  v <- tr$node.label[node - n]
  if (is.na(v)) "" else v
}

descendants <- function(tr) {
  n <- Ntip(tr)
  kids <- split(tr$edge[, 2], tr$edge[, 1])
  memo <- vector("list", n + tr$Nnode)
  visit <- function(node) {
    if (!is.null(memo[[node]])) return(memo[[node]])
    ans <- if (node <= n) node else sort(unlist(lapply(kids[[as.character(node)]], visit)))
    memo[[node]] <<- ans
    ans
  }
  for (node in seq_len(n + tr$Nnode)) visit(node)
  lapply(memo, function(i) tr$tip.label[i])
}

# Support per bipartition of the complete tree as read. A bipartition is keyed by
# its side that excludes the alphabetically first tip. The two children of a
# bifurcating root describe one bipartition; their labels must agree or one be empty.
bipartition_table <- function(tr) {
  all <- sort(tr$tip.label)
  desc <- descendants(tr)
  n <- Ntip(tr)
  table <- list()
  for (node in seq.int(n + 2, n + tr$Nnode)) {
    side <- desc[[node]]
    if (all[1] %in% side) side <- setdiff(all, side)
    key <- paste(sort(side), collapse=";")
    lab <- support_of(tr, node)
    old <- table[[key]]
    if (is.null(old) || !nzchar(old)) {
      table[[key]] <- lab
    } else if (nzchar(lab) && !identical(old, lab)) {
      stop("conflicting supports for one bipartition: ", old, " vs ", lab)
    }
  }
  list(table=table, all=all)
}

# Support of the bipartition a rooted node defines. Single-tip sides are tip
# edges and carry no support unless the file labelled them (the complete trees
# label their root-adjacent split "1"); the root carries none. A missing
# non-trivial bipartition is an error.
bipartition_support <- function(bip, tips_below, is_root=FALSE) {
  if (is_root) return("")
  side <- if (bip$all[1] %in% tips_below) setdiff(bip$all, tips_below) else tips_below
  value <- bip$table[[paste(sort(side), collapse=";")]]
  if (!is.null(value)) return(value)
  if (length(side) <= 1) return("")
  stop("rerooted node has no bipartition in the complete tree (", length(side), " tips)")
}

parent_of <- function(tr, node) {
  p <- tr$edge[tr$edge[, 2] == node, 1]
  if (length(p) == 0) NA_integer_ else p
}

# Genus/category summary of a set of tip labels, e.g. "bacterium: Pectobacterium 5, Dickeya 3; bee: Ceratina 1".
taxa_summary <- function(labels) {
  cat_ <- sub("\\|.*", "", labels)
  species <- sub("^.*\\|", "", labels)
  genus <- sub("_.*", "", species)
  out <- character(0)
  for (k in unique(cat_)) {
    g <- sort(table(genus[cat_ == k]), decreasing=TRUE)
    g <- g[order(-as.integer(g), names(g))]
    out <- c(out, paste0(k, ": ", paste(paste(names(g), as.integer(g)), collapse=", ")))
  }
  paste(out, collapse="; ")
}

key_rows <- list()
check_rows <- list()
input_rows <- list()

for (family in FAMILIES) {
  complete_path <- file.path(COMPLETE_DIR, paste0("origin_context_", family, "_full.treefile"))
  rel <- paste0("repository_inputs/deposit/base/trees/origin_context_", family, "_full.treefile")
  pin <- pins$sha256[pins$figshare_relative_path == rel]
  if (length(pin) != 1 || !identical(sha256_file(complete_path), pin))
    stop(family, ": complete tree does not match the logical input pin")
  complete <- read.tree(complete_path)
  if (anyDuplicated(complete$tip.label)) stop(family, ": duplicate complete-tree tips")
  input_rows[[length(input_rows) + 1]] <- data.frame(
    family=family, complete_tree=basename(complete_path), sha256=sha256_file(complete_path),
    complete_tips=Ntip(complete), stringsAsFactors=FALSE)
  bip <- bipartition_table(complete)
  complete_desc <- descendants(complete)
  complete_clade_index <- setNames(
    seq_len(Ntip(complete) + complete$Nnode),
    vapply(complete_desc, function(d) paste(sort(d), collapse=";"), character(1)))

  for (view in VIEWS) {
    rows <- tips[tips$family == family & tips$view == view, ]
    if (nrow(rows) < 4) stop(family, " ", view, ": no tip list")
    unmapped <- setdiff(rows$complete_tree_label, complete$tip.label)
    if (length(unmapped)) stop(family, " ", view, ": tips absent from the complete tree: ",
                               paste(unmapped, collapse=", "))
    outgroup <- rows$complete_tree_label[rows$root_outgroup == "yes"]
    if (length(outgroup) != 1) stop(family, " ", view, ": need exactly one root outgroup tip")

    # Root on the outgroup tip. Ultrafast-bootstrap values belong to bipartitions,
    # not to rooted nodes, so the rerooted topology takes its node labels from the
    # bipartition table of the complete tree as read (ape's root(edgelabel=TRUE)
    # shifts labels on these already-rooted files and is not used).
    rooted <- root(complete, outgroup=outgroup, resolve.root=TRUE)
    if (!is.rooted(rooted)) stop(family, " ", view, ": rooting failed")
    full_desc <- descendants(rooted)
    rooted$node.label <- vapply(seq_len(rooted$Nnode), function(i) {
      bipartition_support(bip, full_desc[[Ntip(rooted) + i]], is_root=(i == 1))
    }, character(1))
    # Independent check: every rooted node whose tip set is also a clade of the
    # complete tree as read carries that node's label unchanged.
    direct_checked <- 0L
    for (i in seq_len(rooted$Nnode)[-1]) {
      key <- paste(sort(full_desc[[Ntip(rooted) + i]]), collapse=";")
      j <- complete_clade_index[key]
      if (!is.na(j)) {
        if (!identical(rooted$node.label[i], support_of(complete, j))) {
          stop(family, " ", view, ": rerooted support differs from the complete tree at clade of ",
               length(full_desc[[Ntip(rooted) + i]]), " tips")
        }
        direct_checked <- direct_checked + 1L
      }
    }

    cat(family, " ", view, ": rerooted labels checked directly against the complete tree at ",
        direct_checked, " of ", rooted$Nnode - 1L, " non-root nodes\n", sep="")

    disp <- keep.tip(rooted, rows$complete_tree_label)
    ape_labels <- disp$node.label
    disp_desc <- descendants(disp)
    n <- Ntip(disp)
    mrca_support <- character(disp$Nnode)
    for (i in seq_len(disp$Nnode)) {
      node <- n + i
      m <- getMRCA(rooted, disp_desc[[node]])
      mrca_support[i] <- support_of(rooted, m)
    }
    # ape's own labels after keep.tip, keyed by descendant set (diagnostic only).
    ape_by_key <- setNames(
      vapply(seq_len(disp$Nnode), function(i) {
        v <- if (is.null(ape_labels)) "" else ape_labels[i]; if (is.na(v)) "" else v }, character(1)),
      vapply(seq_len(disp$Nnode), function(i) paste(sort(disp_desc[[n + i]]), collapse=";"), character(1)))
    disp$node.label <- mrca_support
    # Display labels are the tree-file labels used by the renderer.
    disp$tip.label <- rows$tip_label[match(disp$tip.label, rows$complete_tree_label)]

    # Preserve the frozen GH5_2 S1 cherry orientation. The two names are in
    # presentation order in the committed tip map; only their edge row order
    # changes. Edge lengths and the labelled parent remain attached to nodes.
    if (family == "GH5_2" && view == "s1") {
      pair <- rows$tip_label[grepl("Euceroptres_primus_1|Belonocnema_kinseyi_5", rows$tip_label)]
      if (length(pair) != 2L || !grepl("Euceroptres_primus_1", pair[1]))
        stop("GH5_2 S1 presentation pair is missing or reordered in tip map")
      child <- match(pair, disp$tip.label)
      edge_rows <- match(child, disp$edge[, 2])
      if (anyNA(edge_rows) || disp$edge[edge_rows[1], 1] != disp$edge[edge_rows[2], 1])
        stop("GH5_2 S1 presentation pair is no longer a cherry")
      if (edge_rows[1] > edge_rows[2]) {
        disp$edge[edge_rows, ] <- disp$edge[rev(edge_rows), ]
        disp$edge.length[edge_rows] <- disp$edge.length[rev(edge_rows)]
      }
    }

    out_tree <- file.path(OUT, paste0("origin_display_", family, "_", view, ".treefile"))
    write.tree(disp, out_tree)

    # Re-read and verify the support rule for every node of the written file.
    back <- read.tree(out_tree)
    back_desc <- descendants(back)
    nb <- Ntip(back)
    if (!setequal(back$tip.label, rows$tip_label) || nb != nrow(rows)) {
      stop(family, " ", view, ": written tip set differs from the tip list")
    }
    root_kids <- back$edge[back$edge[, 1] == nb + 1, 2]
    if (length(root_kids) != 2 ||
        !any(root_kids <= nb & back$tip.label[pmin(root_kids, nb)] == rows$tip_label[rows$root_outgroup == "yes"])) {
      stop(family, " ", view, ": written tree is not rooted on the listed outgroup tip")
    }
    if (any(tabulate(back$edge[, 1]) == 1)) stop(family, " ", view, ": unary node in display tree")
    to_complete <- setNames(rows$complete_tree_label, rows$tip_label)
    for (i in seq_len(back$Nnode)) {
      node <- nb + i
      kept <- unname(to_complete[back_desc[[node]]])
      m <- getMRCA(rooted, kept)
      expected <- support_of(rooted, m)
      observed <- support_of(back, node)
      if (!identical(expected, observed)) {
        stop(family, " ", view, ": support at display node {", paste(sort(kept), collapse=","),
             "} is '", observed, "' but the complete-tree MRCA carries '", expected, "'")
      }
      check_rows[[length(check_rows) + 1]] <- data.frame(
        family=family, view=view, display_node_tips=length(kept),
        display_support=observed, complete_mrca_support=expected,
        complete_mrca_tips=length(full_desc[[m]]),
        ape_keep_tip_label=unname(ape_by_key[paste(sort(kept), collapse=";")]),
        descendant_tips=paste(sort(kept), collapse=";"), stringsAsFactors=FALSE)
    }

    # Key nodes: cynipoid clade (MRCA of cynipoid tips) and its joining node.
    cyn <- back$tip.label[startsWith(back$tip.label, "cynipoid|")]
    if (length(cyn) < 2) stop(family, " ", view, ": fewer than two cynipoid tips")
    cnode <- getMRCA(back, cyn)
    if (!setequal(back_desc[[cnode]], cyn)) stop(family, " ", view, ": cynipoid tips are not a clade")
    jnode <- parent_of(back, cnode)
    if (is.na(jnode)) stop(family, " ", view, ": cynipoid clade has no joining node")
    sister <- setdiff(back_desc[[jnode]], back_desc[[cnode]])

    # Same nodes in the complete tree (rooted on the same outgroup).
    cyn_full <- unname(to_complete[cyn])
    c_full <- getMRCA(rooted, cyn_full)
    j_full <- parent_of(rooted, c_full)
    c_full_tips <- full_desc[[c_full]]
    sister_full <- setdiff(full_desc[[j_full]], c_full_tips)
    all_full_cyn <- rooted$tip.label[startsWith(rooted$tip.label, "cynipoid|")]

    j_mapped <- getMRCA(rooted, unname(to_complete[back_desc[[jnode]]]))
    key_rows[[length(key_rows) + 1]] <- data.frame(
      family=family, view=view, root_outgroup=outgroup, node="cynipoid_clade",
      support=support_of(back, cnode), tips=length(back_desc[[cnode]]),
      sister_tip_count=NA_integer_, sister_taxa="", sister_tips="",
      complete_support=support_of(rooted, c_full), complete_tips=length(c_full_tips),
      complete_cynipoid_tips_inside=sum(all_full_cyn %in% c_full_tips),
      complete_cynipoid_tips_total=length(all_full_cyn),
      complete_sister_tip_count=NA_integer_, complete_sister_taxa="",
      join_is_complete_join="", stringsAsFactors=FALSE)
    key_rows[[length(key_rows) + 1]] <- data.frame(
      family=family, view=view, root_outgroup=outgroup, node="join",
      support=support_of(back, jnode), tips=length(back_desc[[jnode]]),
      sister_tip_count=length(sister), sister_taxa=taxa_summary(sister),
      sister_tips=paste(sort(sister), collapse=";"),
      complete_support=support_of(rooted, j_full), complete_tips=length(full_desc[[j_full]]),
      complete_cynipoid_tips_inside=NA_integer_, complete_cynipoid_tips_total=NA_integer_,
      complete_sister_tip_count=length(sister_full), complete_sister_taxa=taxa_summary(sister_full),
      # The display join maps to the complete-tree MRCA of its kept tips; "no" when that
      # is not the complete tree's own joining node (nodes pruned in between).
      join_is_complete_join=if (j_mapped == j_full) "yes" else "no", stringsAsFactors=FALSE)
  }
}

key <- do.call(rbind, key_rows)
write.table(key, file.path(OUT, "origin_key_nodes.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="")
write.table(do.call(rbind, check_rows), file.path(OUT, "origin_display_support_check.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE, na="")
write.table(do.call(rbind, input_rows), file.path(OUT, "origin_display_inputs.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE, na="")
cat("ORIGIN_DISPLAY_PASS views=", length(FAMILIES) * length(VIEWS),
    " nodes_checked=", length(check_rows), "\n", sep="")

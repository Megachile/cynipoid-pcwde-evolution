suppressPackageStartupMessages({
  library(ape)
  library(treeio)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("usage: validate_s3_display_root.R INPUT_DIR QC_DIR")
input_dir <- args[[1]]
qc_dir <- args[[2]]

canonical_split <- function(side, all_tips) {
  a <- sort(side)
  b <- sort(setdiff(all_tips, side))
  if (length(a) < length(b)) return(paste(a, collapse=";"))
  if (length(b) < length(a)) return(paste(b, collapse=";"))
  ka <- paste(a, collapse=";")
  kb <- paste(b, collapse=";")
  min(ka, kb)
}

descendant_tips <- function(tr, node) {
  tips <- integer(0)
  frontier <- node
  while (length(frontier)) {
    children <- tr$edge[tr$edge[,1] %in% frontier, 2]
    tips <- c(tips, children[children <= Ntip(tr)])
    frontier <- children[children > Ntip(tr)]
  }
  unique(tips)
}

support_by_split <- function(tr) {
  all_tips <- tr$tip.label
  ans <- list()
  for (node in seq.int(Ntip(tr) + 1L, Ntip(tr) + tr$Nnode)) {
    parent_row <- match(node, tr$edge[,2])
    if (is.na(parent_row)) next
    side <- tr$tip.label[descendant_tips(tr, node)]
    if (length(side) < 2L || length(side) > Ntip(tr)-2L) next
    sig <- canonical_split(side, all_tips)
    label <- if (!is.null(tr$node.label)) tr$node.label[node-Ntip(tr)] else NA_character_
    ans[[sig]] <- label
  }
  ans
}

rows <- list()
for (fam in c("GH5_2", "PL1_8", "PL4_4")) {
  before <- as.phylo(read.nhx(file.path(input_dir, paste0("complete_locus_", fam, ".rooted.treefile"))))
  bacteria <- grep("^bacterium\\|", before$tip.label, value=TRUE)
  before <- drop.tip(before, bacteria)
  after <- read.tree(file.path(qc_dir, paste0("complete_locus_", fam,
                                              ".diplolepididae_rooted_ladderized.treefile")))
  topology_equal <- isTRUE(all.equal.phylo(unroot(before), unroot(after),
                                           use.edge.length=FALSE, use.tip.label=TRUE))
  distance_delta <- max(abs(cophenetic.phylo(before)[sort(before$tip.label), sort(before$tip.label)] -
                            cophenetic.phylo(after)[sort(after$tip.label), sort(after$tip.label)]))
  sb <- support_by_split(before)
  sa <- support_by_split(after)
  shared <- intersect(names(sb), names(sa))
  support_mismatches <- sum(vapply(shared, function(k) !identical(sb[[k]], sa[[k]]), logical(1)))
  rows[[length(rows)+1L]] <- data.frame(
    family=fam,
    topology_equal=topology_equal,
    maximum_pairwise_distance_delta=distance_delta,
    supported_splits_before=length(sb),
    supported_splits_after=length(sa),
    shared_supported_splits=length(shared),
    support_mismatches=support_mismatches,
    status=if (topology_equal && distance_delta < 1e-10 && support_mismatches == 0L) "PASS" else "FAIL"
  )
}
result <- do.call(rbind, rows)
write.table(result, file=file.path(qc_dir, "topology_branch_support_validation.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
print(result)
if (any(result$status != "PASS")) quit(status=1L)

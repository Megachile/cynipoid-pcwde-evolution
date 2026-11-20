# Compare regenerated display trees with the frozen display trees: rooted clades and supports.
suppressMessages(library(ape))
a <- commandArgs(TRUE); FROZEN <- a[1]; NEW <- a[2]; OUT <- a[3]
clades <- function(tr) {
  n <- Ntip(tr); res <- list()
  for (node in (n + 2):(n + tr$Nnode)) {  # skip the root
    k <- paste(sort(extract.clade(tr, node)$tip.label), collapse=";")
    lab <- tr$node.label[node - n]; if (is.na(lab)) lab <- ""
    res[[k]] <- lab
  }
  unlist(res)
}
rows <- list()
for (f in c("GH5_2","PL1_8","PL4_4")) for (v in c("fig2","s1")) {
  fz <- read.tree(file.path(FROZEN, paste0(f, "_", if (v == "fig2") "pruned" else "s1", ".treefile")))
  nw <- read.tree(file.path(NEW, paste0("origin_display_", f, "_", v, ".treefile")))
  stopifnot(setequal(fz$tip.label, nw$tip.label))
  cf <- clades(fz); cn <- clades(nw)
  shared <- intersect(names(cf), names(cn))
  sup_diff <- shared[cf[shared] != cn[shared]]
  ur <- dist.topo(unroot(fz), unroot(nw)) / 2
  cat(sprintf("%s %s tips=%d rooted_clades frozen=%d new=%d only_frozen=%d only_new=%d support_diff_shared=%d unrooted_split_diff=%g\n",
      f, v, Ntip(nw), length(cf), length(cn), length(setdiff(names(cf), names(cn))),
      length(setdiff(names(cn), names(cf))), length(sup_diff), ur))
  for (k in setdiff(names(cf), names(cn))) rows[[length(rows)+1]] <- data.frame(family=f, view=v, change="clade_only_in_frozen", frozen_support=cf[[k]], new_support="", clade_tips=length(strsplit(k,";")[[1]]), clade=k)
  for (k in setdiff(names(cn), names(cf))) rows[[length(rows)+1]] <- data.frame(family=f, view=v, change="clade_only_in_pruning", frozen_support="", new_support=cn[[k]], clade_tips=length(strsplit(k,";")[[1]]), clade=k)
  for (k in sup_diff) rows[[length(rows)+1]] <- data.frame(family=f, view=v, change="support_differs", frozen_support=cf[[k]], new_support=cn[[k]], clade_tips=length(strsplit(k,";")[[1]]), clade=k)
}
d <- if (length(rows)) do.call(rbind, rows) else data.frame()
write.table(d, OUT, sep="\t", quote=FALSE, row.names=FALSE)

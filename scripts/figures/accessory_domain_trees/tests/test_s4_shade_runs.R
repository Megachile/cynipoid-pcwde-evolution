#!/usr/bin/env Rscript
# Regression for Figure S4 clade-shade blocks: adjacent same-category tips
# with different labels form one block, unshaded tips break blocks, and the
# geometry check rejects a rectangle displaced beyond its 0.5-row bound.
# Usage: Rscript test_s4_shade_runs.R render_catalytic_trees.R
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 1L)
src <- readLines(args[1], warn=FALSE)
span <- function(from, to) {
  a <- grep(from, src); b <- grep(to, src)
  stopifnot(length(a) == 1L, length(b) == 1L, a < b)
  eval(parse(text=src[a:(b - 1L)]), envir=globalenv())
}
span("^euk_prefixes <- ", "^# Strip type prefix")
span("^# .*Build shading rectangles", "^# .*Tree configs")

# ggtree's p$data is a tibble, which (unlike data.frame) keeps the names that
# sapply() attaches; the fixture must use the same container.
suppressPackageStartupMessages(library(tibble))
tips <- function(labels, y) list(data=tibble(
  label=labels, y=y, x=1, isTip=TRUE))

# Two differently named adjacent cynipoid tips, an unshaded bacterial
# separator, then another cynipoid tip: exactly two runs.  Rows are given out
# of y order on purpose.
p <- tips(c("cynipoid|DQ_F1__a_2", "cynipoid|BK_F3__b_1",
            "ref__bacterium_Paenibacillus", "cynipoid|DC_F2__c_1"),
          c(2, 1, 3, 4))
shade <- build_shade_rects(p)
cat("  PRODUCER_BLOCKS fixture", nrow(shade), "cat_names:",
    if (is.null(names(shade$cat))) "none" else paste(names(shade$cat), collapse=","), "\n")
stopifnot(nrow(shade) == 2L, is.null(names(shade$cat)),
          identical(shade$cat, c("cynipoid", "cynipoid")),
          isTRUE(all.equal(shade$ymin, c(0.51, 3.51))),
          isTRUE(all.equal(shade$ymax, c(2.49, 4.49))))
stopifnot(exists("check_shade_geometry"), exists("shade_runs"))
runs <- shade_runs(p$data)
stopifnot(length(runs) == 2L,
          identical(runs[[1]]$tips, c("cynipoid|BK_F3__b_1", "cynipoid|DQ_F1__a_2")),
          identical(runs[[2]]$tips, "cynipoid|DC_F2__c_1"))
g <- check_shade_geometry(shade, p$data, "fixture")
stopifnot(max(g$mismatch_rows) < 1e-9)

# A category change without a separator also splits the block.
p2 <- tips(c("cynipoid|A_F1__x_1", "x__fungus_Aspergillus", "y__fungus_Fusarium"),
           1:3)
s2 <- build_shade_rects(p2)
stopifnot(nrow(s2) == 2L, identical(s2$cat, c("cynipoid", "fungus")))
check_shade_geometry(s2, p2$data, "fixture_change")

fails <- function(expr) inherits(tryCatch(expr, error=function(e) e), "error")
# Displacement beyond the 0.5-row bound must fail; within it must pass.
bad <- shade; bad$ymax[1] <- bad$ymax[1] + 0.6
stopifnot(fails(check_shade_geometry(bad, p$data, "displaced_0.6")))
bad <- shade; bad$ymin[2] <- bad$ymin[2] - 0.51
stopifnot(fails(check_shade_geometry(bad, p$data, "displaced_0.51")))
ok <- shade; ok$ymax[1] <- ok$ymax[1] + 0.3
check_shade_geometry(ok, p$data, "displaced_0.3")
# A missing block or a wrong category must fail.
stopifnot(fails(check_shade_geometry(shade[1, ], p$data, "dropped")))
bad <- shade; bad$cat[2] <- "fungus"
stopifnot(fails(check_shade_geometry(bad, p$data, "wrong_class")))
# A named category vector (the old failure) must still compare by value.
named <- shade; named$cat <- setNames(named$cat, c("t1", "t2"))
check_shade_geometry(named, p$data, "named_cat")

cat("S4_SHADE_RUNS_PASS\n")

#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 1L)
source_lines <- readLines(args[1], warn=FALSE)
first <- grep("^plotmath_quote <- function", source_lines)
last <- grep("^# .*Parse domain FAA", source_lines)
stopifnot(length(first) == 1L, length(last) == 1L, first < last)
eval(parse(text=source_lines[first:(last - 1L)]))

cases <- c(
  "Archangium" = 'italic("Archangium")',
  "Corallococcus" = 'italic("Corallococcus")',
  "Archangium gephyra" = 'italic("Archangium gephyra")*plain("")',
  "Archangium sp." = 'italic("Archangium")*plain(" sp.")',
  "Archangium gephyra - WP_123.1 - domain 1" =
    'italic("Archangium gephyra")*plain(" - WP_123.1 - domain 1")',
  "WP_123.1 - domain 1" = 'plain("WP_123.1 - domain 1")',
  "unclassified Corallococcus" =
    'plain("unclassified ")*italic("Corallococcus")*plain("")',
  "unclassified Streptomyces - WP_1.1" =
    'plain("unclassified ")*italic("Streptomyces")*plain(" - WP_1.1")',
  "Ceroptres sp. ex Druon ignotum" =
    'italic("Ceroptres")*plain(" sp. ex ")*italic("Druon ignotum")*plain("")',
  "whitefly LGT BtaF20 - XP_018906980" =
    'plain("whitefly LGT BtaF20 - XP_018906980")'
)
actual <- vapply(names(cases), italicize_s4_label, character(1))
stopifnot(identical(unname(actual), unname(cases)))
cat("S4_LABEL_ITALICS_PASS", length(cases), "cases\n")

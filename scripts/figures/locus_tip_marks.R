locus_tip_marks <- function(family, tips) {
  path <- Sys.getenv("PCWDE_LOCUS_TIP_CROSSWALK", unset = "")
  if (!nzchar(path)) stop("PCWDE_LOCUS_TIP_CROSSWALK is required")
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  x <- x[x$family == family, ]
  if (anyDuplicated(x$display_tree_tip)) stop("duplicate display tip in crosswalk: ", family)
  i <- match(tips, x$display_tree_tip)
  if (anyNA(i)) stop("missing tip in crosswalk: ", family, ": ", paste(tips[is.na(i)], collapse = ", "))
  kind <- x$source_kind[i]
  chr <- x$corrected_chromosome[i]
  if (any(!kind %in% c("genome_copy", "transcript_assembly"))) stop("unknown source kind: ", family)
  if (any(kind == "genome_copy" & is.na(chr))) stop("missing genome placement: ", family)
  mark <- rep("", length(tips))
  mark[kind == "transcript_assembly"] <- "T"
  mark[kind == "genome_copy" & chr == "unplaced"] <- "U"
  mark
}

locus_tip_chromosomes <- function(family, tips) {
  # BK chromosome per displayed tip from the built tip table (genome copies; "" for transcripts).
  path <- Sys.getenv("PCWDE_LOCUS_TIP_CROSSWALK", unset = "")
  if (!nzchar(path)) stop("PCWDE_LOCUS_TIP_CROSSWALK is required")
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  x <- x[x$family == family, ]
  i <- match(tips, x$display_tree_tip)
  if (anyNA(i)) stop("missing tip in crosswalk: ", family, ": ", paste(tips[is.na(i)], collapse = ", "))
  chr <- x$corrected_chromosome[i]
  chr[is.na(chr)] <- ""
  chr
}

write_locus_chr_labels <- function(dir, family, view, tips, labels) {
  # Record the chromosome-column text actually drawn, for the build's consistency check.
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  write.table(
    data.frame(family = family, view = view, display_tree_tip = tips, label = unname(labels)),
    file.path(dir, paste0("chromosome_labels_", family, "_", view, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

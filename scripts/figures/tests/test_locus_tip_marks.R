script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
source(file.path(dirname(dirname(normalizePath(script))), "locus_tip_marks.R"))
tab <- data.frame(family="GH5_2", display_tree_tip=c("g1","g2","t1"),
                  source_kind=c("genome_copy","genome_copy","transcript_assembly"),
                  corrected_chromosome=c("chr6","unplaced",""))
path <- tempfile(fileext=".tsv")
write.table(tab,path,sep="\t",quote=FALSE,row.names=FALSE)
Sys.setenv(PCWDE_LOCUS_TIP_CROSSWALK=path)
stopifnot(identical(locus_tip_marks("GH5_2", c("g1","g2","t1")), c("","U","T")))
stopifnot(inherits(try(locus_tip_marks("GH5_2", "missing"),silent=TRUE),"try-error"))
unlink(path)

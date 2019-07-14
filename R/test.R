# library(data.table)
# library(dplyr)
# library(doParallel)
# library(gtools)
#
#
#
# pedigree <- convertDMATRIX(sample.info=sample.info,
#                            branch.points=branch.points,
#                            dmatrix=dmatrix,
#                            design="sibling")
# # #
# props.name <- read.table(paste("inst/extdata/dm/", "rc.csv", sep=""), sep="\t", header=T, stringsAsFactors = FALSE)
# ### File 2: input file containing information on generation times and pedigree lineages
# sample.info <- read.table(paste("inst/extdata/dm/", "sampleInfo.csv", sep=""), sep="\t", header=T, stringsAsFactors = FALSE)
# ### File 3: input file containing lineage branch points
# branch.points <- read.table(paste("inst/extdata/dm/", "branchPoints.csv", sep=""), sep="\t", header=T, stringsAsFactors = FALSE)
# ### File 4: input file containing 5mC divergence values for each sample pair
# dmatrix <- read.table(paste("inst/extdata/dm/", "d-matrix-CG.csv", sep=""), sep="\t", header=T, stringsAsFactors = FALSE)
# context <- "CG"
#
# outliers <- "none"
# dmatrix <- dmatrix[which(dmatrix[,1] != outliers), ]
# dmatrix <- dmatrix[which(dmatrix[,2] != outliers), ]
# pedigree <- pedigree[c(as.numeric(rownames(dmatrix))),]
#
# props <- props.name[which(as.character(props.name[,2]) == context),]
# props <- props.name[which(!is.element(props.name[,1], outliers) == TRUE),]
# p0uu_in <- 1-mean(as.numeric(as.character(props[,3])))
# p0uu_in
#
#
#
#
# file1 <- "inst/extdata/dm/pedigree.csv"
# pedigree <- as.matrix(read.table(file1,sep=",", header=TRUE, stringsAsFactors = FALSE))
# p0uu_in <- 0.7435074
# eqp.weight <- 1
# Nstarts <- 5
# output.data.dir <- paste0( getwd(),"/")
# out.name <- "CG_global_estimates_ABneutral"
# b<-as.matrix(read.table("/home/yadi/SERVER/Rpackages/alphabeta/inst/extdata/dm/pedigree.csv",sep=",", header=TRUE, stringsAsFactors = FALSE))
# out1 <- ABnull(pedigree.data = pedigree,
#
#                   out.dir=output.data.dir,
#                   out.name="CG_global_estimates_ABneutral")

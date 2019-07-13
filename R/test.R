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
# props.name <- read.table(paste("/home/yadi/TESTING/test-alphabeta/", "global_MA1_1_CG.csv", sep=""), sep="\t", header=T, stringsAsFactors = FALSE)
# ### File 2: input file containing information on generation times and pedigree lineages
# sample.info <- read.table(paste("/home/yadi/TESTING/test-alphabeta/", "Sample_info.csv", sep=""), sep="\t", header=T, stringsAsFactors = FALSE)
# ### File 3: input file containing lineage branch points
# branch.points <- read.table(paste("/home/yadi/TESTING/test-alphabeta/", "branchPoints.csv", sep=""), sep="\t", header=T, stringsAsFactors = FALSE)
# ### File 4: input file containing 5mC divergence values for each sample pair
# dmatrix <- read.table(paste("/home/yadi/TESTING/test-alphabeta/", "d-matrix-global-MA1_1-CG.csv", sep=""), sep="\t", header=T, stringsAsFactors = FALSE)
#


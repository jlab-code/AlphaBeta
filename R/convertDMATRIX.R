#' convertDMATRIX
#'
#' calculate divergence times of the pedigree
#'
#' @param sample.info input file containing information on generation times and pedigree lineages
#' @param branch.points input file containing lineage branch points
#' @param dmatrix input file containing 5mC divergence values for each sample pair
#' @param design "sibling" or "direct"
#' @return pedigree
#' @export
#' @examples
#'## Get some toy data
#' file1 <- system.file("extdata/dm/","sampleInfo.csv", package="alphabeta")
#' file2<-system.file("extdata/dm/","branchPoints.csv", package="alphabeta")
#' file3<-system.file("extdata/dm/","d-matrix-CG.csv", package="alphabeta")
#' sample.info <-read.table(file1,sep="\t", header=TRUE, stringsAsFactors = FALSE)
#' branch.points <-read.table(file2,sep="\t", header=TRUE, stringsAsFactors = FALSE)
#' dmatrix <-read.table(file3,sep="\t", header=TRUE, stringsAsFactors = FALSE)
#' pedigree <- convertDMATRIX(sample.info=sample.info,
#'  branch.points=branch.points,dmatrix=dmatrix,design="sibling")
#' head(pedigree)



convertDMATRIX<-function(sample.info, branch.points, dmatrix, design)
{

  #sample.info=sample.info
  #branch.points=branch.points
  #dmatrix=dmatrix
  #design="direct"

  si<-sample.info
  bp<-branch.points

  time0<-NULL
  time1<-NULL
  time2<-NULL

  for (a in seq_len(nrow(dmatrix)))
  {

      dtemp<-dmatrix[a,]
      si1<-si[which(as.character(si[,1]) == as.character(dtemp[,1])),]
      si2<-si[which(as.character(si[,1]) == as.character(dtemp[,2])),]

        if (as.character(si1[,3]) != as.character(si2[,3]))
        {
            time1[a]<-as.numeric(si1[,2])
            time2[a]<-as.numeric(si2[,2])
            time0[a]<-0
        }
        if (as.character(si1[,3]) == as.character(si2[,3]))
        {
            min.si<-min(as.numeric(si1[,2]), as.numeric(si2[,2]))
            bptemp<-bp[which(bp[,3] == si1[,3]),]

            if (design == "sibling")
            {
                 bptemp<-bptemp[which(bptemp[,2] < min.si),]
                 bptemp<-bptemp[which.max(bptemp[,2]),]
                 time1[a]<-as.numeric(si1[,2])
                 time2[a]<-as.numeric(si2[,2])
                 time0[a]<-as.numeric(bptemp[,2])
            }
            if (design == "direct")
            {
                bptemp<-bptemp[which(bptemp[,2] <= min.si),]
                bptemp<-bptemp[which.max(bptemp[,2]),]
                time1[a]<-as.numeric(si1[,2])
                time2[a]<-as.numeric(si2[,2])
                time0[a]<-as.numeric(bptemp[,2])
            }
        }
  }




   dmatrixout<-cbind(time0, time1, time2, dmatrix[,3])
   colnames(dmatrixout)[4]<-"D.value"

  return(dmatrixout)

}

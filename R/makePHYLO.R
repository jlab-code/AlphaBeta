#' Generate pedigree
#' @param tall pedigree data.
#' @param pedigree dmatrix.
#' @param sample.info sample.info
#' @return pedigree data.
#' @export
#' @examples
#'## Get some toy data
#' file1 <- system.file("extdata/soma/","D-matrix_CG.csv", package="AlphaBeta")
#' file2 <- system.file("extdata/soma/","sample_info.csv", package="AlphaBeta")
#' d.matrix <- read.table(file1, sep="\t", header=TRUE, stringsAsFactors = FALSE)
#' sample.info <- read.table(file2, sep="\t", header=TRUE, stringsAsFactors = FALSE)
#' out <- makePHYLO(tall=330, pedigree = d.matrix, sample.info = sample.info)
#'


makePHYLO<-function(tall, pedigree, sample.info)
{

  t13<-tall
  t14<-tall
  s13<-sample.info[which(sample.info[,"Stem"] == 13),]
  s13[,"Branchpoint_date"]<-t13- s13[,"Branchpoint_date"]
  s14<-sample.info[which(sample.info[,"Stem"] == 14),]
  s14[,"Branchpoint_date"]<-t14 - s14[,"Branchpoint_date"]
  sample.info<-rbind(s13, s14)

  time1<-NULL
  time2<-NULL
  time0<-NULL
  delta.t<-NULL

  for (a in 1:nrow(pedigree))
  {

    temp0<-pedigree[a,]
    pair1<-sample.info[which(as.character(sample.info[,1]) == as.character(temp0[,1])),]
    pair2<-sample.info[which(as.character(sample.info[,1]) == as.character(temp0[,2])),]

    if (pair1[,"Stem"] == pair2[,"Stem"])
    {
      t1<-pair1[,2]+pair1[,3]
      t2<-pair2[,2]+pair2[,3]
      t0<-min(pair1[,3], pair2[,3])
      time1[a]<-t1
      time2[a]<-t2
      time0[a]<-t0

    }

    if (pair1[,"Stem"] != pair2[,"Stem"])
    {
      t1<-pair1[,2]+pair1[,3]
      t2<-pair2[,2]+pair2[,3]
      t0<-min(pair1[,3], pair2[,3])
      time1[a]<-t1
      time2[a]<-t2
      time0[a]<-0
    }


  }

  delta.t<-time1 + time2 - 2*time0
  pdata<-cbind(pedigree, time1, time2, time0, delta.t)

  pslim<-cbind(pdata[,"time0"],  pdata[,"time1"], pdata[,"time2"], pdata[,"D.value"])
  colnames(pslim)<-c("time0", "time1", "time2", "D.value")

  out<-list(pslim, pdata)

} # End of function






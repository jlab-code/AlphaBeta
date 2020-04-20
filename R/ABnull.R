#' Run model that considers no accumulation of epimutations (ABnull)
#' @param pedigree.data Generation table name, you can find sample file in
#' @param out.dir outputdirectory
#' @param out.name name of file
#' @return ABnull RData file.
#' @export
#' @examples
#' #Get some toy data
#' inFile <- readRDS(system.file("extdata/dm/","output.rds", package="AlphaBeta"))
#' pedigree <- inFile$Pdata
#' out.name <- "CG_global_estimates_ABnull"
#' out <- ABnull(pedigree.data = pedigree,
#'                   out.dir=getwd(),
#'                   out.name=out.name)
#'
#' summary(out)
#'


ABnull<-function(pedigree.data, out.dir, out.name)
{

    pedigree<-pedigree.data
    delta.t<-pedigree[,2] + pedigree[,3] - 2*pedigree[,1]
    div.pred<-mean(pedigree[,"D.value"])
    residual<-pedigree[,"D.value"] - div.pred
    pedigree<-cbind(pedigree, delta.t, div.pred, residual)
    colnames(pedigree)[4:7]<-c("div.obs", "delta.t", "div.pred", "residual")

    model<-"ABnull.R"
    final.1<-div.pred
    final.2<-NULL
    info.out<-NULL
    pedigree.new<-cbind(pedigree[,1:3], div.pred, delta.t)
    colnames(pedigree.new)<-c("time0", "time1", "time2", "div.sim", "delta.t")
    abfree.out<-list(final.1, final.2, pedigree, info.out, model, pedigree.new)
    names(abfree.out)<-c("estimates", "estimates.flagged", "pedigree", "settings", "model", "for.fit.plot")

    ## Ouputting result datasets
    dput(abfree.out, paste0(out.dir,"/", out.name, ".Rdata", sep=""))
    return(abfree.out)

} #End of function


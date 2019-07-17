#' Comparison of different models and selection of best model
#' @param pedigree.select pedigree model.
#' @param pedigree.null ABnull pedigree.
#' @importFrom  stats pf
#' @return outfinal compared data.
#' @export
#' @examples
#'## Get some toy data
#' file1 <- system.file("extdata/models/","CG_global_estimates_ABneutral.Rdata", package="AlphaBeta")
#' file2 <- system.file("extdata/models/","CG_global_estimates_ABnull.Rdata", package="AlphaBeta")
#' out <- FtestRSS(pedigree.select=file1,
#'                 pedigree.null=file2)


FtestRSS<-function(pedigree.select, pedigree.null)
{

    # Reading data
    est<-dget(pedigree.select)
    estN<-dget(pedigree.null)

    if (estN$model == "ABnull.R")
    {
          # Testing
          RSSf<-est$estimates[1,"value"]
          RSSr<-sum((estN$pedigree[,"residual"])^2)
          Npara_r<-1
          Npara_f<-5
          dfF<-length(est$pedigree[,"residual"])-5
          dfR<-length(estN$pedigree[,"residual"])-1
          dfN<-dfR-dfF
          Fvalue<-((RSSr - RSSf)/(Npara_f-Npara_r))/(RSSf/dfF)
          pvalue<-pf(Fvalue, dfN, dfF, lower.tail=FALSE)
          output<-c(RSSf, RSSr, dfF, dfR, Fvalue, pvalue)
          names(output)<-c("RSS_F", "RSS_R", "df_F", "df_R", "Fvalue", "pvalue")
          outfinal<-list(output, est$estimates, estN$estimates)
          names(outfinal)<-c("Ftest", "est.selection", "est.neutral")
    }

    if (estN$model != "ABnull.R")
    {
      # Testing
      RSSf<-est$estimates[1,"value"]
      RSSr<-estN$estimates[1,"value"]
      Npara_r<-4
      Npara_f<-5
      dfF<-length(est$pedigree[,"residual"])-5
      dfR<-length(estN$pedigree[,"residual"])-4
      dfN<-dfR-dfF
      Fvalue<-((RSSr - RSSf)/(Npara_f-Npara_r))/(RSSf/dfF)
      pvalue<-pf(Fvalue, dfN, dfF, lower.tail=FALSE)
      output<-c(RSSf, RSSr, dfF, dfR, Fvalue, pvalue)
      names(output)<-c("RSS_F", "RSS_R", "df_F", "df_R", "Fvalue", "pvalue")
      outfinal<-list(output, est$estimates, estN$estimates)
      names(outfinal)<-c("Ftest", "est.selection", "est.neutral")

    }

    outfinal
}

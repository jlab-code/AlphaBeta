#' Plotting estimates
#'
#' Plotting Estimating epimutation
#' @param pedigree.names Models output AB*.Rdata
#' @param output.dir output directory
#' @param out.name filename
#' @param alpha ggplot parameters
#' @param geom.point.size ggplot parameters
#' @param geom.line.size ggplot parameters
#' @param plot.height ggplot parameters
#' @param plot.width ggplot parameters
#' @param plot.type type of plot (data.only, fit.only, both)
#' @param lsq.line Least Square Regression line (theory or pred)
#' @param intract  to see intarctive plot. (useing plotly)
#' @import  ggplot2
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom stats aggregate
#' @importFrom stats  median
#' @importFrom plotly ggplotly
#' @return plot
#' @export
#' @examples
#' # Get some toy data
#' file <- system.file("extdata/dm/","Col_CG_global_estimates_ABneutral.Rdata", package="AlphaBeta")
#' ABplot(pedigree.names=file, output.dir=getwd(), out.name="ABneutral")



ABplot <- function(pedigree.names, output.dir, out.name, alpha=0.5, geom.point.size=2,
  geom.line.size=0.9, plot.height=8, plot.width=11, plot.type="both", lsq.line="theory", intract=FALSE) {

  #initialize vectors
  dfDiv <- c()
  predDfFits <- c()
  theoryDfFits <- c()
  gg <- ggplot()

  #NCOL and NROW do the same, treating a vector as 1-column matrix.

    #for(j in 1:NCOL(pedigree.names)) {

    datain <- dget(pedigree.names)
    name <- gsub("(_).*", "", basename(as.character(pedigree.names)))
    category <- gsub(pattern=paste0(name, "_|\\.Rdata$"), "", basename(as.character(pedigree.names)))

    div <- as.data.frame(datain$pedigree[,c("delta.t","div.obs")])
    name.div <- cbind(category, name, div)
    dfDiv <- rbind(dfDiv, name.div)
    intercept <- datain$estimates[1, "intercept"]

    #predictive fit
    pred.fit.data <- aggregate(datain$pedigree[,"div.pred"], by=list(datain$pedigree[,"delta.t"]), median)
    pred.fits <- c(intercept, pred.fit.data[,2])
    pred.fit.t <- c(0, pred.fit.data[,1])

    predDfFits <- rbind(predDfFits, data.frame(category, name, pred.fit.t, pred.fits))

    #theoritical fit
    theory.fit.data <- datain$for.fit.plot
    theory.fits <- c(intercept, theory.fit.data[,"div.sim"])
    theory.fit.t <- c(0, theory.fit.data[,"delta.t"])

    theoryDfFits <- rbind(theoryDfFits, data.frame(category, name, theory.fit.t, theory.fits))

  # }

  if ((plot.type == "data.only")) {
    gg <- gg + geom_point(data=dfDiv, aes(x=delta.t, y=div.obs, colour=name), alpha = alpha, size=geom.point.size)
  }

  if ((plot.type == "fit.only") && (lsq.line == "pred")) {
    gg <- gg + geom_line(data=predDfFits, aes(x=pred.fit.t, y=pred.fits, color=name), size=geom.line.size)
  }

  if ((plot.type == "fit.only") && (lsq.line == "theory")) {
    gg <- gg + geom_line(data=theoryDfFits, aes(x=theory.fit.t, y=theory.fits, color=name), size=geom.line.size)
  }

  if ((plot.type == "both") && (lsq.line == "pred")) {
    gg <- gg + geom_point(data=dfDiv, aes(x=delta.t, y=div.obs, colour=name), alpha = alpha, size=geom.point.size) +
      geom_line(data=predDfFits, aes(x=pred.fit.t, y=pred.fits, color=name), size=geom.line.size)
  }

  if ((plot.type == "both") && (lsq.line == "theory")) {
    gg <- gg + geom_point(data=dfDiv, aes(x=delta.t, y=div.obs, colour=name), alpha = alpha, size=geom.point.size) +
      geom_line(data=theoryDfFits, aes(x=theory.fit.t, y=theory.fits, color=name), size=geom.line.size)
  }

  #adjust x, y-axis
  ymax <- max(dfDiv$div.obs) + 0.005
  xmax <- max(dfDiv$delta.t) + 10


  gg <- gg + labs(x="delta t (generations)", y="methylation divergence", color="category") +
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_x_continuous(expand = c(0, 0), limits=c(0, xmax)) +
    scale_y_continuous(expand = c(0, 0), limits=c(0, ymax)) +
    theme(legend.position="none")
  gg <- gg  + theme_test()

  if(!intract){
  pdf(paste0(output.dir,"/", out.name, ".pdf", sep=""), colormodel = 'cmyk', width = plot.width, height = plot.height)
  print(gg)
  dev.off()
  }else{
    pl <- ggplotly(gg)
    pl
  }

}


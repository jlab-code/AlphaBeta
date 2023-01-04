#' Bootstrap analysis with the best model
#' @param pedigree.data pedigree data.
#' @param Nboot number of boot.
#' @param out.dir output directory.
#' @param out.name output file name.
#' @import optimx
#' @import expm
#' @importFrom  stats quantile
#' @importFrom stats runif
#' @importFrom stats sd
#' @return bootstrap result.
#' @export
#' @examples
#'## Get some toy data
#' inFile <- system.file("extdata/models/","ABneutral_CG_global_estimates.Rdata", package="AlphaBeta")
#' Nboot <- 4
#' out.name <-"Boot_CG_global_estimates_ABneutral"
#' Bout <- BOOTmodel(pedigree.data=inFile,
#'                 Nboot=Nboot,
#'                 out.dir=getwd(),
#'                 out.name=out.name)
#' summary(Bout)




BOOTmodel<-function(pedigree.data, Nboot, out.dir, out.name)
{
  pedigree.data<-dget(pedigree.data)
  model<-pedigree.data$model

#################################
########### ABselectMM ##########
#################################
  if (model == "ABselectMM.R")
  {

      ## Reading the dataset for bootstrapping and extracting the parameter settings
      settings<-pedigree.data$settings
      est<-pedigree.data$estimates
      eqp<-as.numeric(as.character(settings[which(settings[,1] == "eqp"),2]))
      eqp.weight<-as.numeric(as.character(settings[which(settings[,1] == "eqp.weight"),2]))
      optim.method<-as.character(settings[which(settings[,1] == "optim.method"),2])
      #p0mm<-round(as.numeric(as.character(settings[which(settings[,1] == "p0mm"),2])),16)
      p0uu<-round(as.numeric(as.character(settings[which(settings[,1] == "p0uu"),2])),16)
      #p0um<-round(as.numeric(as.character(settings[which(settings[,1] == "p0um"),2])),15)
      pedigree<-pedigree.data$pedigree
      p0mm = 1-p0uu
      p0um = 0


    	if(sum(c(p0mm, p0um, p0uu), na.rm =TRUE) != 1) {stop("The initial state probabilities don't sum to 1")}

    	## Defining the divergence function
    	divergence <- function(pedigree, p0mm, p0um, p0uu, param)
    	{

    	  ## Initializing parameters
    	  PrMM <- p0mm
    	  PrUM <- p0um
    	  PrUU <- p0uu
    	  alpha <- param[1]
    	  bet   <- param[2]
    	  weight <-param[3]
    	  sel    <-param[4]

    	  ## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
    	  svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

    	  element11<-((1-alpha)^2)*sel
    	  element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
    	  element13<-(alpha^2)
    	  rowtotal1<-element11 + element12 + element13

    	  element21<-(1/4*(bet + 1 - alpha)^2)*sel
    	  element22<-(1/2*(bet + 1 - alpha)*(alpha + (1 - bet)))*(1/2*(1+sel))
    	  element23<-(1/4*(alpha + (1 - bet))^2)
    	  rowtotal2<-element21 + element22 + element23

    	  element31<-(bet^2)*sel
    	  element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
    	  element33<-((1-bet))^2
    	  rowtotal3<-element31 + element32 + element33

    	  ## Defining the generation (or transition) matrix
    	  Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
    	                        element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
    	                        element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

    	  ## Calculating the expected divergence for every observed pair in 'pedigree.txt'
    	  Dt1t2<-NULL

    	  for (p in seq_len(NROW(pedigree)))
    	  {

    	    ## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
    	    svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
    	    svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
    	    svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
    	    svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
    	    svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
    	    svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
    	    svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

    	    ## Conditional divergences
    	    dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
    	                        svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

    	    dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
    	                        svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

    	    dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
    	                        svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

    	    ## Total (weighted) divergence
    	    Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


    	  }

    	  ## Pr(UU) at equilibrium given alpha and beta; Note: this only approximates the equilibrium values
    	  puuinf.est<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
    	  puuinf.est<- puuinf.est[1,1]
    	  divout<-list(puuinf.est, Dt1t2)

    	  return(divout)

    	}

    	## Defining the Least Square function to be minimized
    	LSE_intercept<-function(param_int)
    	{
    	  d = divergence(pedigree, p0mm, p0um, p0uu, param_int[1:4])
    	  sum((pedigree[,4] - param_int[5] - divergence(pedigree, p0mm, p0um, p0uu, param_int[1:4])[[2]])^2) +
    	    eqp.weight*nrow(pedigree)*((divergence(pedigree, p0mm, p0um, p0uu, param_int[1:4])[[1]]- eqp)^2)
    	}


      ## Initializing
    	final<-NULL
    	counter<-0


    	## Defining starting values
    	alpha.start  <-est[1,1]
    	beta.start   <-est[1,2]
    	weight.start <-est[1,3]
    	sel.start<-est[1,4]
    	intercept.start <-est[1,5]
    	param_int0 = c(alpha.start, beta.start, weight.start, sel.start, intercept.start)


    	    ## Start of boostrap loops
        	for (booting in seq_len(Nboot))
        	{

        	  opt.out<-NULL
        		pedigree[,"div.obs"]<-pedigree[,"div.pred"]+sample(pedigree[,"residual"], nrow(pedigree), replace=TRUE)

        			counter<-counter+1

        			message("Bootstrap interationsssss: ", counter/Nboot, "\n")

        			opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
        			alphafinal<-as.numeric(opt.out[1])
        			betfinal<-as.numeric(opt.out[2])
        			weightfinal<-as.numeric(opt.out[3])
        			selfinal<-as.numeric(opt.out[4])

        			## Calculating equilibrium frequencies based on the model estimates
        			svGzero   <- c(p0uu, (weightfinal)*p0mm, (1-weightfinal)*p0mm)

        			element11<-((1-alphafinal)^2)*selfinal
        			element12<-(2*(1-alphafinal)*alphafinal)*(1/2*(1+selfinal))
        			element13<-(alphafinal^2)
        			rowtotal1<-element11 + element12 + element13

        			element21<-(1/4*(betfinal + 1 - alphafinal)^2)*selfinal
        			element22<-(1/2*(betfinal + 1 - alphafinal)*(alphafinal + (1 - betfinal)))*(1/2*(1+selfinal))
        			element23<-(1/4*(alphafinal + (1 - betfinal))^2)
        			rowtotal2<-element21 + element22 + element23

        			element31<-(betfinal^2)*selfinal
        			element32<-(2*((1-betfinal))*betfinal)*(1/2*(1+selfinal))
        			element33<-((1-betfinal))^2
        			rowtotal3<-element31 + element32 + element33

        			## Defining the generation (or transition) matrix
        			Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
        			                      element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
        			                      element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

        			## Note: This is an approximation to the equilibrium values
        			pinf.vec<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
        			PrMMinf<- pinf.vec[1,3]
        			PrUMinf<- pinf.vec[1,2]
        			PrUUinf<- pinf.vec[1,1]
        			opt.out <-cbind(opt.out, PrMMinf, PrUMinf, PrUUinf, alpha.start, beta.start, weight.start,
        			                sel.start, intercept.start)

        			## Collecting the results
        			final<-rbind(final, opt.out)


        	} # End of Bootstrap loops


    	 colnames(final)[1:5]<-c("alpha", "beta", "weight", "sel.coef", "intercept")
    	 colnames(final)[14:16]<-c("PrMMinf", "PrUMinf", "PrUUinf")

    	 SE.alpha<-sd(final[,1],na.rm=TRUE)
    	 SE.beta<-sd(final[,2],na.rm=TRUE)
    	 SE.beta.alpha<-sd(final[,2]/final[,1], na.rm=TRUE)
    	 SE.weight<-sd(final[,3],na.rm=TRUE)
    	 SE.sel.coef<-sd(final[,4], na.rm=TRUE)
    	 SE.inter<-sd(final[,5],na.rm=TRUE)
    	 SE.PrMMinf<-sd(final[,14],na.rm=TRUE)
    	 SE.PrUMinf<-sd(final[,15],na.rm=TRUE)
    	 SE.PrUUinf<-sd(final[,16],na.rm=TRUE)

    	 CI.alpha<-quantile(final[,1],probs=c(0.025, 0.975))
    	 CI.beta<-quantile(final[,2],probs=c(0.025, 0.975))
    	 CI.beta.alpha<-quantile(final[,2]/final[,1], probs=c(0.025, 0.97))
    	 CI.weight<-quantile(final[,3],probs=c(0.025, 0.975))
    	 CI.sel.coef<-quantile(final[,4],probs=c(0.025, 0.975))
    	 CI.inter<-quantile(final[,5],probs=c(0.025, 0.975))
    	 CI.PrMMinf<-quantile(final[,14],probs=c(0.025, 0.975))
    	 CI.PrUMinf<-quantile(final[,15],probs=c(0.025, 0.975))
    	 CI.PrUUinf<-quantile(final[,16],probs=c(0.025, 0.975))

    	 SE<-c(SE.alpha, SE.beta, SE.beta.alpha, SE.weight, SE.sel.coef, SE.inter, SE.PrMMinf, SE.PrUMinf, SE.PrUUinf)
    	 CI<-rbind(CI.alpha, CI.beta, CI.beta.alpha, CI.weight, CI.sel.coef, CI.inter, CI.PrMMinf, CI.PrUMinf, CI.PrUUinf)

    	 SE.out<-cbind(SE, CI)
    	 colnames(SE.out)[1]<-"SE"
    	 rownames(SE.out)<-c("alpha", "beta", "beta/alpha", "weight", "sel.coef", "intercept", "PrMMinf", "PrUMinf", "PrUUinf")

    	 final<-data.frame(final)
    	 good.boots<-length(which(is.na(final[,"alpha"]) == FALSE))

    	 SE.out<-list(SE.out, est[1,], settings, Nboot, good.boots, final, model)
    	 names(SE.out)<-c("standard.errors", "boot.base", "settings", "N.boots", "N.good.boots", "boot.results", "model")

    	 ## Ouputting result datasets
    	 dput(SE.out, paste0(out.dir,"/", out.name, ".Rdata", sep=""))

    	 return(SE.out)

  } # End of "ABselectMM" if statement


  #################################
  ########### ABselectUU ##########
  #################################
  if (model == "ABselectUU.R")
  {

    ## Reading the dataset for bootstrapping and extracting the parameter settings
    settings<-pedigree.data$settings
    est<-pedigree.data$estimates
    eqp<-as.numeric(as.character(settings[which(settings[,1] == "eqp"),2]))
    eqp.weight<-as.numeric(as.character(settings[which(settings[,1] == "eqp.weight"),2]))
    optim.method<-as.character(settings[which(settings[,1] == "optim.method"),2])
    p0uu<-round(as.numeric(as.character(settings[which(settings[,1] == "p0uu"),2])),16)
    pedigree<-pedigree.data$pedigree
    p0mm = 1-p0uu
    p0um = 0


    if(sum(c(p0mm, p0um, p0uu), na.rm =TRUE) != 1) {stop("The initial state probabilities don't sum to 1")}

    ## Defining the divergence function
    divergence <- function(pedigree, p0mm, p0um, p0uu, param)
    {

      ## Initializing parameters
      PrMM <- p0mm
      PrUM <- p0um
      PrUU <- p0uu
      alpha <- param[1]
      bet   <- param[2]
      weight <-param[3]
      sel    <-param[4]

      ## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
      svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

      element11<-((1-alpha)^2)
      element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
      element13<-(alpha^2)*sel
      rowtotal1<-element11 + element12 + element13

      element21<-(1/4*(bet + 1 - alpha)^2)
      element22<-(1/2*(bet + 1 - alpha)*(alpha + (1 - bet)))*(1/2*(1+sel))
      element23<-(1/4*(alpha + (1 - bet))^2)*sel
      rowtotal2<-element21 + element22 + element23

      element31<-(bet^2)
      element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
      element33<-(((1-bet))^2)*sel
      rowtotal3<-element31 + element32 + element33

      ## Defining the generation (or transition) matrix
      Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
                            element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
                            element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

      ## Calculating the expected divergence for every observed pair in 'pedigree.txt'
      Dt1t2<-NULL

      for (p in seq_len(NROW(pedigree)))
      {

        ## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
        svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
        svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
        svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
        svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
        svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
        svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
        svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

        ## Conditional divergences
        dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
                            svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

        dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
                            svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

        dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
                            svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

        ## Total (weighted) divergence
        Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


      }

      ## Pr(UU) at equilibrium given alpha and beta; Note: this only approximates the equilibrium values
      puuinf.est<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
      puuinf.est<- puuinf.est[1,1]
      divout<-list(puuinf.est, Dt1t2)

      return(divout)

    }

    ## Defining the Least Square function to be minimized
    LSE_intercept<-function(param_int)
    {
      sum((pedigree[,4] - param_int[5] - divergence(pedigree, p0mm, p0um, p0uu, param_int[1:4])[[2]])^2) +
        eqp.weight*nrow(pedigree)*((divergence(pedigree, p0mm, p0um, p0uu, param_int[1:4])[[1]]- eqp)^2)
    }



    ## Initializing
    final<-NULL
    counter<-0


    ## Defining starting values
    alpha.start  <-est[1,1]
    beta.start   <-est[1,2]
    weight.start <-est[1,3]
    sel.start<-est[1,4]
    intercept.start <-est[1,5]
    param_int0 = c(alpha.start, beta.start, weight.start, sel.start, intercept.start)


    ## Start of boostrap loops
    for (booting in seq_len(Nboot))
    {

      opt.out<-NULL
      pedigree[,"div.obs"]<-pedigree[,"div.pred"]+sample(pedigree[,"residual"], nrow(pedigree), replace=TRUE)

      counter<-counter+1

      message("Bootstrap interation: ", counter/Nboot, "\n")

      opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
      alphafinal<-as.numeric(opt.out[1])
      betfinal<-as.numeric(opt.out[2])
      weightfinal<-as.numeric(opt.out[3])
      selfinal<-as.numeric(opt.out[4])

      ## Calculating equilibrium frequencies based on the model estimates
      svGzero   <- c(p0uu, (weightfinal)*p0mm, (1-weightfinal)*p0mm)

      element11<-((1-alphafinal)^2)
      element12<-(2*(1-alphafinal)*alphafinal)*(1/2*(1+selfinal))
      element13<-(alphafinal^2)*selfinal
      rowtotal1<-element11 + element12 + element13

      element21<-(1/4*(betfinal + 1 - alphafinal)^2)
      element22<-(1/2*(betfinal + 1 - alphafinal)*(alphafinal + (1 - betfinal)))*(1/2*(1+selfinal))
      element23<-(1/4*(alphafinal + (1 - betfinal))^2)*selfinal
      rowtotal2<-element21 + element22 + element23

      element31<-(betfinal^2)
      element32<-(2*((1-betfinal))*betfinal)*(1/2*(1+selfinal))
      element33<-(((1-betfinal))^2)*selfinal
      rowtotal3<-element31 + element32 + element33

      ## Defining the generation (or transition) matrix
      Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
                            element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
                            element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

      ## Note: This is an approximation to the equilibrium values
      pinf.vec<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
      PrMMinf<- pinf.vec[1,3]
      PrUMinf<- pinf.vec[1,2]
      PrUUinf<- pinf.vec[1,1]
      opt.out <-cbind(opt.out, PrMMinf, PrUMinf, PrUUinf, alpha.start, beta.start, weight.start,
                      sel.start, intercept.start)
      final<-rbind(final, opt.out)


    } # End of Bootstrap loops


    colnames(final)[1:5]<-c("alpha", "beta", "weight", "sel.coef", "intercept")
    colnames(final)[14:16]<-c("PrMMinf", "PrUMinf", "PrUUinf")

    SE.alpha<-sd(final[,1],na.rm=TRUE)
    SE.beta<-sd(final[,2],na.rm=TRUE)
    SE.beta.alpha<-sd(final[,2]/final[,1], na.rm=TRUE)
    SE.weight<-sd(final[,3],na.rm=TRUE)
    SE.sel.coef<-sd(final[,4], na.rm=TRUE)
    SE.inter<-sd(final[,5],na.rm=TRUE)
    SE.PrMMinf<-sd(final[,14],na.rm=TRUE)
    SE.PrUMinf<-sd(final[,15],na.rm=TRUE)
    SE.PrUUinf<-sd(final[,16],na.rm=TRUE)

    CI.alpha<-quantile(final[,1],probs=c(0.025, 0.975))
    CI.beta<-quantile(final[,2],probs=c(0.025, 0.975))
    CI.beta.alpha<-quantile(final[,2]/final[,1], probs=c(0.025, 0.97))
    CI.weight<-quantile(final[,3],probs=c(0.025, 0.975))
    CI.sel.coef<-quantile(final[,4],probs=c(0.025, 0.975))
    CI.inter<-quantile(final[,5],probs=c(0.025, 0.975))
    CI.PrMMinf<-quantile(final[,14],probs=c(0.025, 0.975))
    CI.PrUMinf<-quantile(final[,15],probs=c(0.025, 0.975))
    CI.PrUUinf<-quantile(final[,16],probs=c(0.025, 0.975))

    SE<-c(SE.alpha, SE.beta, SE.beta.alpha, SE.weight, SE.sel.coef, SE.inter, SE.PrMMinf, SE.PrUMinf, SE.PrUUinf)
    CI<-rbind(CI.alpha, CI.beta, CI.beta.alpha, CI.weight, CI.sel.coef, CI.inter, CI.PrMMinf, CI.PrUMinf, CI.PrUUinf)

    SE.out<-cbind(SE, CI)
    colnames(SE.out)[1]<-"SE"
    rownames(SE.out)<-c("alpha", "beta", "beta/alpha", "weight", "sel.coef", "intercept", "PrMMinf", "PrUMinf", "PrUUinf")

    final<-data.frame(final)
    good.boots<-length(which(is.na(final[,"alpha"]) == FALSE))

    SE.out<-list(SE.out, est[1,], settings, Nboot, good.boots, final, model)
    names(SE.out)<-c("standard.errors", "boot.base", "settings", "N.boots", "N.good.boots", "boot.results", "model")

    ## Ouputting result datasets
    dput(SE.out, paste0(out.dir,"/", out.name, ".Rdata", sep=""))

    return(SE.out)

  } # End of "ABselectUU" if statement





  ###################################
  ########### ABfixselectMM #########
  ###################################
  if (model == "ABfixselectMM.R")
  {

    ## Reading the dataset for bootstrapping and extracting the parameter settings
    settings<-pedigree.data$settings
    est<-pedigree.data$estimates
    eqp<-as.numeric(as.character(settings[which(settings[,1] == "eqp"),2]))
    eqp.weight<-as.numeric(as.character(settings[which(settings[,1] == "eqp.weight"),2]))
    optim.method<-as.character(settings[which(settings[,1] == "optim.method"),2])
    p0uu<-round(as.numeric(as.character(settings[which(settings[,1] == "p0uu"),2])),16)
    pedigree<-pedigree.data$pedigree
    p0mm = 1-p0uu
    p0um = 0


    if(sum(c(p0mm, p0um, p0uu), na.rm =TRUE) != 1) {stop("The initial state probabilities don't sum to 1")}

    ##### Defining the divergence function
    divergence <- function(pedigree, sel, p0mm, p0um, p0uu, param)
    {

      ## Initializing parameters
      PrMM <- p0mm
      PrUM <- p0um
      PrUU <- p0uu
      sel  <- sel
      alpha <- param[1]
      bet   <- param[2]
      weight <-param[3]


      ## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
      svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

      element11<-((1-alpha)^2)*sel
      element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
      element13<-(alpha^2)
      rowtotal1<-element11 + element12 + element13

      element21<-(1/4*(bet + 1 - alpha)^2)*sel
      element22<-(1/2*(bet + 1 - alpha)*(alpha + (1 - bet)))*(1/2*(1+sel))
      element23<-(1/4*(alpha + (1 - bet))^2)
      rowtotal2<-element21 + element22 + element23

      element31<-(bet^2)*sel
      element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
      element33<-((1-bet))^2
      rowtotal3<-element31 + element32 + element33

      ## Defining the generation (or transition) matrix
      Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
                            element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
                            element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)


      ## Calculating the expected divergence for every observed pair in 'pedigree.txt'
      Dt1t2<-NULL

      for (p in seq_len(NROW(pedigree)))
      {

        ## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
        svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
        svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
        svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
        svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
        svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
        svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
        svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

        ## Conditional divergences
        dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
                            svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

        dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
                            svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

        dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
                            svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

        ## Total (weighted) divergence
        Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


      }

      ## Pr(UU) at equilibrium given alpha and beta; Note: this only approximates the equilibrium values
      puuinf.est<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
      puuinf.est<- puuinf.est[1,1]
      divout<-list(puuinf.est, Dt1t2)

      return(divout)

    }

    ###### Defining the Least Square function to be minimized
    ###### Note the equilibrium constraint, which can be made as small as desired.

    LSE_intercept<-function(param_int)
    {
      sum((pedigree[,4] - param_int[4] - divergence(pedigree, sel, p0mm, p0um, p0uu, param_int[1:3])[[2]])^2) +
        eqp.weight*nrow(pedigree)*((divergence(pedigree, sel, p0mm, p0um, p0uu, param_int[1:3])[[1]]- eqp)^2)
    }



    ## Initializing
    final<-NULL
    counter<-0


    ## Defining starting values
    alpha.start  <-est[1,1]
    beta.start   <-est[1,2]
    weight.start <-est[1,3]
    intercept.start <-est[1,4]
    param_int0 = c(alpha.start, beta.start, weight.start, intercept.start)


    ## Start of boostrap loops
    for (booting in NROW(Nboot))
    {

      opt.out<-NULL
      pedigree[,"div.obs"]<-pedigree[,"div.pred"]+sample(pedigree[,"residual"], nrow(pedigree), replace=TRUE)

      counter<-counter+1

      message("Bootstrap interation: ", counter/Nboot, "\n")

      opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
      alphafinal<-as.numeric(opt.out[1])
      betfinal<-as.numeric(opt.out[2])
      weightfinal<-as.numeric(opt.out[3])
      selfinal<-sel


      ## Here we want to calculate the equilibrium freqs based on the estimates
      svGzero   <- c(p0uu, (weightfinal)*p0mm, (1-weightfinal)*p0mm)

      element11<-((1-alphafinal)^2)*selfinal
      element12<-(2*(1-alphafinal)*alphafinal)*(1/2*(1+selfinal))
      element13<-(alphafinal^2)
      rowtotal1<-element11 + element12 + element13

      element21<-(1/4*(betfinal + 1 - alphafinal)^2)*selfinal
      element22<-(1/2*(betfinal + 1 - alphafinal)*(alphafinal + (1 - betfinal)))*(1/2*(1+selfinal))
      element23<-(1/4*(alphafinal + (1 - betfinal))^2)
      rowtotal2<-element21 + element22 + element23

      element31<-(betfinal^2)*selfinal
      element32<-(2*((1-betfinal))*betfinal)*(1/2*(1+selfinal))
      element33<-((1-betfinal))^2
      rowtotal3<-element31 + element32 + element33

      ## Defining the generation (or transition) matrix
      Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
                            element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
                            element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)


      ## Note: This is an approximation to the equilibrium values
      pinf.vec<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
      PrMMinf<- pinf.vec[1,3]
      PrUMinf<- pinf.vec[1,2]
      PrUUinf<- pinf.vec[1,1]
      opt.out <-cbind(opt.out, PrMMinf, PrUMinf, PrUUinf, alpha.start, beta.start, weight.start, intercept.start)
      final<-rbind(final, opt.out)


    } # End of Bootstrap loops


    colnames(final)[1:4]<-c("alpha", "beta", "weight", "intercept")
    colnames(final)[13:15]<-c("PrMMinf", "PrUMinf", "PrUUinf")

    SE.alpha<-sd(final[,1],na.rm=TRUE)
    SE.beta<-sd(final[,2],na.rm=TRUE)
    SE.beta.alpha<-sd(final[,2]/final[,1], na.rm=TRUE)
    SE.weight<-sd(final[,3],na.rm=TRUE)
    SE.inter<-sd(final[,4],na.rm=TRUE)
    SE.PrMMinf<-sd(final[,13],na.rm=TRUE)
    SE.PrUMinf<-sd(final[,14],na.rm=TRUE)
    SE.PrUUinf<-sd(final[,15],na.rm=TRUE)

    CI.alpha<-quantile(final[,1],probs=c(0.025, 0.975))
    CI.beta<-quantile(final[,2],probs=c(0.025, 0.975))
    CI.beta.alpha<-quantile(final[,2]/final[,1], probs=c(0.025, 0.97))
    CI.weight<-quantile(final[,3],probs=c(0.025, 0.975))
    CI.inter<-quantile(final[,4],probs=c(0.025, 0.975))
    CI.PrMMinf<-quantile(final[,13],probs=c(0.025, 0.975))
    CI.PrUMinf<-quantile(final[,14],probs=c(0.025, 0.975))
    CI.PrUUinf<-quantile(final[,15],probs=c(0.025, 0.975))

    SE<-c(SE.alpha, SE.beta, SE.beta.alpha, SE.weight, SE.inter, SE.PrMMinf, SE.PrUMinf, SE.PrUUinf)
    CI<-rbind(CI.alpha, CI.beta, CI.beta.alpha, CI.weight, CI.inter, CI.PrMMinf, CI.PrUMinf, CI.PrUUinf)

    SE.out<-cbind(SE, CI)
    colnames(SE.out)[1]<-"SE"
    rownames(SE.out)<-c("alpha", "beta", "beta/alpha", "weight", "intercept", "PrMMinf", "PrUMinf", "PrUUinf")

    final<-data.frame(final)
    good.boots<-length(which(is.na(final[,"alpha"]) == FALSE))

    SE.out<-list(SE.out, est[1,], settings, Nboot, good.boots, final, model)
    names(SE.out)<-c("standard.errors", "boot.base", "settings", "N.boots", "N.good.boots", "boot.results", "model")

    ## Ouputting result datasets
    dput(SE.out, paste0(out.dir,"/", out.name, ".Rdata", sep=""))

    return(SE.out)

  } # End of "ABfixselectMM" if statement






  ###################################
  ########### ABfixselectUU #########
  ###################################
  if (model == "ABfixselectUU.R")
  {

    ## Reading the dataset for bootstrapping and extracting the parameter settings
    settings<-pedigree.data$settings
    est<-pedigree.data$estimates
    eqp<-as.numeric(as.character(settings[which(settings[,1] == "eqp"),2]))
    eqp.weight<-as.numeric(as.character(settings[which(settings[,1] == "eqp.weight"),2]))
    optim.method<-as.character(settings[which(settings[,1] == "optim.method"),2])
    p0uu<-round(as.numeric(as.character(settings[which(settings[,1] == "p0uu"),2])),16)
    pedigree<-pedigree.data$pedigree
    p0mm = 1-p0uu
    p0um = 0


    if(sum(c(p0mm, p0um, p0uu), na.rm =TRUE) != 1) {stop("The initial state probabilities don't sum to 1")}

    ##### Defining the divergence function
    divergence <- function(pedigree, sel, p0mm, p0um, p0uu, param)
    {

      ## Initializing parameters
      PrMM <- p0mm
      PrUM <- p0um
      PrUU <- p0uu
      sel  <- sel
      alpha <- param[1]
      bet   <- param[2]
      weight <-param[3]


    ## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
    svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

    element11<-((1-alpha)^2)
    element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
    element13<-(alpha^2)*sel
    rowtotal1<-element11 + element12 + element13

    element21<-(1/4*(bet + 1 - alpha)^2)
    element22<-(1/2*(bet + 1 - alpha)*(alpha + (1 - bet)))*(1/2*(1+sel))
    element23<-(1/4*(alpha + (1 - bet))^2)*sel
    rowtotal2<-element21 + element22 + element23

    element31<-(bet^2)
    element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
    element33<-(((1-bet))^2)*sel
    rowtotal3<-element31 + element32 + element33

    ## Defining the generation (or transition) matrix
    Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
                          element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
                          element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

    ## Calculating the expected divergence for every observed pair in 'pedigree.txt'
    Dt1t2<-NULL

    for (p in seq_len(NROW(pedigree)))
    {

      ## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
      svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
      svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
      svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
      svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

      ## Conditional divergences
      dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
                          svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

      dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
                          svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

      dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
                          svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

      ## Total (weighted) divergence
      Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


    }

    ## Pr(UU) at equilibrium given alpha and beta; Note: this only approximates the equilibrium values
    puuinf.est<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
    puuinf.est<- puuinf.est[1,1]
    divout<-list(puuinf.est, Dt1t2)

    return(divout)

  }

    ###### Defining the Least Square function to be minimized
    ###### Note the equilibrium constraint, which can be made as small as desired.

    LSE_intercept<-function(param_int)
    {
      sum((pedigree[,4] - param_int[4] - divergence(pedigree, sel, p0mm, p0um, p0uu, param_int[1:3])[[2]])^2) +
        eqp.weight*nrow(pedigree)*((divergence(pedigree, sel, p0mm, p0um, p0uu, param_int[1:3])[[1]]- eqp)^2)
    }



    ## Initializing
    final<-NULL
    counter<-0


    ## Defining starting values
    alpha.start  <-est[1,1]
    beta.start   <-est[1,2]
    weight.start <-est[1,3]
    intercept.start <-est[1,4]
    param_int0 = c(alpha.start, beta.start, weight.start, intercept.start)


    ## Start of boostrap loops
    for (booting in seq_len(Nboot))
    {

      opt.out<-NULL
      pedigree[,"div.obs"]<-pedigree[,"div.pred"]+sample(pedigree[,"residual"], nrow(pedigree), replace=TRUE)

      counter<-counter+1

      message("Bootstrap interation: ", counter/Nboot, "\n")

      opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
      alphafinal<-as.numeric(opt.out[1])
      betfinal<-as.numeric(opt.out[2])
      weightfinal<-as.numeric(opt.out[3])
      selfinal<-sel


      ## Here we want to calculate the equilibrium freqs based on the estimates
      svGzero   <- c(p0uu, (weightfinal)*p0mm, (1-weightfinal)*p0mm)

      element11<-((1-alphafinal)^2)
      element12<-(2*(1-alphafinal)*alphafinal)*(1/2*(1+selfinal))
      element13<-(alphafinal^2)*selfinal
      rowtotal1<-element11 + element12 + element13

      element21<-(1/4*(betfinal + 1 - alphafinal)^2)
      element22<-(1/2*(betfinal + 1 - alphafinal)*(alphafinal + (1 - betfinal)))*(1/2*(1+selfinal))
      element23<-(1/4*(alphafinal + (1 - betfinal))^2)*selfinal
      rowtotal2<-element21 + element22 + element23

      element31<-(betfinal^2)
      element32<-(2*((1-betfinal))*betfinal)*(1/2*(1+selfinal))
      element33<-(((1-betfinal))^2)*selfinal
      rowtotal3<-element31 + element32 + element33

      ## Defining the generation (or transition) matrix
      Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
                            element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
                            element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)


      ## Note: This is an approximation to the equilibrium values
      pinf.vec<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
      PrMMinf<- pinf.vec[1,3]
      PrUMinf<- pinf.vec[1,2]
      PrUUinf<- pinf.vec[1,1]
      opt.out <-cbind(opt.out, PrMMinf, PrUMinf, PrUUinf, alpha.start, beta.start, weight.start, intercept.start)
      final<-rbind(final, opt.out)



    } # End of Bootstrap loops


    colnames(final)[1:4]<-c("alpha", "beta", "weight", "intercept")
    colnames(final)[13:15]<-c("PrMMinf", "PrUMinf", "PrUUinf")

    SE.alpha<-sd(final[,1],na.rm=TRUE)
    SE.beta<-sd(final[,2],na.rm=TRUE)
    SE.beta.alpha<-sd(final[,2]/final[,1], na.rm=TRUE)
    SE.weight<-sd(final[,3],na.rm=TRUE)
    SE.inter<-sd(final[,4],na.rm=TRUE)
    SE.PrMMinf<-sd(final[,13],na.rm=TRUE)
    SE.PrUMinf<-sd(final[,14],na.rm=TRUE)
    SE.PrUUinf<-sd(final[,15],na.rm=TRUE)

    CI.alpha<-quantile(final[,1],probs=c(0.025, 0.975))
    CI.beta<-quantile(final[,2],probs=c(0.025, 0.975))
    CI.beta.alpha<-quantile(final[,2]/final[,1], probs=c(0.025, 0.97))
    CI.weight<-quantile(final[,3],probs=c(0.025, 0.975))
    CI.inter<-quantile(final[,4],probs=c(0.025, 0.975))
    CI.PrMMinf<-quantile(final[,13],probs=c(0.025, 0.975))
    CI.PrUMinf<-quantile(final[,14],probs=c(0.025, 0.975))
    CI.PrUUinf<-quantile(final[,15],probs=c(0.025, 0.975))

    SE<-c(SE.alpha, SE.beta, SE.beta.alpha, SE.weight, SE.inter, SE.PrMMinf, SE.PrUMinf, SE.PrUUinf)
    CI<-rbind(CI.alpha, CI.beta, CI.beta.alpha, CI.weight, CI.inter, CI.PrMMinf, CI.PrUMinf, CI.PrUUinf)

    SE.out<-cbind(SE, CI)
    colnames(SE.out)[1]<-"SE"
    rownames(SE.out)<-c("alpha", "beta", "beta/alpha", "weight", "intercept", "PrMMinf", "PrUMinf", "PrUUinf")

    final<-data.frame(final)
    good.boots<-length(which(is.na(final[,"alpha"]) == FALSE))

    SE.out<-list(SE.out, est[1,], settings, Nboot, good.boots, final, model)
    names(SE.out)<-c("standard.errors", "boot.base", "settings", "N.boots", "N.good.boots", "boot.results", "model")

    ## Ouputting result datasets
    dput(SE.out, paste0(out.dir,"/", out.name, ".Rdata", sep=""))

    return(SE.out)

  } # End of "ABfixselectUU" if statement







  ###################################
  ########### ABneutral      ########
  ###################################

  if (model == "ABneutral.R")
  {
  ## Reading the dataset for bootstrapping and extracting the parameter settings
    settings<-pedigree.data$settings
    est<-pedigree.data$estimates
    eqp<-as.numeric(as.character(settings[which(settings[,1] == "eqp"),2]))
    eqp.weight<-as.numeric(as.character(settings[which(settings[,1] == "eqp.weight"),2]))
    optim.method<-as.character(settings[which(settings[,1] == "optim.method"),2])
    p0uu<-round(as.numeric(as.character(settings[which(settings[,1] == "p0uu"),2])),16)
    pedigree<-pedigree.data$pedigree
    p0mm = 1-p0uu
    p0um = 0


  if(sum(c(p0mm, p0um, p0uu), na.rm =TRUE) != 1) {stop("The initial state probabilities don't sum to 1")}

  ##### Defining the divergence function
  divergence <- function(pedigree, p0mm, p0um, p0uu, param)
  {

    ## Initializing parameters
    PrMM <- p0mm
    PrUM <- p0um
    PrUU <- p0uu
    alpha <- param[1]
    bet <- param[2]
    weight <- param[3]


    ## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
    svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)



    ## Defining the generation (or transition) matrix
    Genmatrix <- matrix(c((1-alpha)^2, 2*(1-alpha)*alpha, alpha^2,
                          1/4*(bet + 1 - alpha)^2, 1/2*(bet + 1 - alpha)*(alpha + 1 - bet), 1/4*(alpha + 1 - bet)^2,
                          bet^2, 2*(1-bet)*bet, (1-bet)^2), nrow=3, byrow=TRUE)


    ## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
    Dt1t2<-NULL

    for (p in seq_len(NROW(pedigree)))
    {

      ## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
      svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
      svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
      svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
      svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

      ## Conditional divergences
      dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
                          svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

      dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
                          svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

      dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
                          svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

      ## Total (weighted) divergence
      Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


    }

    # Pr(UU) at equilibrium given alpha and beta
    puuinf.est<- (bet* ((1-bet)^2 - (1-alpha)^2 -1))/((alpha + bet)*((alpha + bet -1)^2 - 2))
    divout<-list(puuinf.est, Dt1t2)

    return(divout)

  }

  ###### Defining the Least Square function to be minimized
  ###### Note the equilibrium constraint, which can be made as small as desired.

  LSE_intercept<-function(param_int)
  {
     d = divergence(pedigree, p0mm, p0um, p0uu, param_int[1:3])
     sum((pedigree[,4] - param_int[4] - d[[2]])^2) +
      eqp.weight*nrow(pedigree)*((d[[1]]- eqp)^2)
  }



  ## Initializing
  final<-NULL
  counter<-0


  ## Defining starting values
  alpha.start  <-est[1,1]
  beta.start   <-est[1,2]
  weight.start <-est[1,3]
  intercept.start <-est[1,4]
  param_int0 = c(alpha.start, beta.start, weight.start, intercept.start)


  ## Start of boostrap loops
  for (booting in seq_len(Nboot))
  {

    opt.out<-NULL
    pedigree[,"div.obs"]<-pedigree[,"div.pred"]+sample(pedigree[,"residual"], nrow(pedigree), replace=TRUE)

    counter<-counter+1

    message("Bootstrap interationssss: ", counter/Nboot, "\n")


    opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
    alphafinal<-opt.out[1]
    betfinal<-opt.out[2]
    PrMMinf <- (alphafinal* ((1-alphafinal)^2 - (1-betfinal)^2 -1))/((alphafinal + betfinal)*((alphafinal + betfinal -1)^2 - 2))
    PrUMinf <- (4*alphafinal*betfinal*(alphafinal + betfinal -2))/((alphafinal + betfinal)*((alphafinal + betfinal -1)^2 -2))
    PrUUinf <- (betfinal* ((1-betfinal)^2 - (1-alphafinal)^2 -1))/((alphafinal + betfinal)*((alphafinal + betfinal -1)^2 - 2))
    opt.out <-cbind(opt.out, PrMMinf, PrUMinf, PrUUinf, alpha.start, beta.start, weight.start, intercept.start)
    final<-rbind(final, opt.out)


  } # End of Bootstrap loops


  colnames(final)[1:5]<-c("alpha", "beta", "weight", "sel.coef", "intercept")
  colnames(final)[13:15]<-c("PrMMinf", "PrUMinf", "PrUUinf")

  SE.alpha<-sd(final[,1],na.rm=TRUE)
  SE.beta<-sd(final[,2],na.rm=TRUE)
  SE.beta.alpha<-sd(final[,2]/final[,1], na.rm=TRUE)
  SE.weight<-sd(final[,3],na.rm=TRUE)
  SE.inter<-sd(final[,4],na.rm=TRUE)
  SE.PrMMinf<-sd(final[,13],na.rm=TRUE)
  SE.PrUMinf<-sd(final[,14],na.rm=TRUE)
  SE.PrUUinf<-sd(final[,15],na.rm=TRUE)

  CI.alpha<-quantile(final[,1],probs=c(0.025, 0.975))
  CI.beta<-quantile(final[,2],probs=c(0.025, 0.975))
  CI.beta.alpha<-quantile(final[,2]/final[,1], probs=c(0.025, 0.97))
  CI.weight<-quantile(final[,3],probs=c(0.025, 0.975))
  CI.inter<-quantile(final[,4],probs=c(0.025, 0.975))
  CI.PrMMinf<-quantile(final[,13],probs=c(0.025, 0.975))
  CI.PrUMinf<-quantile(final[,14],probs=c(0.025, 0.975))
  CI.PrUUinf<-quantile(final[,15],probs=c(0.025, 0.975))

  SE<-c(SE.alpha, SE.beta, SE.beta.alpha, SE.weight, SE.inter, SE.PrMMinf, SE.PrUMinf, SE.PrUUinf)
  CI<-rbind(CI.alpha, CI.beta, CI.beta.alpha, CI.weight, CI.inter, CI.PrMMinf, CI.PrUMinf, CI.PrUUinf)

  SE.out<-cbind(SE, CI)
  colnames(SE.out)[1]<-"SE"
  rownames(SE.out)<-c("alpha", "beta", "beta/alpha", "weight", "intercept", "PrMMinf", "PrUMinf", "PrUUinf")

  final<-data.frame(final)
  good.boots<-length(which(is.na(final[,"alpha"]) == FALSE))

  SE.out<-list(SE.out, est[1,], settings, Nboot, good.boots, final, model)
  names(SE.out)<-c("standard.errors", "boot.base", "settings", "N.boots", "N.good.boots", "boot.results", "model")

  ## Ouputting result datasets
  dput(SE.out, paste0(out.dir,"/", out.name, ".Rdata", sep=""))

  return(SE.out)

} # End of "ABneutral if statement







  ###################################
  ########### ABneutralSOMA     #####
  ###################################

 if (model == "ABneutralSOMA.R")
  {

allow.neg.intercept="yes"


  ##### Reading the dataset for bootstrapping and extracting the parameter settings
    settings<-pedigree.data$settings
    est<-pedigree.data$estimates
    eqp<-as.numeric(as.character(settings[which(settings[,1] == "eqp"),2]))
    eqp.weight<-as.numeric(as.character(settings[which(settings[,1] == "eqp.weight"),2]))
    optim.method<-as.character(settings[which(settings[,1] == "optim.method"),2])
    p0uu<-round(as.numeric(as.character(settings[which(settings[,1] == "p0uu"),2])),16)
    pedigree<-pedigree.data$pedigree
    p0mm = 1-p0uu
    p0um = 0
  #props<-pedigree.data$prop.data.used


  if(sum(c(p0mm, p0um, p0uu), na.rm =TRUE) != 1) {stop("The initial state probabilities don't sum to 1")}





##### Defining the divergence function
  divergence <- function(pedigree, p0mm, p0um, p0uu, param)
  {

    ## Initializing parameters
    PrMM <- p0mm
    PrUM <- p0um
    PrUU <- p0uu
    alpha <- param[1]
    bet <- param[2]
    weight <- param[3]


## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
    svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)



  ## Defining the generation (or transition) matrix for the mitotic case
    Genmatrix <- matrix(c((1-alpha)^2, 2*(1-alpha)*alpha,alpha^2,
                          bet*(1-alpha), (1-alpha)*(1-bet)+alpha*bet, alpha*(1-bet),
              bet^2, 2*(1-bet)*bet, (1-bet)^2),nrow=3, byrow=TRUE)


  ## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
    Dt1t2<-NULL

      for (p in seq_len(NROW(pedigree)) )
      {

      ## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
      svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree[p,1]))
      svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
      svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))
      svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,2] - pedigree[p,1]))
      svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree[p,3] - pedigree[p,1]))

      ## Conditional divergences
      dt1t2.MM  <- 1/2*(svt1.MM[,1] * svt2.MM[,2] + svt1.MM[,2] * svt2.MM[,1] + svt1.MM[,2] * svt2.MM[,3] +
                svt1.MM[,3] * svt2.MM[,2]) + 1*(svt1.MM[,1] * svt2.MM[,3]  + svt1.MM[,3] * svt2.MM[,1])

      dt1t2.UM  <- 1/2*(svt1.UM[,1] * svt2.UM[,2] + svt1.UM[,2] * svt2.UM[,1] + svt1.UM[,2] * svt2.UM[,3] +
                svt1.UM[,3] * svt2.UM[,2]) + 1*(svt1.UM[,1] * svt2.UM[,3] +  svt1.UM[,3] * svt2.UM[,1])

      dt1t2.UU  <- 1/2*(svt1.UU[,1] * svt2.UU[,2] + svt1.UU[,2] * svt2.UU[,1] + svt1.UU[,2] * svt2.UU[,3] +
                svt1.UU[,3] * svt2.UU[,2]) + 1*(svt1.UU[,1] * svt2.UU[,3] + svt1.UU[,3] * svt2.UU[,1])

      ## Total (weighted) divergence
      Dt1t2[p]<- svt0[,1]*dt1t2.UU + svt0[,2]*dt1t2.UM + svt0[,3]*dt1t2.MM


      }

    # Pr(MM) at equilibrium given alpha and beta
    puuinf.est<-(bet^2)/((alpha+bet)^2)
    divout<-list(puuinf.est, Dt1t2)

    return(divout)

  }


  ###### Defining the Least Square function to be minimized
  ###### Note the equilibrium constraint, which can be made as small as desired.

  LSE_intercept<-function(param_int)
  {
    sum((pedigree[,4] - param_int[4] - divergence(pedigree, p0mm, p0um, p0uu, param_int[1:3])[[2]])^2) +
      eqp.weight*nrow(pedigree)*((divergence(pedigree, p0mm, p0um, p0uu, param_int[1:3])[[1]]- eqp)^2)
  }






##### Initializing
  final<-NULL
  counter<-0
  opt.out<-NULL

  ## Defining starting values
  alpha.start  <-est[1,1]
  beta.start   <-est[1,2]
  weight.start <-est[1,3]
  intercept.start <-est[1,4]
  param_int0 = c(alpha.start, beta.start, weight.start, intercept.start)


  for (booting in seq_len(Nboot))
  {
    pedigree[,"div.obs"]<-pedigree[,"div.pred"]+sample(pedigree[,"residual"], nrow(pedigree), replace=TRUE)

      ## Initializing
      counter<-counter+1

      message("Bootstrap interation: ", counter/Nboot, "\n")

            opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
            alphafinal<-opt.out[1]
            betfinal<-opt.out[2]
            PrMMinf <- (alphafinal^2)/((alphafinal+betfinal)^2)
            PrUMinf <- (2*alphafinal*betfinal)/((alphafinal+betfinal)^2)
            PrUUinf <- (betfinal^2)/((alphafinal+betfinal)^2)
            opt.out <-cbind(opt.out, PrMMinf, PrUMinf, PrUUinf, alpha.start, beta.start, weight.start, intercept.start)
            final<-rbind(final, opt.out)


  } # End of Bootstrap loops


   colnames(final)[1:4]<-c("alpha", "beta", "weight", "intercept")
   colnames(final)[13:15]<-c("PrMMinf", "PrUMinf", "PrUUinf")

   if (allow.neg.intercept == "yes")
   { index<-which(final["alpha"] > 0 & final["beta"] > 0 & final["convcode"] == 0)}

   if (allow.neg.intercept == "no")
   {index<-which(final["alpha"] > 0 & final["beta"] > 0 & final["intercept"] > 0 & final["convcode"] == 0)}

   final<-final[index,]

   SE.alpha<-sd(final[,1],na.rm=TRUE)
   SE.beta<-sd(final[,2],na.rm=TRUE)
   SE.beta.alpha<-sd(final[,2]/final[,1], na.rm=TRUE)
   SE.weight<-sd(final[,3],na.rm=TRUE)
   SE.inter<-sd(final[,4],na.rm=TRUE)
   SE.PrMMinf<-sd(final[,13],na.rm=TRUE)
   SE.PrUMinf<-sd(final[,14],na.rm=TRUE)
   SE.PrUUinf<-sd(final[,15],na.rm=TRUE)

   CI.alpha<-quantile(final[,1],probs=c(0.025, 0.975))
   CI.beta<-quantile(final[,2],probs=c(0.025, 0.975))
   CI.beta.alpha<-quantile(final[,2]/final[,1], probs=c(0.025, 0.97))
   CI.weight<-quantile(final[,3],probs=c(0.025, 0.975))
   CI.inter<-quantile(final[,4],probs=c(0.025, 0.975))
   CI.PrMMinf<-quantile(final[,13],probs=c(0.025, 0.975))
   CI.PrUMinf<-quantile(final[,14],probs=c(0.025, 0.975))
   CI.PrUUinf<-quantile(final[,15],probs=c(0.025, 0.975))

   SE<-c(SE.alpha, SE.beta, SE.beta.alpha, SE.weight, SE.inter, SE.PrMMinf, SE.PrUMinf, SE.PrUUinf)
   CI<-rbind(CI.alpha, CI.beta, CI.beta.alpha, CI.weight, CI.inter, CI.PrMMinf, CI.PrUMinf, CI.PrUUinf)

   SE.out<-cbind(SE, CI)
   colnames(SE.out)[1]<-"SE"
   rownames(SE.out)<-c("alpha", "beta", "beta/alpha", "weight", "intercept", "PrMMinf", "PrUMinf", "PrUUinf")

    final<-data.frame(final)
      good.boots<-length(which(is.na(final[,"alpha"]) == FALSE))

      SE.out<-list(SE.out, est[1,], settings, Nboot, good.boots, final, model)
      names(SE.out)<-c("standard.errors", "boot.base", "settings", "N.boots", "N.good.boots", "boot.results", "model")

      ## Ouputting result datasets
      dput(SE.out, paste0(out.dir,"/", out.name, ".Rdata", sep=""))

  return(SE.out)

  }

} #End of function








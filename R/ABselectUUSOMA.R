#' Model with selection against spontaneous loss of methylation (outselectUU)
#' @param pedigree.data pedigree data.
#' @param p0uu initial proportion of unmethylated cytosines.
#' @param eqp equilibrium proportion of unmethylated cytosines.
#' @param eqp.weight weight assigned to equilibrium function.
#' @param Nstarts iterations for non linear LSQ optimization.
#' @param out.dir output directory.
#' @param out.name output file name.
#' @import optimx
#' @import expm
#' @importFrom stats runif
#' @return ABneutralSoma data.
#' @export
#' @examples
#'## Get some toy data
#' file1 <- system.file("extdata/soma/","pedigreeSoma.csv", package="AlphaBeta")
#' pedigree <- as.matrix(read.table(file1, sep=",", header=TRUE, stringsAsFactors = FALSE))
#' p0uu_in <- 0.54755
#' eqp.weight <- 0.001
#' Nstarts <- 2
#' output.data.dir <- paste0( getwd(),"/")
#' out.name <- "ABselectUUSOMA_CG_estimates"
#' out <- ABselectUUSOMA(pedigree.data = pedigree,
#'                   p0uu=p0uu_in,
#'                   eqp=p0uu_in,
#'                   eqp.weight=eqp.weight,
#'                   Nstarts=Nstarts,
#'                   out.dir=output.data.dir,
#'                   out.name=out.name)
#'
#' summary(out)
#'



ABselectUUSOMA<-function(pedigree.data, p0uu, eqp, eqp.weight, Nstarts, out.dir, out.name)
{


 allow.neg.intercept="no"

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
      sel    <-param[4]


	## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
	  svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

	  element11<-(1-alpha)^2
	  element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
	  element13<-(alpha^2)*sel
	  rowtotal1<-element11 + element12 + element13

	  element21<-(bet*(1-alpha))
	  element22<-((1-alpha)*(1-bet)+alpha*bet)*(1/2*(1+sel))
	  element23<-alpha*(1-bet)*sel
	  rowtotal2<-element21 + element22 + element23

	  element31<-(bet^2)
	  element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
	  element33<-((1-bet))^2*sel
	  rowtotal3<-element31 + element32 + element33

	## Defining the generation (or transition) matrix
	  Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
	                        element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
	                        element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)



	## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
	  Dt1t2<-NULL

		  for (p in 1:nrow(pedigree))
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
	  puuinf.est<- t(svGzero)  %*% ((Genmatrix)%^% 10000)
	  puuinf.est<- puuinf.est[1,1]
  	  divout<-list(puuinf.est, Dt1t2)

	  return(divout)

	}


###### Defining the Least Square function to be minimized
###### Note the equilibrium constraint, which can be made as small as desired.

		LSE_intercept<-function(param_int)
		{
			sum((pedigree[,4] - param_int[5] - divergence(pedigree, p0mm, p0um, p0uu, param_int[1:4])[[2]])^2) +
			eqp.weight*nrow(pedigree)*((divergence(pedigree, p0mm, p0um, p0uu, param_int[1:4])[[1]]-eqp)^2)
		}



###### Calculating the initial proportions
###### We always assume that:
		# 1. p0mm is larger than actually observed. This means if p0um is available from measurements,
		#    we will just add it to p0mm.
		# 2. As a consequence of (1.) we also assume that p0um = 0.

		p0uu<-p0uu
		p0mm<-1-p0uu
		p0um<-0


   if(is.null(p0mm ==TRUE | is.null(eqp)==TRUE))
   {stop("Both eqp value AND p0mm have to be supplied")}

   if(sum(c(p0mm, p0um, p0uu), na.rm =TRUE) != 1)
  {stop("The initial state probabilities don't sum to 1")}




##### Initializing
	optim.method<-"Nelder-Mead"
	final<-NULL
	counter<-0
	opt.out<-NULL
	pedigree<-pedigree.data


		for (s in 1:Nstarts)
		{

			## Draw random starting values
			alpha.start  <-10^(runif(1, log10(10^-9), log10(10^-2)))
			beta.start   <-10^(runif(1, log10(10^-9), log10(10^-2)))
	   		weight.start <-runif(1,0,0.5)
	     	sel.start <-runif(1,0.1,1)
	    	intercept.start <-runif(1,0,max(pedigree[,4]))
	    	param_int0 = c(alpha.start, beta.start, weight.start, sel.start, intercept.start)

			## Initializing
			counter<-counter+1

			cat("Progress: ", counter/Nstarts, "\n")


						opt.out  <- suppressWarnings(optimx(par = param_int0, fn = LSE_intercept, method=optim.method))
						alphafinal<-opt.out[1]
						betfinal<-opt.out[2]
						alphafinal<-as.numeric(opt.out[1])
    					betfinal<-as.numeric(opt.out[2])
    					weightfinal<-as.numeric(opt.out[3])
    					selfinal<-as.numeric(opt.out[4])

    					## Calculating equilibrium frequencies based on the model estimates
    					svGzero   <- c(p0uu, (weightfinal)*p0mm, (1-weightfinal)*p0mm)


							  element11<-(1-alphafinal)^2
							  element12<-(2*(1-alphafinal)*alphafinal)*(1/2*(1+selfinal))
							  element13<-(alphafinal^2)*selfinal
							  rowtotal1<-element11 + element12 + element13

							  element21<-(betfinal*(1-alphafinal))
							  element22<-((1-alphafinal)*(1-betfinal)+alphafinal*betfinal)*(1/2*(1+selfinal))
							  element23<-alphafinal*(1-betfinal)*selfinal
							  rowtotal2<-element21 + element22 + element23

							  element31<-(betfinal^2)
							  element32<-(2*((1-betfinal))*betfinal)*(1/2*(1+selfinal))
							  element33<-((1-betfinal))^2*selfinal
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

		} # End of Nstarts loop

	 colnames(final)[1:5]<-c("alpha", "beta", "weight", "sel.coef", "intercept")
	colnames(final)[14:16]<-c("PrMMinf", "PrUMinf", "PrUUinf")




##### Calculating the least square of the first part of the minimized function
	 lsqpart<-NULL

	 for (l in 1:nrow(final))
	 {
			  PrMM <- p0mm
			  PrUM <- p0um
	          PrUU <- p0uu
			  alpha  <- final[l, "alpha"]
			  bet    <- final[l, "beta"]
			  weight <- final[l, "weight"]
			  sel <-final[l, "sel.coef"]
			  intercept<-final[l,"intercept"]


			## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
			  svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

			  element11<-(1-alpha)^2
			  element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
			  element13<-(alpha^2)*sel
			  rowtotal1<-element11 + element12 + element13

			  element21<-(bet*(1-alpha))
			  element22<-((1-alpha)*(1-bet)+alpha*bet)*(1/2*(1+sel))
			  element23<-alpha*(1-bet)*sel
			  rowtotal2<-element21 + element22 + element23

			  element31<-(bet^2)
			  element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
			  element33<-((1-bet))^2*sel
			  rowtotal3<-element31 + element32 + element33

			## Defining the generation (or transition) matrix
			  Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
			                        element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
			                        element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

			  ## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
			  Dt1t2<-NULL

				  for (p in 1:nrow(pedigree))
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


			 ## Calculating the least square part
			 lsqpart[l]<-sum((pedigree[,4] - intercept - Dt1t2)^2)
		}

	 ## Collecting results and filtering them
	 final<-cbind(final, lsqpart)
	 colnames(final)[ncol(final)]<-c("value.part")
	 final<-final[order(final[,"value"]),]
	 #index.1<-which(final["alpha"] > 0 & final["beta"] > 0 & final["convcode"] == 0 & final["sel.coef"] >= 0 & final["sel.coef"] <= 1)
	 index.1<-which(final["alpha"] > 0 & final["beta"] > 0 & final["sel.coef"] >= 0 & final["sel.coef"] <= 1 & final["intercept"] > 0)
	 index.2<-setdiff(1:nrow(final), index.1)
	 final.1<-final[index.1,]
	 final.2<-final[index.2,]


  ## Calculting the predicted values based on the 'best' model (i.e. that with the lowest least square)
	 PrMM <- p0mm
	 PrUM <- p0um
	 PrUU <- p0uu
	 alpha  <- final.1[1, "alpha"]
	 bet    <- final.1[1, "beta"]
	 weight <- final.1[1, "weight"]
	 intercept<-final.1[1,"intercept"]
	 sel<-final.1[1,"sel.coef"]

	## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
	 svGzero   <- c(PrUU, (weight)*PrMM, (1-weight)*PrMM)

	  		  element11<-(1-alpha)^2
			  element12<-(2*(1-alpha)*alpha)*(1/2*(1+sel))
			  element13<-(alpha^2)*sel
			  rowtotal1<-element11 + element12 + element13

			  element21<-(bet*(1-alpha))
			  element22<-((1-alpha)*(1-bet)+alpha*bet)*(1/2*(1+sel))
			  element23<-alpha*(1-bet)*sel
			  rowtotal2<-element21 + element22 + element23

			  element31<-(bet^2)
			  element32<-(2*((1-bet))*bet)*(1/2*(1+sel))
			  element33<-((1-bet))^2*sel
			  rowtotal3<-element31 + element32 + element33

			## Defining the generation (or transition) matrix
			  Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
			                        element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
			                        element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)

			  ## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
			  Dt1t2<-NULL
			  Residual<-NULL

				  for (p in 1:nrow(pedigree))
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

			 ## Calculating the least square part
			 Residual<-(pedigree[,4] - intercept - Dt1t2)



## Augmenting pedigree
	delta.t<-pedigree[,2] + pedigree[,3] - 2*pedigree[,1]
	pedigree<-cbind(pedigree, delta.t, Dt1t2 + intercept, Residual)
	colnames(pedigree)[c(4,5,6,7)]<-c("div.obs", "delta.t","div.pred", "residual")

## Making info about settings
	info<-c("p0mm", "p0um", "p0uu", "eqp", "eqp.weight", "Nstarts", "optim.method")
	info2<-c(p0mm, p0um, p0uu, eqp, eqp.weight, Nstarts, optim.method)
	info.out<-data.frame(info, info2)
	colnames(info.out)<-c("Para", "Setting")


## Generating theoretical fit

		## Reading in pedigree
			obs<-pedigree[,"div.obs"]
			dtime<-pedigree[,"delta.t"]

		## Reading in parameter estimates
			est <-final.1
			alpha <-as.numeric(est[1,1])
	    beta<-as.numeric(est[1,2])
		  weight<-as.numeric(est[1,3])
		  sel<-as.numeric(est[1,4])
	    intercept<-as.numeric(est[1,5])

		## Reading initial state vector
			settings<-info.out
			PrMM<-p0mm<-as.numeric(as.character(settings[1,2]))
			PrUM<-p0um<-as.numeric(as.character(settings[2,2]))
			PrUU<-p0uu<-as.numeric(as.character(settings[3,2]))
			time1<- seq(1,max(c(pedigree[,2], pedigree[,3])))
			time2<- seq(1,max(c(pedigree[,2], pedigree[,3])))
			time.out<-expand.grid(time1,time2)
			#time0<- rep(min(pedigree[,1]), nrow(time.out))
			time0<- rep(0, nrow(time.out))
			pedigree.new<-as.matrix(cbind(time0,time.out))
			pedigree.new<-cbind(pedigree.new, c(pedigree.new[,2] + pedigree.new[,3] - 2*pedigree.new[,1]))
			pedigree.new<-pedigree.new[!duplicated(pedigree.new[,4]), ]
			pedigree.new<-pedigree.new[,1:3]

		## State probabilities at G0; first element = PrUU, second element = PrUM, third element = PrMM
			svGzero   <- c(PrUU, PrMM*weight, (1-weight)*PrMM)

							alphafinal<-alpha
							betfinal<-beta
							selfinal<-sel
							interceptfinal<-intercept

							element11<-(1-alphafinal)^2
							element12<-(2*(1-alphafinal)*alphafinal)*(1/2*(1+selfinal))
							element13<-(alphafinal^2)*selfinal
							rowtotal1<-element11 + element12 + element13

							element21<-(betfinal*(1-alphafinal))
							element22<-((1-alphafinal)*(1-betfinal)+alphafinal*betfinal)*(1/2*(1+selfinal))
							element23<-alphafinal*(1-betfinal)*selfinal
							rowtotal2<-element21 + element22 + element23

							element31<-(betfinal^2)
							element32<-(2*((1-betfinal))*betfinal)*(1/2*(1+selfinal))
							element33<-((1-betfinal))^2*selfinal
							rowtotal3<-element31 + element32 + element33

							## Defining the generation (or transition) matrix
							  Genmatrix <- matrix(c(element11/rowtotal1, element12/rowtotal1, element13/rowtotal1,
							                        element21/rowtotal2,  element22/rowtotal2,  element23/rowtotal2,
							                        element31/rowtotal3,  element32/rowtotal3,  element33/rowtotal3), nrow=3, byrow=TRUE)


							## Calculating theoretical divergence for every observed pair in 'pedigree.txt'
							Dt1t2<-NULL

								for (p in 1:nrow(pedigree.new))
								{

									## Define state vectors for t1,t2 and t0 from pedigree using matrix multiplications from library(expm)
									svt0      <- t(svGzero)  %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,1]))
									svt1.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
									svt2.MM   <- t(c(0,0,1)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))
									svt1.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
									svt2.UM   <- t(c(0,1,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))
									svt1.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,2] - pedigree.new[p,1]))
									svt2.UU   <- t(c(1,0,0)) %*% ((Genmatrix)%^% as.numeric(pedigree.new[p,3] - pedigree.new[p,1]))

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

			pedigree.new<-cbind(pedigree.new, Dt1t2+interceptfinal, c(pedigree.new[,2] + pedigree.new[,3] - 2*pedigree.new[,1]))
			colnames(pedigree.new)<-c("time0", "time1", "time2", "div.sim", "delta.t")
			pedigree.new<-pedigree.new[order(pedigree.new[,5]),]

			model<-"ABselectUUSOMA.R"

    	abfreeS.out<-list(final.1, final.2, pedigree, info.out, model, pedigree.new)
    	names(abfreeS.out)<-c("estimates", "estimates.flagged", "pedigree", "settings", "model", "for.fit.plot")

## Ouputting result datasets
	dput(abfreeS.out, paste(out.dir, out.name, ".Rdata", sep=""))
	return(abfreeS.out)



} #End of function


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3992612.svg)](https://doi.org/10.5281/zenodo.3992612)


# ![jlab-code/AlphaBeta](vignettes/alphabeta.png)


### Computational inference of epimutation rates and spectra from high-throughput DNA methylation data in plants


------------------------------------------------------------------------

**AlphaBeta** is a computational method for estimating epimutation rates and spectra from high-throughput DNA methylation data in plants.

The method has been specifically designed to:

**1.** Analyze 'germline' epimutations in the context of multi-generational mutation accumulation lines (MA-lines).

**2.** Analyze 'somatic' epimutations in the context of plant development and aging.

Heritable changes in cytosine methylation can arise stochastically in plant genomes independently of DNA sequence alterations. 
These so-called ‘spontaneous epimutations’ appear to be a byproduct of imperfect DNA methylation maintenance during mitotic and meitotic cell divisions. 

Accurate estimates of the rate and spectrum of these stochastic events are necessary to be able to quantify how epimutational processes shape methylome diversity in the context of plant evolution, development and aging. 

Here we describe AlphaBeta, a computational method for estimating epimutation rates and spectra from pedigree-based high-throughput DNA methylation data in plants. 

The method requires that the topology of the pedigree is known, which is typically the case in the construction of mutation accumulation lines (MA-lines) in sexually or clonally reproducing plant species. 

However, the method also works for inferring somatic epimutations in long-lived perrenials, such as trees, using leaf methylomes and coring data as input. In this case, AlphaBeta treats the tree branching structure as an intra-organismal phylogeny of somatic lineages that carry information about the epimutational history of each branch.  

------------------------------------------------------------------------

### Installation


##### Installing from the Github

To install from GitHub (development version), follow the steps given below. 

##### Step 1 — Install a last version of R (>=3.6)

##### Step 2 — In R, please install all dependencies and execute the following commands:
 - 2.1.  install.packages("devtools")
 - 2.2.  library(devtools)
 - 2.3.  install_github("jlab-code/AlphaBeta")

------------------------------------------------------------------------

### How to use AlphaBeta

Please open the [vignette](https://github.com/jlab-code/AlphaBeta/blob/master/vignettes/AlphaBeta.pdf) file.

------------------------------------------------------------------------
### Contributors

- Yadi Shahryary  -  y.shahryary@tum.de
- Rashmi Hazarika -  hazarika.rr@gmail.com  
- Frank Johannes  -  frank@johanneslab.org 



rm(list=ls())
basedir <- '~/Dropbox/Neurodegeneration/TauSpread/tau-spread/'
setwd(basedir)
params <- list(basedir=basedir,               
               grps = c('NTG','G20'),
               #injection.site=c('iDG','iVISam'), # define injection sites (using ABA nomenclature)
               injection.site = c('iDG', 'iCA1', 'iCA3', 'iVISam', 'iRSPagl'), # new sites
               tps=c(1,3,6,9), # define time points post injection
               c.min = 1e-5, # empirically conservative minimum for time constant
               c.max = 0.2, # empirically conservative maximum for time constant
               c.n = 100) # number of time constants to try
source('code/misc/miscfxns.R')
params$source.save <- source.save
params$opdir <- paste('TauDiffusion032520_Inject',paste0(params$injection.site,collapse='-'),'_CMax',params$c.max,'/',sep='')
dir.create(params$opdir,recursive = T)

#################################################
### Load packages & create output directories ###
#################################################

source('code/misc/packages.R')

##############################
### Process pathology data ###
##############################

source('code/process/transpose_data.R')
source('code/process/process.R')
#source('code/process/testCNDRtoABAmatching.R')

########################################
### Process ABA gene expression data ###
########################################

source('code/aba/process_ontology.R')
source('code/aba/download_gene_expression.R') # specify variable gois in this file to get expression for whatever genes you want

#######################
### Diffusion model ###
#######################

# retrograde model with additive Mapt expression
injection.site <- params$injection.site # use all seed sites at once
for(grp in params$grps){
  source('code/diffmodel/analyzespread_CNDRspace.R')
  goi <- 'Mapt'
  probe <- 'RP_071204_01_D02'
  source('code/diffmodel/plotCNDRspacefit.R')
}

# retrograde model with additive Mapt expression
for(injection.site in params$injection.site){ # use each seed site separately
  for(grp in params$grps){
    source('code/diffmodel/analyzespread_CNDRspace.R')
    goi <- 'Mapt'
    probe <- 'RP_071204_01_D02'
    source('code/diffmodel/plotCNDRspacefit.R')
  }
}


# bidirectional, independent, additive diffusion model: anterograde and retrograde additive and independent
injection.site <- params$injection.site # use all seed sites at once
for(grp in params$grps){
  source('code/diffmodel/analyzebidirectionalspread_CNDRspace.R')
  source('code/diffmodel/plotCNDRspacebidirectionalfit.R')
}

for(injection.site in params$injection.site){ # use each seed site separately
  for(grp in params$grps){
    source('code/diffmodel/analyzebidirectionalspread_CNDRspace.R')
    source('code/diffmodel/plotCNDRspacebidirectionalfit.R')
  }
}

for(grp in params$grps){
  source('code/diffmodel/seedspec.R')
  source('code/diffmodel/plotseedspec.R')
}

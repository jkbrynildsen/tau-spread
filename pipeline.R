rm(list=ls())
basedir <- '~/Box/PD_Project/Quant_pathology_dataset/'
setwd(basedir)
params <- list(basedir=basedir,
               grps = c('NTG', 'G20'), # NTG are WT mice; G20 are mice with LRRK2 mutation
               injection.site = c('iDG', 'iCA1', 'iCA3', 'iVISam', 'iRSPagl'), # injection sites (using ABA nomenclature)
               tps=c(3,6), # define time points post injection
               treatments=c('0 MLi-2','75 MLi-2', '450 MLi-2'), # treatment group - 0, 75, or 450 MLi-2
               c.min = 1e-10, # minimum for time constant
               c.max = 0.001, # maximum for time constant
               c.n = 100) # number of time constants to try
injection.site <- params$injection.site
source('code/misc/miscfxns.R')
params$source.save <- source.save
params$opdir <- paste('TauDiffusion030524_Inject',paste0(params$injection.site,collapse='-'),'_CMax',params$c.max,'/',sep='')
dir.create(params$opdir,recursive = T)

#################################################
### Load packages & create output directories ###
#################################################

source('code/misc/packages.R')

#############################################
### Process pathology and structural data ###
#############################################

source('code/process/pathology_data_prep.R')
source('code/process/process_struct.R')

#######################
### Diffusion model ###
#######################

# separately fit retrograde and bidirectional diffusion models for each genotype and treatment group
injection.site <- params$injection.site
for(grp in params$grps){
  for(treatment in params$treatments){
    source('code/diffmodel/analyze_retrogradespread.R')
    source('code/diffmodel/analyzebidirectionalspread.R') # run this to initialize parameters for optimization
    source('code/diffmodel/plot_fit.R')
}}

# bidirectional diffusion model with different a single linear model weighting anterograde and retrograde for each genotype and treatment group
for(grp in params$grps){
  for(treatment in params$treatments){
    source('code/diffmodel/optim_bidirectionalspread_onelm.R')
    source('code/diffmodel/plot_bidirectionalonelmfit.R')
}}

###################################################
### Get euclidean distances between ABA regions ###
###################################################

source('code/aba/process_ontology.R')
source('code/aba/atlas_structures.R') # required for most scripts in code/nullmodels/

################################
### Quality control analyses ###
################################
# these are very time consuming
# for the manuscript, model comparisons were performed on pooled data since results were comparable across genotype and treatment groups

# compare model fits based on random vs actual seed regions      
injection.site <- params$injection.site
# refit bidirectional model on all data (pooled across genotypes and treatment groups)
source('code/diffmodel/analyzebidirectionalspread_all.R')
source('code/diffmodel/optim_bidirectionalspread_onelm_all.R')
# compare fit to random seeds
source('code/nullmodels/seedspec_multi_onelm.R')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               source('code/nullmodels/seedspec_multi_onelm.R')
source('code/nullmodels/plotseedspec_multi.R')

# compare out of sample fits between models (anterograde, retrograde, euclidean, bidirectional) on pooled data
injection.site <- params$injection.site; 
source('code/modelcomparison/modelcomparison_traintest.R')
source('code/modelcomparison/plot_modelcomparison_testset.R')

############################
### Bootstrapping analyses to compare model parameters across groups ###
############################

# bootstrap models for treatment groups to compare parameters
for(grp in params$grps){
  for (treatment in params$treatments){
    source('code/modelcomparison/bootstrap_optimspread_bidirectional.R')
  }}

for(grp in params$grps){
  group1 <- paste0(grp,"_0 MLi-2")
  group2 <- paste0(grp,"_75 MLi-2")
  group3 <- paste0(grp,"_450 MLi-2")
  source('code/modelcomparison/Groupwise_bootstrap_onelm_3grp_comparison.R')
}
source('code/modelcomparison/Groupwise_bootstrap_onelm_3grp_comparison.R')

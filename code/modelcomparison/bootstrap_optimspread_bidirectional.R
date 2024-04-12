#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','treatment','injection.site')))

print(grp)
print(treatment)
basedir <- params$basedir
setwd(basedir)
savedir <- paste0(params$opdir,'modelcomparison/bootstrap/',paste0(injection.site,collapse='-'),'/')
dir.create(savedir,recursive=T)

source('code/misc/modelcomparisonfxns.R')
load(paste0(params$opdir,'processed/pathdata.RData'))  # load path data and ROI names

# get mean pathology for each time point
tps = 6 # use only data from 6 MPI to maintain consistency across groups
Mice <- lapply(tps, function(tp) path.data[path.data$Genotype == grp & path.data$Treatment == treatment & path.data$MPI == tp, path.names])
# bootstrap log pathology
nboot <- 500
log.path.boot <- bootstrap.path.tps(Mice,nboot)

# load structural networks
W <- read.csv(paste(params$opdir,'processed/W.csv', sep=''), header = FALSE) # subbed this and the below line for the above line
W <- as.matrix(W)
L.out.retro <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='retro') # compute out-degree Laplacian for retrograde connectivity only (not weighted by Snca)
L.out.antero <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='antero') # compute out-degree Laplacian for anterograde connectivity only (not weighted by Snca)

#########################################################################################
### Use optimization to jointly fit two time scaling parameter on average of all mice ###
#########################################################################################

# set up controls for all models
Xo <- get.Xo(region.names,injection.site) # seed pathology
load(file=paste0(params$opdir,'diffmodel/bidirectional/',paste0(injection.site,collapse='-'),'_independentfit/',grp,'_',treatment,'IndependentBidirectionalFit_params.RData'))
ctrl.optim <- list(fnscale=-1) # set controls for optim function -- maximize the objective function instead of minimizing (default)
ctrl <- list(Retrograde=L.out.retro,Anterograde=L.out.antero,
             Xo=Xo,fxn=reticulate::import('scipy.sparse.linalg')$expm,tps=tps,
             c.r.a.init=c(c.Grp.retro,c.Grp.antero),
             ctrl.optim=ctrl.optim,
             one.lm=FALSE)

# Loop through model types, train on train sets, evaluate on train and test sets
mdl.names <- c('BidirectionalOneLM')
results <- list()
for(mdl.name in mdl.names){
  print(paste('training and testing',mdl.name))
  results[[mdl.name]] <- list(train=NULL)
  for(REP in 1:nboot){
    print(paste('REP',REP))
    X.train <- log.path.boot[[REP]]
    m.out <- spread.model.train(X.train,mdl.name,ctrl)
    results[[mdl.name]]$train[[REP]] <- spread.model.eval(m.out,X.train)
  }
}

# save data
save(results,Mice,file = paste0(savedir,grp,'_',treatment,'ModelComparison_Bootstrap.RData'))

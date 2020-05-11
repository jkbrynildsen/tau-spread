#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','injection.site')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/bidirectional_bootstrap/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/optimfxns.R')
source('code/misc/miscfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

# get mean pathology for each time point
tps <- params$tps
Mice <- lapply(tps, function(tp) path.data[path.data$Condition == grp & path.data$Month == tp,path.names])
# bootstrap log pathology
log.path.boot <- bootstrap.path.tps(Mice,nboot = 500)

# load structural networks
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
L.out.retro <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='retro') # compute out-degreee Laplacian for retrograde connectivity only (not weighted by Snca)
L.out.antero <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='antero') # compute out-degreee Laplacian for anterograde connectivity only (not weighted by Snca)

#########################################################################################
### Use optimization to jointly fit two time scaling parameter on average of all mice ###
#########################################################################################

Xo <- get.Xo(region.names,injection.site) # seed pathology in iCPu
scipy.linalg <- reticulate::import('scipy.linalg') # import scipy matrix exponential function because it's faster
#params.opt <- c(0.006070303,0.02223111) # c.retro, c.antero: use values from independent fit to initialize
#params.opt <- c(0.01,0.01) # c.retro, c.antero: use values from independent fit to initialize
load(file=paste0(params$opdir,'diffmodel/bidirectional/',paste0(injection.site,collapse='-'),'_independentfit/',grp,'CNDRSpaceIndependentBidirectionalFit_params.RData'))
params.opt <- c(c.Grp.retro,c.Grp.antero)
ctrl <- list(fnscale=-1) # maximize the objective function instead of minimizing (default) # c's must be > 0

results <- list()
nboot <- length(log.path.boot)
df.resid.init <- data.frame(matrix(ncol=length(tps),nrow=n.regions.CNDR, dimnames=list(path.names, paste(tps,'MPI'))),check.names = FALSE) # initialize data frame to hold residuals for every region (even if some are NaN b/c of 0 pathology)
for(REP in 1:nboot){
  print(paste('REP',REP))
  X.train <- log.path.boot[[REP]] # optimize model on training data
  params.opt.fit <- optim(params.opt,c.CNDRspace.objective,control = ctrl,method="L-BFGS-B",lower=c(1e-7,1e-7), # optimization. c's must be > 0
                          log.path=X.train,tps=tps,L.out.retro=L.out.retro,L.out.antero=L.out.antero,
                          Xo=Xo,ABA.to.CNDR.key=ABA.to.CNDR.key,fxn =scipy.linalg$expm) # static inputs  
  # extract parameters from optimization output
  c.train.retro <- params.opt.fit$par[1]
  c.train.antero <- params.opt.fit$par[2]
  
  # retrieve linear models for antero and retro
  Xt.Grp.retro <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.retro,Xo,c.train.retro,t,fxn=scipy.linalg$expm),path.names,ABA.to.CNDR.key)), base = 10))) # predict pathology using connectivity, time constant, and seed
  Xt.Grp.antero <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.antero,Xo,c.train.antero,t,fxn=scipy.linalg$expm),path.names,ABA.to.CNDR.key)), base = 10))) # predict pathology using connectivity, time constant, and seed
  df.train <- lapply(1:length(tps), function(t) data.frame(path = X.train[[t]], pred.retro = Xt.Grp.retro[,t,drop=FALSE], pred.antero = Xt.Grp.antero[,t,drop=FALSE]))
  m.train <- lapply(df.train, function(df.i) lm(path~pred.retro+pred.antero,data=inf.nan.mask(df.i))) # fit linear regression at each time point to combine anterograde and retrograde
  train.fits.r <- sapply(m.train, function(m.i) cor(m.i$fitted.values,m.i$model$path)) # extract fits as pearson r
  train.fits.mse <- sapply(m.train, function(m.i) mean(residuals(m.i)^2))
  resid <- sapply(m.train, function(m.i) residuals(m.i)) # get residuals of each model
  df.resid <- df.resid.init
  for(tp in 1:length(tps)){df.resid[names(resid[[tp]]),paste(tps[tp],'MPI')] <- resid[[tp]]}
  results[[REP]] <- list(c.train.retro=c.train.retro,c.train.antero=c.train.antero,resid=df.resid,
                         m.train=m.train,train.fits.r=train.fits.r,train.fits.mse=train.fits.mse)
}

save(results,file = paste(savedir,grp,'CNDRSpaceBidirectionalOptim_Bootstrap.RData',sep=''))

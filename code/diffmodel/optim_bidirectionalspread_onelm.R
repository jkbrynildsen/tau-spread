#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','treatment','injection.site')))
print(grp)
print(treatment)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/bidirectional_onelm/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/optimfxns.R')
source('code/misc/miscfxns.R')
load(paste0(params$opdir,'processed/pathdata.RData'))  # load path data and ROI names

# get mean pathology for each time point
if (treatment == '75 MLi-2'){tps = 6} else {tps=c(3,6)} # assign time points (3 and 6 mo for 0 and 450 MLi-2 treatment; 6 mo only for 75 MLi-2 treatment)
Mice <- lapply(tps, function(tp) path.data[path.data$Genotype == grp & path.data$Treatment == treatment & path.data$MPI == tp,path.names])
Grp.mean <- lapply(Mice,function(X) colMeans(X,na.rm = T))
log.path <- lapply(Grp.mean, function(X) log(X,base=10))

# load structural networks
W <- read.csv(paste(params$opdir,'processed/W.csv', sep=''), header = FALSE)
W <- as.matrix(W)

L.out.retro <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='retro') # compute out-degreee Laplacian for retrograde connectivity only
L.out.antero <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='antero') # compute out-degreee Laplacian for anterograde connectivity only

#########################################################################################
### Use optimization to jointly fit two time scaling parameter on average of all mice ###
#########################################################################################

Xo <- get.Xo(region.names,injection.site) # seed pathology
scipy.sparse.linalg <- reticulate::import('scipy.sparse.linalg') # import scipy matrix exponential function because it's faster
load(file=paste0(params$opdir,'diffmodel/bidirectional/',paste0(injection.site,collapse='-'),'_independentfit/',grp,'_',treatment,'IndependentBidirectionalFit_params.RData'))
params.opt <- c(c.Grp.retro,c.Grp.antero)
ctrl <- list(fnscale=-1) # minimize objective function

params.opt.fit <- optim(params.opt,c.objective,control = ctrl, method = 'L-BFGS-B',lower=c(10e-7,10e-7), # optimization. c's must be > 0
                        log.path=log.path,tps=tps,L.out.retro=L.out.retro,L.out.antero=L.out.antero,
                        Xo=Xo,fxn=scipy.sparse.linalg$expm,one.lm=TRUE,excl.inj=NULL) # static inputs

# extract time constants
c.Grp.retro <- params.opt.fit$par[1]
c.Grp.antero <- params.opt.fit$par[2]

# save data
Xt.Grp.retro <- do.call('cbind',lapply(tps, function(t) log(quiet(predict.Lout(L.out.retro,Xo,c.Grp.retro,t,fxn=scipy.sparse.linalg$expm)),base=10))) # predict pathology using connectivity, time constant, and seed
Xt.Grp.antero <- do.call('cbind',lapply(tps, function(t) log(quiet(predict.Lout(L.out.antero,Xo,c.Grp.antero,t,fxn=scipy.sparse.linalg$expm)),base=10))) # predict pathology using connectivity, time constant, and seed
lm.mask.results <- lm.mask.ant.ret.all(log.path,10^Xt.Grp.retro,10^Xt.Grp.antero) # undo log10 because this function automatically computes it
m <- lm.mask.results$m; e <- lm.mask.results$e; m.fits <- lm.mask.results$m.fits; df <- lm.mask.results$df
save(df,c.Grp.antero,c.Grp.retro,m,m.fits,file = paste(savedir,grp,'_',treatment,'BidirectionalOptimOneLM_data.RData',sep=''))

#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','treatment','injection.site')))
print(grp)
print(treatment)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/retrograde/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
load(paste0(params$opdir,'processed/pathdata.RData'))  # load path data and ROI names

# get mean pathology for each time point
if (treatment == '75 MLi-2'){tps = 6} else {tps=c(3,6)} # assign time points (3 and 6 mo for 0 and 450 MLi-2 treatment; 6 mo only for 75 MLi-2 treatment)
Mice <- lapply(tps, function(tp) path.data[path.data$Genotype == grp & path.data$Treatment == treatment & path.data$MPI == tp, path.names])
Grp.mean <- lapply(Mice,function(X) colMeans(X,na.rm = T))

W <- read.csv(paste(params$opdir,'processed/W.csv', sep=''), header = FALSE)
W <- as.matrix(W)

L.out <- get.Lout(W,rep(1,n.regions.ABA)) # compute out-degree Laplacian for connectivity only

# Fit time scaling parameter on average of all mice
c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n) # scaling parameter
log.path <- lapply(Grp.mean, function(X) log(X,base=10))
Xo <- get.Xo(region.names,injection.site) # seed pathology in injection sites
c.output <- c.fit(log.path,tps,L.out,Xo,c.rng,excl.inj = params$excl.inj)
c.Grp <- c.output$c.best; Xt.sweep <- c.output$Xt.sweep # store time constant and r values

Xt.Grp <- do.call('cbind',lapply(tps, function(t) log(quiet(predict.Lout(L.out,Xo,c.Grp,t)), base = 10))) # predict pathology using connectivity, time constant, and seed
df <- lapply(1:length(tps), function(t) data.frame(path = log.path[[t]], pred = Xt.Grp[,t,drop=FALSE]))
save(df,c.Grp,Xt.sweep,file = paste(savedir,grp,'_',treatment,'RetroFit_data.RData',sep=''))
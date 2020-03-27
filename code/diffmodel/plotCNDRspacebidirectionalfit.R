#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','injection.site')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/bidirectional/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
source('code/misc/plottingfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
load(paste(savedir,grp,'CNDRSpaceBidirectionalFit_data.RData',sep=''))
tps <- params$tps

# exclude regions with 0 pathology at each time point for purposes of computing fit
lapply(1:length(tps), function(t) paste0(t,' MPI: ',sum(df[[t]]$path != -Inf & !is.na(df[[t]]$pred)),'/',nrow(df[[t]]),' regions left'))

##########################################################
### add together anterograde and retrograde prediction ###
##########################################################

m <- lapply(df, function(df.i) m <- lm(path~pred.retro+pred.antero,data=inf.nan.mask(df.i)))
vulnerability <- lapply(m,residuals)
pred <- lapply(m, function(m.i) m.i$fitted.values)

p <- lapply(1:length(tps), function(t) 
  p.xy(x=m[[t]]$fitted.values,y=m[[t]]$model$path,ylab='Actual',xlab='Predicted',
       ttl=paste0(grp,': ',tps[t],' MPI'),col='#007257',alpha=0.7))
p <- plot_grid(plotlist=p,align='hv',nrow=1)
ggsave(p,filename = paste(savedir,grp,'CNDRSpaceFit_bidirectionaladditivemodel.pdf',sep=''),
       units = 'cm',height = 3.75,width = 3.75*length(tps))
#lapply(1:length(tps),function(t) write.csv(data.frame(vuln=vulnerability[[t]]),paste0(savedir,grp,'vulnerability',tps[t],'MPI',goi,probe,'.csv')))
#lapply(1:length(tps),function(t) write.csv(data.frame(pred=pred[[t]]),paste0(savedir,grp,'predictedpath',tps[t],'MPI',goi,probe,'.csv')))

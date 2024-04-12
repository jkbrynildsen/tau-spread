#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','treatment','injection.site')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/bidirectional/',paste0(injection.site,collapse='-'),'_independentfit/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
source('code/misc/plottingfxns.R')
load(paste0(params$opdir,'processed/pathdata.RData'))  # load path data and ROI names
load(paste0(savedir,grp,'_',treatment,'BidirectionalFit_data.RData'))
if (treatment == '75 MLi-2'){tps = 6} else {tps=c(3,6)} # assign time points (3 and 6 mo for 0 and 450 MLi-2 treatment; 6 mo only for 75 MLi-2 treatment)

###########################
### Plot retrograde fit ###
###########################

# exclude regions with 0 pathology at each time point for purposes of computing fit
lapply(1:length(tps), function(t) paste0(t,' MPI: ',sum(df[[t]]$path != -Inf & !is.na(df[[t]]$pred.retro)),'/',nrow(df[[t]]),' regions left'))
# plot for each time point, using p.xy function 
p <- lapply(1:length(tps), function(t) 
  p.xy(x=df[[t]]$pred.retro,y=df[[t]]$path,ylab='Log(pathology)',xlab='Log(predicted)',
       ttl=paste0(grp,': ',tps[t],' MPI'),col='#007257',alpha=0.7))
p <- plot_grid(plotlist=p,align='hv',nrow=1)

ggsave(p,filename = paste(savedir,grp,'_',treatment,'retroFit_basemodel.pdf',sep=''),
       units = 'cm',height = 3.75,width = 3.75*length(tps))

c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n)
p <- plot.Xt(Xt = t(Xt.sweep.retro),t. = c.rng) + geom_vline(xintercept = c.Grp.retro,color='grey50',linetype='dashed') +
  annotate(geom='text',x=c.Grp.retro,y=max(Xt.sweep.retro),hjust=-0.25,label='Optimal',size=2,color='grey50')+
  scale_color_manual(values = as.character(1:length(tps)),labels=paste(tps,'MPI'),name='') + 
  xlab('c') + ylab('Pearson r with\nPathology') +ggtitle(paste0(grp,': ',paste0(injection.site,collapse = '-'))) +
  theme_bw() + theme(text=element_text(size=8),plot.title = element_text(size=8,hjust=0.5),legend.key.size = unit(0.1,'cm'))
p
ggsave(p,filename = paste(savedir,grp,'_',treatment,'retroCSweepByTimePoint_basemodel.pdf',sep=''),
       units = 'cm',height = 4,width = 9)

############################
### Plot anterograde fit ###
############################
savedir <- paste(params$opdir,'diffmodel/anterograde/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

# exclude regions with 0 pathology at each time point for purposes of computing fit
lapply(1:length(tps), function(t) paste0(t,' MPI: ',sum(df[[t]]$path != -Inf & !is.na(df[[t]]$pred.antero)),'/',nrow(df[[t]]),' regions left'))
# plot for each time point, using p.xy function 
p <- lapply(1:length(tps), function(t) 
  p.xy(x=df[[t]]$pred.antero,y=df[[t]]$path,ylab='Log(pathology)',xlab='Log(predicted)',
       ttl=paste0(grp,': ',tps[t],' MPI'),col='#007257',alpha=0.7))
p <- plot_grid(plotlist=p,align='hv',nrow=1)

ggsave(p,filename = paste(savedir,grp,'_',treatment,'anteroFit_basemodel.pdf',sep=''),
       units = 'cm',height = 3.75,width = 3.75*length(tps))

c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n)
p <- plot.Xt(Xt = t(Xt.sweep.antero),t. = c.rng) + geom_vline(xintercept = c.Grp.antero,color='grey50',linetype='dashed') +
  annotate(geom='text',x=c.Grp.retro,y=max(Xt.sweep.antero),hjust=-0.25,label='Optimal',size=2,color='grey50')+
  scale_color_manual(values = as.character(1:length(tps)),labels=paste(tps,'MPI'),name='') + 
  xlab('c') + ylab('Pearson r with\nPathology') +ggtitle(paste0(grp,': ',paste0(injection.site,collapse = '-'))) +
  theme_bw() + theme(text=element_text(size=8),plot.title = element_text(size=8,hjust=0.5),legend.key.size = unit(0.1,'cm'))
p
ggsave(p,filename = paste(savedir,grp,'_',treatment,'anteroCSweepByTimePoint_basemodel.pdf',sep=''),
       units = 'cm',height = 4,width = 9)
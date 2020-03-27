#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','goi','probe','injection.site')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
source('code/misc/plottingfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
load(paste(savedir,grp,'CNDRSpaceFit_data.RData',sep=''))
tps <- params$tps

# exclude regions with 0 pathology at each time point for purposes of computing fit
lapply(1:length(tps), function(t) paste0(t,' MPI: ',sum(df[[t]]$path != -Inf & !is.na(df[[t]]$pred)),'/',nrow(df[[t]]),' regions left'))
# plot for each time point, using p.xy function 
p <- lapply(1:length(tps), function(t) 
  p.xy(x=df[[t]]$pred,y=df[[t]]$path,ylab='Actual',xlab='Predicted',
       ttl=paste0(grp,': ',tps[t],' MPI'),col='#007257',alpha=0.7))
p <- plot_grid(plotlist=p,align='hv',nrow=1)

ggsave(p,filename = paste(savedir,grp,'CNDRSpaceFit_basemodel.pdf',sep=''),
       units = 'cm',height = 3.75,width = 3.75*length(tps))

###############################################
### Save predicted values and vulnerability ###
###############################################

m <- lapply(df, function(df.i) m <- lm(path~pred,data=inf.nan.mask(df.i)))
vulnerability <- lapply(m,residuals)

lapply(1:length(tps),function(t) write.csv(data.frame(vuln=vulnerability[[t]]),paste0(savedir,grp,'vulnerability',tps[t],'MPI.csv')))
lapply(1:length(tps),function(t) write.csv(data.frame(pred=df[[t]][,'pred',drop=FALSE]),paste0(savedir,grp,'predictedpath',tps[t],'MPI.csv')))
                                                                             
#####################
### Add gene data ###
#####################

# load gene data specified by goi, probe input variables
gene.exp <- read.csv(paste0(basedir,'data/aba/expression/',goi,probe,'ExpressionCNDR.csv'),row.names = 1)
df.gene <- lapply(df, function(x) merge(x,gene.exp,by=0)) # merge
for(j in 1:length(df.gene)){rownames(df.gene[[j]]) <- df.gene[[j]]$Row.names} # remove row names
for(j in 1:length(df.gene)){df.gene[[j]] <- df.gene[[j]][,-1]}

##########################################################################
### show how gene expression correlates with residuals of spread model ###
##########################################################################

vuln.gene <- lapply(vulnerability, function(v) merge(as.data.frame(v),gene.exp,by=0))
p <- lapply(1:length(tps), function(t) 
  p.xy(x=vuln.gene[[t]][,paste0(goi,'_',probe)],y=vuln.gene[[t]]$v,ylab='Vulnerability',xlab=paste0(goi,' Expression'),
       ttl=paste0(grp,': ',tps[t],' MPI'),col='#007257',alpha=0.7))
p <- plot_grid(plotlist=p,align='hv',nrow=1)
ggsave(p,filename = paste(savedir,grp,'CNDRSpaceVulnerability_basemodel_vs',goi,probe,'.pdf',sep=''),
       units = 'cm',height = 3.75,width = 3.75*length(tps))

##############################
### Include genes in model ###
##############################

m <- lapply(df.gene, function(df.i) m <- lm(path~.,data=inf.nan.mask(df.i)))
vulnerability <- lapply(m,residuals)
pred <- lapply(m, function(m.i) m.i$fitted.values)

p <- lapply(1:length(tps), function(t) 
  p.xy(x=m[[t]]$fitted.values,y=m[[t]]$model$path,ylab='Actual',xlab='Predicted',
       ttl=paste0(grp,': ',tps[t],' MPI'),col='#007257',alpha=0.7))
p <- plot_grid(plotlist=p,align='hv',nrow=1)
ggsave(p,filename = paste(savedir,grp,'CNDRSpaceFit_basemodel+',goi,probe,'.pdf',sep=''),
       units = 'cm',height = 3.75,width = 3.75*length(tps))
lapply(1:length(tps),function(t) write.csv(data.frame(vuln=vulnerability[[t]]),paste0(savedir,grp,'vulnerability',tps[t],'MPI',goi,probe,'.csv')))
lapply(1:length(tps),function(t) write.csv(data.frame(pred=pred[[t]]),paste0(savedir,grp,'predictedpath',tps[t],'MPI',goi,probe,'.csv')))

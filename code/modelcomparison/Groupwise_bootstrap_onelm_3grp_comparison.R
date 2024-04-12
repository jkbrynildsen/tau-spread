#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','injection.site','group1','group2','group3')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'modelcomparison/bootstrap/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/optimfxns.R')
source('code/misc/miscfxns.R')
load(paste0(params$opdir,'processed/pathdata.RData'))  # load path data and ROI names

grps <- c(group1,group2, group3) # create list of the three groups being compared

###################################################
### load bootstrap data for bidirectional model ###
###################################################

mdl.name <- 'BidirectionalOneLM'

# load group 1
load(file = paste(savedir,group1,'ModelComparison_Bootstrap.RData',sep=''))
results.group1 <- results[[mdl.name]]$train
rm(results)

# load group 2
load(file = paste(savedir,group2,'ModelComparison_Bootstrap.RData',sep=''))
results.group2 <- results[[mdl.name]]$train
rm(results)

# load group 3
load(file = paste(savedir,group3,'ModelComparison_Bootstrap.RData',sep=''))
results.group3 <- results[[mdl.name]]$train 
rm(results)

MPI.names <- paste(6,'MPI') 
group.colors <- get3GroupColors(grps)

# 1. compare model fits between group1, group2 and group3
cfg <- list(fits.mse='MSE',fits.r='Pearson r') # loop through pearson r and MSE
p.fit <- list()
for(f.met in names(cfg)){ # compare fits for MSE and pearson r
  group1.r <- as.data.frame(do.call(rbind,lapply(results.group1, function(R) R[[f.met]])))
  group2.r <- as.data.frame(do.call(rbind,lapply(results.group2, function(R) R[[f.met]])))
  group3.r <- as.data.frame(do.call(rbind,lapply(results.group3, function(R) R[[f.met]])))
  df.all <- rbind(cbind(group1.r,Group=paste(group1),stringsAsFactors=F),cbind(group2.r,Group=paste(group2),stringsAsFactors=F),cbind(group3.r,Group=paste(group3),stringsAsFactors=F)) #exchange Treatment for Group
  colnames(df.all) <- c(MPI.names,'Group') #exchange Treatment for Group
  df.plt <- collapse.columns(df.all,cnames = MPI.names,groupby = 'Group') #exchange Treatment for Group
  df.plt$group <- factor(df.plt$group,levels = c(group1, group2, group3),ordered=T)
  #p.lab <- sapply(1:ncol(group1.r), function(tp) pval.2tail.np(0,group1.r[,tp]-group2.r[,tp]))
  p.fit[[cfg[[f.met]]]] <- ggplot(df.plt) + geom_boxplot(aes(x=names,y=values,fill=group),size=0.25,outlier.size=0.25) + theme_classic() +
    #annotate(geom='text',x=MPI.names,y=Inf,label=ifelse(p.lab<0.05,yes='*',no=''),],vjust=1,size=2)+
    ylab(cfg[[f.met]]) + xlab('') + 
    scale_fill_manual(limits=c(group1, group2, group3),values =group.colors,name='') + 
    scale_y_continuous(limits=c(0,NA))+
    theme(text=element_text(size=8),legend.key.size = unit(0.1,'cm'),legend.box.margin = ggplot2::margin(t = 0, unit='cm'),axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) 
  ggsave(p.fit[[cfg[[f.met]]]],filename = paste(savedir,group1,'vs',group2,'vs',group3,'Bootstrap',cfg[[f.met]],'.pdf',sep=''),
         units = 'cm',height = 6,width = 9,useDingbats=FALSE)
}

# 3. compare time constants across the three groups
group1.c <- as.data.frame(t(sapply(results.group1,function(R) c(Anterograde=R$c.train.antero,Retrograde=R$c.train.retro))))
group2.c <- as.data.frame(t(sapply(results.group2,function(R) c(Anterograde=R$c.train.antero,Retrograde=R$c.train.retro))))
group3.c <- as.data.frame(t(sapply(results.group3,function(R) c(Anterograde=R$c.train.antero,Retrograde=R$c.train.retro))))
df.all <- rbind(cbind(group1.c,Group=paste(group1),stringsAsFactors=F),cbind(group2.c,Group=paste(group2),stringsAsFactors=F),cbind(group3.c,Group=paste(group3),stringsAsFactors=F))
df.plt <- collapse.columns(df.all,cnames=c('Anterograde','Retrograde'),groupby = 'Group')
for(grp in grps){
  subset.list <- lapply(c(Anterograde='Anterograde',Retrograde='Retrograde'), function(n) df.plt$values[df.plt$names == n & df.plt$group==grp])
  subset.stats <- lapply(subset.list, function(X) c(mean=mean(X),quantile(X,c(0.025,.975)),median=median(X),cv=sd(X)/mean(X)))
  subset.stats$AnterogradeVsRetrogradeTimeConstant <- pval.np.pub(pval.2tail.np(0,subset.list$Retrograde-subset.list$Anterograde),length(subset.list$Retrograde))
  sink(paste0(savedir,grp,"TimeConstantStats.txt"))
  print(subset.stats)
  sink()
}

#p-value for group 1 vs group 2
p.group1.c.v.group2.c <- sapply(c('Anterograde','Retrograde'), function(d) pval.2tail.np(0,df.plt$values[df.plt$names==d & df.plt$group == group1]-df.plt$values[df.plt$names==d & df.plt$group == group2]))
p.group1.c.v.group2.c.lab <- data.frame(values=paste('p =',signif(p.group1.c.v.group2.c,2)))
p.group1.c.v.group2.c.lab$group <- c('Anterograde','Retrograde')
write.csv(x=p.group1.c.v.group2.c.lab,file=paste0(savedir,group1,'vs',group2,'AnteroAndRetroTimeConstant_Stats.csv'),row.names = F)

#p-value for group 1 vs group 3
p.group1.c.v.group3.c <- sapply(c('Anterograde','Retrograde'), function(d) pval.2tail.np(0,df.plt$values[df.plt$names==d & df.plt$group == group1]-df.plt$values[df.plt$names==d & df.plt$group == group3]))
p.group1.c.v.group3.c.lab <- data.frame(values=paste('p =',signif(p.group1.c.v.group3.c,2)))
p.group1.c.v.group3.c.lab$group <- c('Anterograde','Retrograde')
write.csv(x=p.group1.c.v.group3.c.lab,file=paste0(savedir,group1,'vs',group3,'AnteroAndRetroTimeConstant_Stats.csv'),row.names = F)

#p-value for group 2 vs group 3
p.group2.c.v.group3.c <- sapply(c('Anterograde','Retrograde'), function(d) pval.2tail.np(0,df.plt$values[df.plt$names==d & df.plt$group == group2]-df.plt$values[df.plt$names==d & df.plt$group == group3]))
p.group2.c.v.group3.c.lab <- data.frame(values=paste('p =',signif(p.group2.c.v.group3.c,2)))
p.group2.c.v.group3.c.lab$group <- c('Anterograde','Retrograde')
write.csv(x=p.group2.c.v.group3.c.lab,file=paste0(savedir,group2,'vs',group3,'AnteroAndRetroTimeConstant_Stats.csv'),row.names = F)

df.plt$group <- factor(df.plt$group,levels = grps,ordered=T)
p.c <- ggplot(df.plt) + geom_boxplot(aes(x=names,y=values,fill=group),size=0.25,outlier.size = 0.25) + theme_classic() +
  ylab('Diffusion Rate Constant') + xlab('') + scale_y_continuous(limits=NA)+scale_fill_manual(limits=grps,values=group.colors,name='')+ # used   ylim(0, 0.003)+ to print G20 0 vs 450 plot due to extreme outliers
  ylim(0, 0.01)+# include this only if outliers
  theme(text=element_text(size=8),legend.key.size = unit(0.1,'cm'),legend.box.margin = ggplot2::margin(t = 0, unit='cm'),axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) 
p.c
ggsave(p.c,filename = paste(savedir,group1,'vs',group2,'vs',group3,'TimeConstants.pdf',sep=''),
       units = 'cm',height = 6,width = 9,useDingbats=FALSE)

# 3. compare anterograde and retrograde betas between groups
if (group1 == "NTG_0 MLi-2"){cfg <- list("NTG_0 MLi-2"=results.group1, "NTG_75 MLi-2"=results.group2, "NTG_450 MLi-2"=results.group3)} else {cfg <- list("G20_0 MLi-2"=results.group1, "G20_75 MLi-2"=results.group2, "G20_450 MLi-2"=results.group3)}

for(grp in grps){
  # see code/misc/optimfxns.R, function lm.mask.ant.ret.all: x1 = retrograde, x2 = anterograde
  # ^ -- probably should change it so x1 is just called retro to begin with
  # extract coefficients weighting the importance of anterograde and retrograde betas
  ant.coefs <- as.data.frame(do.call(rbind,lapply(cfg[[grp]], function(R) R$m.train.coefs[[1]]['x2','Estimate'])))
  ret.coefs <-as.data.frame(do.call(rbind,lapply(cfg[[grp]], function(R) R$m.train.coefs[[1]]['x1','Estimate'])))
  
  df.plt <- rbind(cbind(ant.coefs,OI='Anterograde',stringsAsFactors=F),cbind(ret.coefs,OI='Retrograde',stringsAsFactors=F))
  colnames(df.plt) <- c('Beta','OI')
  df.plt <- collapse.columns(df.plt,cnames = 'Beta',groupby = 'OI')
  df.plt$Group <- grp
  df.plt <- df.plt[df.plt$values >= 0, ] # remove negative beta weights, which are uninterpretable
  p.avr <- pval.2tail.np(0,df.plt$values[df.plt$group=='Anterograde'] - df.plt$values[df.plt$group=='Retrograde'])
  p.avr <- p.avr
  p.lab.avr <- paste0('p = ',signif(p.avr,2))
  p.lab.avr[p.avr ==0] <- paste0('p < ',signif(1/length(cfg[[grp]]),2)) # don't say p = 0, say < 1/nboots
  cfg[[grp]] <- list(df=df.plt,p.avr=p.lab.avr,grp=grp) # store plot ready df in looped list
}

df.plt <- do.call(rbind,lapply(cfg,function(X) X$df)) # vertically concatenate df across groups then plot
# p-value comparing ant vs. ret for each group
p.avr.save <- do.call(rbind,lapply(cfg, function(X) c(X$p.avr,Group=X$grp,Pcrit=0.05/length(X$p.avr))))
write.csv(x=p.avr.save,file=paste0(savedir,group1,'_',group2,'_',group3,'AnterovsRetroBetas_Stats.csv'),row.names = F)
# compare group1 to group2
p.group1.v.group2 <- sapply(c('Anterograde','Retrograde'), function(d) pval.2tail.np(0, # loop through direction
                                                                                     cfg[[group1]]$df$values[cfg[[group1]]$df$group==d] - # compare group diff to 0
                                                                                       cfg[[group2]]$df$values[cfg[[group2]]$df$group==d]))
p.group1.v.group2.lab <- data.frame(values=paste('p =',signif(p.group1.v.group2,2)))
p.group1.v.group2.lab$group <- c('Anterograde','Retrograde')
# compare group1 to group3
p.group1.v.group3 <- sapply(c('Anterograde','Retrograde'), function(d) pval.2tail.np(0, # loop through direction
                                                                                     cfg[[group1]]$df$values[cfg[[group1]]$df$group==d] - # compare group diff to 0
                                                                                       cfg[[group3]]$df$values[cfg[[group3]]$df$group==d]))
p.group1.v.group3.lab <- data.frame(values=paste('p =',signif(p.group1.v.group3,2)))
p.group1.v.group3.lab$group <- c('Anterograde','Retrograde')
# compare group2 to group3
p.group2.v.group3 <- sapply(c('Anterograde','Retrograde'), function(d) pval.2tail.np(0, # loop through direction
                  cfg[[group2]]$df$values[cfg[[group2]]$df$group==d] - # compare group diff to 0
                    cfg[[group3]]$df$values[cfg[[group3]]$df$group==d]))
p.group2.v.group3.lab <- data.frame(values=paste('p =',signif(p.group2.v.group3,2)))
p.group2.v.group3.lab$group <- c('Anterograde','Retrograde')
# plot 3-group comparison of beta weights
df.plt$Group <- factor(df.plt$Group,levels = grps,ordered=T)
p <- ggplot() + geom_boxplot(data=df.plt,aes(x=group,y=values,fill=Group),size=0.25,outlier.size=0.25) + 
  theme_classic() +
  ylab('Standardized Beta') + xlab('') +scale_fill_manual(limits=grps,values =group.colors,name='')+
  ylim(0, 1.75)+
  theme(text=element_text(size=8),legend.key.size = unit(0.1,'cm'),legend.box.margin = ggplot2::margin(t = 0, unit='cm'),axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) 
p
ggsave(p,filename = paste(savedir,group1,'vs',group2,'vs',group3,'AnterogradeRetrogradeBetas_Boxplot.pdf',sep=''),
       units = 'cm',height = 6,width = 9,useDingbats=FALSE)
# save p-values from each pairwise comparison
write.csv(x=p.group1.v.group2.lab,file=paste0(savedir,group1,'vs',group2,'AnteroAndRetroBetas_Stats.csv'),row.names = F)
write.csv(x=p.group1.v.group3.lab,file=paste0(savedir,group1,'vs',group3,'AnteroAndRetroBetas_Stats.csv'),row.names = F)
write.csv(x=p.group2.v.group3.lab,file=paste0(savedir,group2,'vs',group3,'AnteroAndRetroBetas_Stats.csv'),row.names = F)

# save the mean and 95% CI of betas and time constants within each group
group1.c <- as.data.frame(t(sapply(results.group1,function(R) c(Anterograde=R$c.train.antero,Retrograde=R$c.train.retro))))
group2.c <- as.data.frame(t(sapply(results.group2,function(R) c(Anterograde=R$c.train.antero,Retrograde=R$c.train.retro))))
group3.c <- as.data.frame(t(sapply(results.group3,function(R) c(Anterograde=R$c.train.antero,Retrograde=R$c.train.retro))))
df.tc <- rbind(cbind(group1.c,Group=group1,stringsAsFactors=F),cbind(group2.c,Group=group2,stringsAsFactors=F),cbind(group3.c,Group=group3,stringsAsFactors=F))
df.tc <- collapse.columns(df.all,cnames=c('Anterograde','Retrograde'),groupby = 'Group')

beta.tc.stats <- data.frame()

for(a.r in c('Anterograde','Retrograde')){
  for(grp in grps){
    mask.beta <- df.plt$group==a.r &df.plt$Group==grp
    beta.tc.stats[paste(grp,a.r,'Beta'),'2.5%'] <- quantile(df.plt$values[mask.beta],0.025)
    beta.tc.stats[paste(grp,a.r,'Beta'),'Median'] <- median(df.plt$values[mask.beta])
    beta.tc.stats[paste(grp,a.r,'Beta'),'97.5%'] <- quantile(df.plt$values[mask.beta],0.975)
  }
}
for(a.r in c('Anterograde','Retrograde')){
  for(grp in grps){
    mask.tc <- df.tc$names == a.r & df.tc$group==grp
    beta.tc.stats[paste(grp,a.r,'Diffusion Time Constant'),'2.5%'] <- quantile(df.tc$values[mask.tc],0.025)
    beta.tc.stats[paste(grp,a.r,'Diffusion Time Constant'),'Median'] <- median(df.tc$values[mask.tc])
    beta.tc.stats[paste(grp,a.r,'Diffusion Time Constant'),'97.5%'] <- quantile(df.tc$values[mask.tc],0.975)
  }
}
write.csv(x=signif(beta.tc.stats,3),file=paste0(savedir,group1,'vs',group2,'vs',group3,'BetaTimeConstantDistributionTable.csv'))

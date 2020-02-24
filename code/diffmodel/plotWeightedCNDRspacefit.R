#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','goi','probe')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste0(params$opdir,'diffmodel/',goi,probe,'_Weighted/')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
load(paste(savedir,grp,goi,probe,'Weighted_CNDRSpaceFit_data.RData',sep=''))

# exclude regions with 0 pathology at each time point for purposes of computing fit
mask <- df$path != -Inf & !is.na(df$pred)
print(paste(sum(mask),'regions')) # number of regions left after exclusion
df <- df[mask,] 
# use linear model to predict path from connectivity and synuclein
c.test <- cor.test(df$pred,df$path)
print(paste(nrow(df),'regions'))

p <- ggplot(df,aes(x=pred,y=path)) + geom_smooth(color = '#007257',method ='lm',size=1) + geom_point(color = '#007257',size = 1,alpha=0.6,stroke=0) +
  annotate(geom='text',x=max(df$pred),y=min(df$path) + 0.1,label = paste('r =',signif(c.test$estimate,2)),size=2.5,hjust=1,vjust=0) +
  #scale_x_continuous(limits=c(min(df$Xt),max(df$Xt))) +
  theme_classic() + xlab('Predicted') + ylab('log(Path.)') + ggtitle(grp) +
  theme(text = element_text(size=8),plot.title = element_text(hjust=0.5,size=8))

ggsave(p,filename = paste(savedir,grp,goi,probe,'Weighted_CNDRSpaceFit_predictedpathcontinuous.pdf',sep=''),units = 'in',height = 1.5,width = 1.5)

###############################################
### Save predicted values and vulnerability ###
###############################################

m <- lm(path~pred,data=df)
vulnerability <- residuals(m)
pred <- m$fitted.values

write.csv(data.frame(vuln=vulnerability),paste(savedir,grp,goi,probe,'Weighted_vulnerability.csv',sep=''))
write.csv(data.frame(pred=pred),paste(savedir,grp,goi,probe,'Weighted_predictedpath.csv',sep=''))

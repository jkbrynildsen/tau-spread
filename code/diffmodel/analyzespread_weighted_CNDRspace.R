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
gene.exp <- read.csv('data/aba/expression/GeneExpressionABA.csv') # load gene expression
gene.exp <- as.matrix(gene.exp[,paste0(goi,'_',probe)])
rownames(gene.exp) <- region.names
missing.mask <- !is.na(gene.exp) # eliminate regions that had no expression data or no correspondence 

# get mean pathology for only time point, 1 month
tp <- 1
Mice <- path.data[path.data$Condition == grp,-1]
Grp.mean <- colMeans(Mice,na.rm = T)

W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
L.out <- get.Lout(W[missing.mask,missing.mask],gene.exp[missing.mask]) # compute out-degree Laplacian weighted by gene.exp, removing regions with no gene expression data

# Fit time scaling parameter on average of all mice
#l1 <- Re(eigen(x = L.out,symmetric = FALSE,only.values = TRUE)$values) # choose maximum c value based on smallest positive eigenvalue of L.out
#c.max <- log(100)/(l1[length(l1)-1]*max(tp))
c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n) # scaling parameter
log.path <- log(Grp.mean,base=10)
Xo <- get.Xo(region.names[missing.mask],'iCP') # seed pathology in iCPu
c.Grp <- c.CNDRspace.fit(log.path,L.out,Xo,c.rng,ABA.to.CNDR.key)

# for the remainder of the script, only evaluating path at 1 month
# so set log path equal to only 1 month

Xt.Grp <- log(map.ABA.to.CNDR(predict.Lout(L.out,Xo,c.Grp),path.names,ABA.to.CNDR.key), base = 10) # predict pathology using connectivity, time constant, and seed
#Xt.Grp <- scale(Xt.Grp,center=T)
df <- data.frame(path = log.path, pred = Xt.Grp)
save(df,c.Grp,file = paste(savedir,grp,goi,probe,'Weighted_CNDRSpaceFit_data.RData',sep=''))

#################
### Load data ###
#################

rm(list=setdiff(ls(),'params'))
basedir <- params$basedir
setwd(basedir)
savedir <- paste0(basedir,'data/')

data <- read.csv(paste0(savedir,'TauPathDataRaw.csv'),stringsAsFactors = F)
data <- t(data)
colnames(data) <- data['Region',]
data <- as.data.frame(data[-which(rownames(data)=='Region'),])
data$Condition <- substr(rownames(data),start=1,stop=3) # extract group label from row names
data$Month <- as.numeric(substr(rownames(data),start=5,stop=5)) # extract MPI
# reorder columns
region.names <- setdiff(colnames(data),c('Condition','Month'))
data <- data[,c('Condition','Month',region.names)]

write.csv(x=data,file = paste0(savedir,'PathData.csv'),row.names = F)

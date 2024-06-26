### Prepare structural connectome matrix for network analyses

rm(list=setdiff(ls(),'params'))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'processed/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')

# read in raw pathology dataset to get region names
path.data <- read.csv(paste0(basedir,'data/TauMLi2_Parent_RemovedFT.csv'), header = TRUE)
path_regions <- path.data$X # create list of region names from pathology dataset

# read in structural connectivity dataset
struct_mat <- read.csv(paste0(basedir,"data/normalized_connection_strength.csv"))
struct_mat2 <- struct_mat[,-1] # eliminate first row, which contains regions names
row.names(struct_mat2) <- struct_mat[,1] # convert first column into row names

# subset connectivity matrix into ipsi and contralateral hemispheres
connectivity.contra <- struct_mat2 %>% 
  select(starts_with("contra"))
colnames(connectivity.contra) <- connectivity.contra[c(1),] # make column names match first row
connectivity.contra <- connectivity.contra[-c(1),] # remove first row
connectivity.contra <- connectivity.contra[, names(connectivity.contra) %in% path_regions] # select only columns that are present in the pathology dataset
connectivity.contra <- connectivity.contra[rownames(connectivity.contra) %in% names(connectivity.contra),] # select only rows that are present in the pathology dataset

connectivity.ipsi <- struct_mat2 %>% 
  select(starts_with("ipsi"))
colnames(connectivity.ipsi) <- connectivity.ipsi[c(1),] # make column names match first row
connectivity.ipsi <- connectivity.ipsi[-c(1),] # remove first row
connectivity.ipsi <- connectivity.ipsi[, names(connectivity.ipsi) %in% names(connectivity.contra)] # select only columns that are present in the pathology dataset
connectivity.ipsi <- connectivity.ipsi[rownames(connectivity.ipsi) %in% names(connectivity.contra),] # select only rows that are present in the pathology dataset
# NOTE: the above excludes 6 regions that are missing contralateral structural edges

conn.names.ipsi <- colnames(connectivity.ipsi)
conn.names.contra <- colnames(connectivity.contra)

# checks 
if(identical(colnames(connectivity.contra),rownames(connectivity.contra))){
  print('contra connectivity colnames and rownames equal')}
if(identical(colnames(connectivity.ipsi),rownames(connectivity.ipsi))){
  print('ipsi connectivity colnames and rownames equal')}

# create elements of full connectivity matrix
upper.left <- connectivity.ipsi
upper.right <- connectivity.contra
lower.left <- connectivity.contra
lower.right <- connectivity.ipsi

W <- rbind(cbind(upper.left,upper.right),cbind(lower.left,lower.right))

rownames(W) <- c(paste('i',rownames(upper.left),sep=''), # add i to ipsilateral regions
                 paste('c',rownames(lower.left),sep='')) # add c to contralateral regions
colnames(W) <- c(paste('i',rownames(upper.left),sep=''), # add i to ipsilateral regions
                 paste('c',rownames(lower.left),sep='')) # add c to contralateral regions

n.regions.ABA <- nrow(W)
n.regions.ABA.hemi <- n.regions.ABA/2

# check if connectivity was tiled into alternating blocks correctly

unit.test(all(W[(n.regions.ABA.hemi+1):n.regions.ABA,(n.regions.ABA.hemi+1):n.regions.ABA] == W[1:n.regions.ABA.hemi,1:n.regions.ABA.hemi]),
          'tiling on-diagonal blocks worked','tiling on-diagonal blocks failed') 
unit.test(all(W[1:n.regions.ABA.hemi,(n.regions.ABA.hemi+1):n.regions.ABA] == W[(n.regions.ABA.hemi+1):n.regions.ABA,1:n.regions.ABA.hemi]),
          'tiling off-diagonal blocks worked','tiling off-diagonal blocks failed')

unit.test(all(colnames(W) == rownames(W)),'row and column names of conn mat are same','ERROR with conn mat names')

PathData <- read.csv(paste(params$opdir,'processed/PathData.csv', sep=''))
filter.path.names <- colnames(PathData)[-c(1:5)] # NOTE that this includes more regions than the connectivity dataset

W <- W[, names(W) %in% filter.path.names] # select only columns that are present in the pathology dataset
W <- W[rownames(W) %in% filter.path.names,] # select only rows that are present in the pathology dataset

W_labeled <- W # store a version of the structural matrix with region labels

n.regions.ABA <- nrow(W)
n.regions.ABA.hemi <- n.regions.ABA/2

region.names <- colnames(W) # store connectivity matrix names

# read in pathology data, match regions with those in the structural dataset
path.data <- PathData[, names(PathData) %in% region.names]
path.data <- cbind(PathData[,1:5], path.data)
path.names <- colnames(path.data)[-c(1:5)]

save(path.data, path.names,region.names, n.regions.ABA, n.regions.ABA.hemi, file = paste(savedir,'pathdata.RData',sep=''))

write.table(W, paste0(savedir,'W.csv'), sep=',', row.names = FALSE, col.names = FALSE)
write.csv(W_labeled, paste0(savedir,'W_labeled.csv'))

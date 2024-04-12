### Create matrix of interregional distance data for null models

rm(list=setdiff(ls(),c('params')))
dist <- read.csv("Ohetal_2014_distance_mat.csv")

dist2 <- dist[,-1]
row.names(dist2) <- dist[,1]

# read in raw pathology dataset to get region names

struct_lab <- read.csv(paste(params$opdir,'processed/W_labeled.csv', sep=''))
region_labels <- struct_lab$X # create list of region names from pathology dataset
# subset contralateral connectivity matrix
dist.contra <- dist2 %>% 
  select(ends_with("_contra")) %>% # pull out contra columns
  filter(grepl("_contra", rownames(dist2))) # pull out contra rows
# rename regions to match other datasets
dist.contra <- dist.contra %>% 
  rename_with(~str_remove(., '_contra'))
names(dist.contra) <- gsub("\\.","",names(dist.contra)) # remove period from names
names(dist.contra) <- c(paste('c',names(dist.contra),sep='')) # add c to contralateral regions

rownames(dist.contra) <- names(dist.contra) # make row names equivalent to column names

# subset ipsilateral connectivity matrix
dist.ipsi <- dist2 %>% 
  select(ends_with("_ipsi")) %>% # pull out ipsi columns
  filter(grepl("_ipsi", rownames(dist2))) # pull out ipsi rows
# rename regions to match other datasets
dist.ipsi <- dist.ipsi %>% 
  rename_with(~str_remove(., '_ipsi'))
names(dist.ipsi) <- gsub("\\.","",names(dist.ipsi)) # remove period from names
names(dist.ipsi) <- c(paste('i',names(dist.ipsi),sep='')) # add i to ipsilateral regions

rownames(dist.ipsi) <- names(dist.ipsi) # make row names equivalent to column names

all_dist <- rbind(cbind(dist.ipsi,dist.contra),cbind(dist.contra,dist.ipsi)) # tile ipsi and contralateral distances into a matrix

all_dist_names <- colnames(all_dist) # save list of all region names in distance

# select rows in the pathology dataset that are present in the distance dataset
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
#added the below chunk
PathData <- read.csv(paste(params$opdir,'processed/PathData.csv', sep=''))
path.data <- PathData[, names(PathData) %in% region.names]
path.data <- cbind(PathData[,1:5], path.data)
path.names <- colnames(path.data)[-c(1:5)]

path_dist_names <- intersect(colnames(all_dist), colnames(path.data)) # get list of names present in both pathology and distance datasets

path_dist <- path.data %>%
  select(all_of(path_dist_names))

# select rows in the distance dataset that are present in the pathology dataset
dist_path <- all_dist[, names(all_dist) %in% path_dist_names] # select only columns that are present in the pathology dataset
dist_path <- dist_path[rownames(dist_path) %in% path_dist_names,] # select only rows that are present in the pathology dataset

#save(dist_path, path.names,region.names, n.regions.ABA, n.regions.ABA.hemi, file = paste(savedir,'pathdata.RData',sep=''))

write.table(dist_path, paste(params$opdir,'processed/Distance.csv', sep=''), sep=",", col.names = FALSE, row.names = FALSE)
write.csv(dist_path, paste(params$opdir,'processed/Distance_labeled.csv', sep=''))

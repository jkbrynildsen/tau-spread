### Prepare pathology datasets for analyses

savedir <- paste0(params$opdir, 'processed/')
dir.create(savedir,recursive = T)

# load and format pathology dataset
parent_tau <- read.csv(paste0(basedir,"data/TauMLi2_Parent_RemovedFT.csv"), row.names = 1,) # sublaminar pathology (360 regions per hemisphere)
parent_tau <- parent_tau %>%
  select(-contains(c("TM1H", "TM1J", "TM4D", "TM4H", "TM5D", "TM6B", "TM6F", "TM1B", "TM3L"))) # remove mice with incorrect genotypes or offsite injections
parent_tau <- as.data.frame(t(parent_tau)) # transpose data frame so mice are rows
parent_tau <- rownames_to_column(parent_tau, "Mouse.ID")

# split row names so hemisphere is its own column
rn <- parent_tau$Mouse.ID
m1 <- do.call(rbind, strsplit(rn, '_')) # split mouse ID labels from hemisphere labels
parent_tau['Mouse.ID'] <- m1[,1] 
parent_tau['Hemisphere'] <- m1[,2] # add column for hemisphere label

# subset pathology data into two data frames by hemisphere
Rparent_tau <- parent_tau %>% filter(Hemisphere == "Right")
Lparent_tau <- parent_tau %>% filter(Hemisphere == "Left")

# load mouse group info
mouse_info <- read.delim(paste0(params$basedir,"data/TauMLi2 Project Animal Key.txt"))
str_sub(mouse_info$Mouse.ID, 4, 4) <- "" # remove dash in ID name to match with pathology data
 
# align mouse group info with pathology data, reorder columns so mouse info comes before data columns
RPathData <- inner_join(Rparent_tau, mouse_info, by = "Mouse.ID")
RPathData <- RPathData[,c(1,362:366,2:361)] # reorder columns so mouse info is listed first
colnames(RPathData)[7:366] <- paste("i", colnames(RPathData[7:366]), sep = "") # add "i" prefix to right hemisphere regions

LPathData <- inner_join(Lparent_tau, mouse_info, by = "Mouse.ID")
LPathData <- LPathData[,c(1,362:366,2:361)]
colnames(LPathData)[7:366] <- paste("c", colnames(LPathData[7:366]), sep = "") # add "c" prefix to left hemisphere regions

PathData <- cbind(RPathData, LPathData[,7:366]) # join data from each hemisphere
names <- c('Hemisphere', 'Genotype', 'Sex', 'Treatment')
PathData[,names] <- lapply(PathData[,names], factor)
levels(PathData$Genotype) <- c("NTG", "G20")
PathData_Final <- subset(PathData, select=-c(2)) # remove Hemisphere column

# remove columns where no pathology was measured (all NA and/or 0)
PathData_Final <- PathData_Final %>% 
  select_if(~!all(is.na(.) | . == 0))

# write out processed data to a CSV file
write.csv(PathData_Final, paste0(savedir, 'PathData.csv'), row.names = F)

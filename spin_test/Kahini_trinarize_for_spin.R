#Title: Trinarize files for spin

### pre: logK text files, 10242 vertices
### post: 6 matrices with 1s if value fdr corrected, 0 if not, and -1 if vertex in medial wall of yeo map
# lh and rh pos and neg together, positive alone, and negative alone
### uses: In order to spin my vertices in a way that allows me to later assess their yeo membership, I need all of my T values indicated, as well as which values are in the medial wall.  This script takes lh and rh fdr corrected vectors generated in previous scripts, booleanizes them into sig/non sig based on T score. Then, -1s are added for vertices in the medial wall. 
### dependencies: and R will do. I used 3.2.5

## set abs and relative paths. Must toggle before running locally/on cluster
homedir = "/Users/eballer/BBL/kahini_spin/"


source(paste0(homedir, '/imco_functions.R'))


### get_parcel_mapping_yeo_function
get_parcel_mapping_yeo <- function(parcel_num){
  
  #pre: input parcel #, either 7 or 17
  #post: list lh and rh yeo networks that map onto code #s
  #uses: easy way to translate the weird numerical maps in fsaverage 5 space into something we are more familiar with
  #dependencies: Any R will do, I used 3.2.5
  
  ## Set Yeo info
  #### set # parcels in case I want to do 7 or 17 or something else in the future
  parcel_type = "Yeo" 
  parcel_num = parcel_num 
  input_parcel_array_length = 10242
  
  # read in yeo fsaverage5 vectors
  parcelID <- read.csv(paste0(homedir, "NetworkIDnumbers", parcel_type, parcel_num, ".csv"), header = F)
  parcelName <- t(read.csv(paste0(homedir, "NetworkNames", parcel_type, parcel_num, ".csv"), header = F))
  
  lh_parcel_nums <- read.csv(paste0(homedir, "/lh_", input_parcel_array_length, "_vertex_nums_", parcel_type, parcel_num, ".csv"), header = F)
  rh_parcel_nums <- read.csv(paste0(homedir, "/rh_", input_parcel_array_length, "_vertex_nums_", parcel_type, parcel_num, ".csv"), header = F)
  
  # map Yeo numbers to parcels
  #make a column of numbers for mapping
  parcelID$network_num <- c(1:dim(parcelID)[1])
  
  #add extra row to parcelID, not clear why this didn't come from Yeo labels, maybe cerebellum?... 8 will equal 65793
  # comment this out if not using yeo 
  parcelID<- rbind(parcelID, c(65793, 8))
  
  #make vector for lh and rh with mapping
  lh_numerical_map <- lh_parcel_nums
  rh_numerical_map <- rh_parcel_nums
  
  #foreach vertex, which contains a bunch of numbers, match it to the appropriate column, and take the network num (i.e. yeo 2, which would correspond to Motor), associated with it
  lh_numerical_map[] <- lapply(lh_parcel_nums, function(x) parcelID$network_num[match(x, parcelID$V1)])
  rh_numerical_map[] <- lapply(rh_parcel_nums, function(x) parcelID$network_num[match(x, parcelID$V1)])
  
  lh_and_rh_numerical_map_list <- list(lh_numerical_map$V1, rh_numerical_map$V1)
  return(lh_and_rh_numerical_map_list)
  
}


#### set # parcels in case I want to do 7 or 17 or something else in the future
parcel_type = "Yeo" 
parcel_num = 7 
input_parcel_array_length = 10242

#set threshold
thresh = 3.09

#will loop through each of these
analyses <- c("logk") 
corrs <- c(thresh) #correction

parcel_mapping <- get_parcel_mapping_yeo(7)
lh_numerical_map <- parcel_mapping[[1]]
rh_numerical_map <- parcel_mapping[[2]]

for (analysis in analyses) {
    for (corr in corrs) {  
      
      ### set results path
    #  stat_path <- paste0("/", analysis, "/")
     # print(stat_path)
      result_path <- paste0(analysis, "_t_", corr)
      print(result_path)
      
      ### set paths
      ## input
      lh_stat_map <- read.csv(paste0(homedir, "L_logk_data.txt"), header = T)
      rh_stat_map <- read.csv(paste0(homedir, "R_logk_data.txt"), header = T)
      
      ## output
      lh_outdir <- paste0(homedir, "spin_results/lh_", result_path, "_", parcel_type, parcel_num, "_1_0_-1.csv")
      rh_outdir <- paste0(homedir, "spin_results/rh_", result_path, "_", parcel_type, parcel_num, "_1_0_-1.csv")
      
      lh_outdir_pos <- paste0(homedir, "spin_results/lh_pos_", result_path, "_", parcel_type, parcel_num, "_1_0_-1.csv")
      rh_outdir_pos <- paste0(homedir, "spin_results/rh_pos_", result_path, "_", parcel_type, parcel_num, "_1_0_-1.csv")
      
      lh_outdir_neg <- paste0(homedir, "spin_results/lh_neg_", result_path, "_", parcel_type, parcel_num, "_1_0_-1.csv")
      rh_outdir_neg <- paste0(homedir, "spin_results/rh_neg_", result_path, "_", parcel_type, parcel_num, "_1_0_-1.csv")
      
      #convert stat map to boolean
      lh_stat_boolean_pos <- ifelse(lh_stat_map >= thresh, 1, 0)
      rh_stat_boolean_pos <- ifelse(rh_stat_map >= thresh, 1, 0)
      
      lh_stat_boolean_neg <- ifelse(lh_stat_map <= -thresh, 1, 0)
      rh_stat_boolean_neg <- ifelse(rh_stat_map <= -thresh, 1, 0)
      
      lh_stat_boolean <- lh_stat_boolean_pos + lh_stat_boolean_neg
      rh_stat_boolean <- rh_stat_boolean_pos + rh_stat_boolean_neg
      
      
      #convert NAs to 0
      lh_stat_boolean[is.na(lh_stat_boolean)] <- 0
      rh_stat_boolean[is.na(rh_stat_boolean)] <- 0
      
      lh_stat_boolean_pos[is.na(lh_stat_boolean_pos)] <- 0
      rh_stat_boolean_pos[is.na(rh_stat_boolean_pos)] <- 0
      
      lh_stat_boolean_neg[is.na(lh_stat_boolean_neg)] <- 0
      rh_stat_boolean_neg[is.na(rh_stat_boolean_neg)] <- 0
      
      
      #substitute -1 for those locations with medial wall stuff
      lh_medial_wall_nums <- which(lh_numerical_map == 8)
      rh_medial_wall_nums <- which(rh_numerical_map == 8)
      
      lh_stat_boolean[lh_medial_wall_nums] <- -1
      rh_stat_boolean[rh_medial_wall_nums] <- -1
      
      lh_stat_boolean_pos[lh_medial_wall_nums] <- -1
      rh_stat_boolean_pos[rh_medial_wall_nums] <- -1
      
      lh_stat_boolean_neg[lh_medial_wall_nums] <- -1
      rh_stat_boolean_neg[rh_medial_wall_nums] <- -1
    
      #write output
      write.table(x = lh_stat_boolean, file = lh_outdir, quote = F, row.names = F, col.names = F)
      write.table(x = rh_stat_boolean, file = rh_outdir, quote = F, row.names = F, col.names = F)
      
      write.table(x = lh_stat_boolean_pos, file = lh_outdir_pos, quote = F, row.names = F, col.names = F)
      write.table(x = rh_stat_boolean_pos, file = rh_outdir_pos, quote = F, row.names = F, col.names = F)
      
      write.table(x = lh_stat_boolean_neg, file = lh_outdir_neg, quote = F, row.names = F, col.names = F)
      write.table(x = rh_stat_boolean_neg, file = rh_outdir_neg, quote = F, row.names = F, col.names = F)
    }
}



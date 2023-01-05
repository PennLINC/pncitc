##############################################
########### Spin Test Distribution for #######
##########Kahinir##########

####### Author: Erica Baller
#### Date: 12/20/22

#######
##pre: right and left 10242 x 1000 matrices from matlab SpinPermuFS, yeo R & L assignments
##post: 2 7 x 1000 matrices (r & l) that contain the proportion of vertices within a network divided by the total number of vertices, and plots
## uses: Takes output of spin test, and calcualted the number of vertices within each of yeo's 7 networks out of the number of total possible vertices within the network
    #### 1) Read in the yeo network assignments and calculate total number of vertices per network
    #### 2) Multiply the yeo networks x the matrices (so every value is 1 -7 if they were within the mask, -1--7 if they were medial wall, and 0 otherwise)
    #### 3) Foreach permutation (r and l separately), and for each network, calculate the (# of vertices with a 1) divided(/) by the (number of total vertices within network minus number of negative vertices
    #### 4) Store
    #### 5) Plot

### dependencies: ggplot2, bigmemory, vroom


#library(bigmemory.sri)
library(ggplot2)
library(tidyr)



#################
### set home directory
#homedir <- "/Users/eballer/BBL/imco/pmacs/PMACS_remote"
homedir <- "/Users/eballer/BBL/kahini_spin/"

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



#initialize
hemis <- c("lh", "rh")
permNum <- 1000
yeo_num <- 7
thresh <- 3.09
models = c("logk", "pos_logk", "neg_logk")

#set flag to 1 if you'd like to calculate a spin for a mean map. It uses the means rather than proportions so it is a little wee bit different
mean_coupling = 0
################
### Read in matrices 
for (model in models) {
  #actual results
  lh_results <- read.table(paste0(homedir, "/spin_results/lh_", model, "_t_", thresh, "_Yeo7_1_0_-1.csv"))
  rh_results <- read.table(paste0(homedir, "/spin_results/rh_", model, "_t_", thresh, "_Yeo7_1_0_-1.csv"))
  
  #spins
  
  #extract suffix
  lh_spin <- t(read.table(paste0(homedir, "/spin_results/lh_spin_test_spin_results_", model,"_", thresh, "_output.csv"), sep = ","))
  rh_spin <- t(read.table(paste0(homedir, "/spin_results/rh_spin_test_spin_results_", model,"_", thresh, "_output.csv"), sep = ","))
                                          
  #bring together, with original values as first column
  lh_act_results_and_spin <- cbind(lh_results, lh_spin)
  rh_act_results_and_spin <- cbind(rh_results, rh_spin)
  
  #grab list of yeo 7 networks in fsaverage5 space
  yeo_networks <- get_parcel_mapping_yeo(yeo_num)
  
  #separate into right and left
  lh_yeo_network <- yeo_networks[[1]]
  rh_yeo_network <- yeo_networks[[2]]
  
  #count up number of vertices per network
  lh_yeo_network_count_table <- table(lh_yeo_network)
  rh_yeo_network_count_table <- table(rh_yeo_network)
  
  #multiply yeo network x spin test
  lh_spinxyeo <- lh_act_results_and_spin*lh_yeo_network
  rh_spinxyeo <- rh_act_results_and_spin*rh_yeo_network
  
  #proportions
  #go through each hemisphere, go through each perm, and go through each network
  
  lh_hemi_spin_proportions <- data.frame(matrix(nrow = yeo_num, ncol = (permNum + 1)))
  rh_hemi_spin_proportions <- data.frame(matrix(nrow = yeo_num, ncol = (permNum + 1)))
  for (hemi in hemis){

    for (perm in 1:(permNum + 1)){
      
      for (network in 1:yeo_num){
        
        #to evaluate
        
        #number of vertices within network that are fdr corrected
        num_pos_to_parse<- paste0("length(which(", hemi, "_spinxyeo[", perm, "] == ", network, "))")
        
        num_vertices_in_spin <- eval(parse(text = as.character(num_pos_to_parse)))
        
        
        #number of vertices within network that are negative (i.e., medial wall)
        num_neg_to_parse <- paste0("length(which(", hemi, "_spinxyeo[", perm, "] == -", network, "))")
      
        num_neg <- eval(parse(text = as.character(num_neg_to_parse)))
        
        
        #total number of vertices in normal network
        total_possible_to_parse <- paste0(hemi, "_yeo_network_count_table[", network, "]")
 
        total_possible <- eval(parse(text = as.character(total_possible_to_parse)))
        
       
        #proportion of vertices within network , with denominator being total possible by # in medial wall
        proportion_potential_vertices <- num_vertices_in_spin/(total_possible - num_neg)
    
        
        #store in matrix
        storing_to_parse <- paste0(hemi, "_hemi_spin_proportions[", network, ",", perm, "] = ", proportion_potential_vertices)

        eval(parse(text = as.character(storing_to_parse)))
      }
    }
  }
  
  write.table(lh_hemi_spin_proportions, file = paste0(homedir, "/spin_results/lh_spin_test_", model, "_", thresh, "_proportions.csv"), sep = ",", col.names = F, row.names = F)
  write.table(rh_hemi_spin_proportions, file = paste0(homedir, "/spin_results/rh_spin_test_", model, "_", thresh, "_proportions.csv"), sep = ",", col.names = F, row.names = F)
#then plot
}


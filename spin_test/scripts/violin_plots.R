####### Violin Plots ###


#### Author: Erica Baller

### 3/16/2021

###pre: requires that you have run the spin tests for your models, including the positive and negative spin tests
###post: 2 violin plots per model. Violin plot x axis - yeo networks, y axis - proportion of vertices within networks. dotted lines - mean from permutations
   #### plot 1 will contain a thick black line for the *actual values from your analysis so you can compare how far above or below your value is from permuted
   #### plot 2 will contain red lines detailing where the actual values from your analysis from the POSITIVE domain, blue for NEGATIVE
   #### You will also get a table that tells you the p value (uncorrected) for each network within each model.
   #### for my final plot, I added a * to these manually
### uses: Creates violin plots and analyses to help make sense of the results from permutation spin tests. The question we are asking is:
    ### what is the likelihood that the number of vertices within a network is significant rather than due to chance
### dependencies: Using R 3.6.3 but any R will do. Libraries to include listed below

library(tidyr)
library(ggplot2)
library(reshape)

### set home directory
homedir <- "/Users/eballer/BBL/kahini_spin/"
thresh <- 3.09

violin_plot <- function (homedir, models, network_names, num_spins){
  for (model in models) {
    #for storing statistics at the end
    lh_spin_df <- data.frame(read.table(paste0(homedir, "/spin_results/lh_spin_test_", model, "_", thresh, "_proportions.csv"), sep = ",")
    )
    rh_spin_df <- data.frame(read.table(paste0(homedir, "/spin_results/rh_spin_test_", model, "_", thresh, "_proportions.csv"), sep = ",")
    )
    
    
    #take means of left and right
    actual_results <- (lh_spin_df[,1] + rh_spin_df[,1])/2
    
    #dataframes for all spins
    spin_without_target_col <- cbind(lh_spin_df[,2:1001],rh_spin_df[,2:1001])
    melted_df <- melt_df_for_violin_plot_yeo7(spin_without_target_col, network_names, num_spins)

    plot_violin <- ggplot(melted_df, aes(x = factor(network_name, level = network_names), y = spin, fill = network_name)) +  
      scale_fill_manual(values=yeo_colors) + 
      geom_violin(trim = TRUE) + 
      xlab("Yeo 7 Network") + ylab(paste0("Proportion")) +
      ylim(0,NA) + 
      theme_classic() + 
      theme(legend.position = "none",
            legend.title = element_blank(),
            axis.text.x = element_text(size = 10, colour = "black"),
            axis.text.y = element_text(size = 10, colour = "black"),
            axis.title.y = element_text(size = 10),
            axis.title.x = element_blank(),
            plot.title = element_text(size = 10)) +
      geom_segment(aes(x = 0.6, y = actual_results[1], xend = 1.4, yend = actual_results[1])) + 
      geom_segment(aes(x = 1.6, y = actual_results[2], xend = 2.4, yend = actual_results[2])) +
      geom_segment(aes(x = 2.6, y = actual_results[3], xend = 3.4, yend = actual_results[3])) +
      geom_segment(aes(x = 3.6, y = actual_results[4], xend = 4.4, yend = actual_results[4])) +
      geom_segment(aes(x = 4.6, y = actual_results[5], xend = 5.4, yend = actual_results[5])) +
      geom_segment(aes(x = 5.6, y = actual_results[6], xend = 6.4, yend = actual_results[6])) +
      geom_segment(aes(x = 6.6, y = actual_results[7], xend = 7.4, yend = actual_results[7]))
    ggsave(plot=plot_violin, filename = paste0(homedir, "/spin_results/violin_", model, "_",thresh, ".pdf"), width = 4.81, height = 4.81)
    plot_violin
  }
}


#################


source(paste0(homedir, '/scripts/imco_functions_kahini.R'))

outdir_name <- "spin_stats"

models_for_stats = c("logk", "pos_logk", "neg_logk")
models_for_plots = c("logk", "pos_logk", "neg_logk")

network_names <- c("VIS", "MOT", "DA", "VA", "LIM", "FP", "DM")
num_spins = 2000
#get Yeo colors from function, these values were set manually - this can be obtained through https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation_Yeo2011
   #I typed the RGB values into a rgb -> hex converter and stored the values here. Works!

yeo_colors <- get_yeo7_colors()

#for storing statistics at the end
stats <- data.frame(matrix(nrow = 7, ncol = length(models_for_stats)))
names(stats) <- models_for_stats
row.names(stats) <- network_names

for (model in models_for_stats) {
  

    lh_spin_df <- data.frame(read.table(paste0(homedir, "/spin_results/lh_spin_test_", model, "_", thresh, "_proportions.csv"), sep = ","))
    rh_spin_df <- data.frame(read.table(paste0(homedir, "/spin_results/rh_spin_test_", model, "_", thresh, "_proportions.csv"), sep = ","))
  
  
  
  #take means of left and right
  actual_results <- (lh_spin_df[,1] + rh_spin_df[,1])/2
  
  #dataframes for all spins
  spin_without_target_col <- cbind(lh_spin_df[,2:1001],rh_spin_df[,2:1001])
  
  ####### Stats ########
  for (i in 1:7){
    #store p values
    #equivalen to: stats$model[i] <- (length(which(spin_without_target_col[i,] > actual_results[i]))/2000)
    eval(parse(text=as.character(paste0("stats$", model, "[", i ,"] <- (length(which(spin_without_target_col[",i, ",] > actual_results[", i, "]))/", num_spins, ")"))))
  }
}

print(stats)
violin_plot(homedir = homedir, models = models_for_plots, network_names = network_names, num_spins = num_spins)
write.csv(stats, paste0(homedir, "/spin_results/", outdir_name, ".csv"))


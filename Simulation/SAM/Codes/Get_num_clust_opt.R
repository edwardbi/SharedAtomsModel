rm(list=ls())

##################################################################
# Command line arguments
##################################################################

args = commandArgs(TRUE)

if(length(args) != 5){
  stop("Usage: CAM_script.R <seed_n> <work_dir> <input_dat> <sim_data_num> <out_dir>")
}

work_dir = args[1]
opt_clust = args[2]
sim_data_num = as.numeric(args[4])
out_dir = args[5]

setwd(work_dir)

##################################################################
# Pcakages and Sources
##################################################################
CAM_file1 <- "SAM1_univ_optcluster_"
seed_num <- c(12345, 1, 4, 1, 2, 7, 4, 3, 12345, 3, 
              8, 5, 9, 6, 8, 2, 1, 5, 12345, 4,
              4, 10, 3, 12345, 10, 8, 4, 2, 1, 8)

for(i in 1:30){
    opt_clust <- paste0(CAM_file1,seed_num[i],"_",i,".RData")
    optmal_cluster <- readRDS(opt_clust)
    clus <- optmal_cluster$cl[1,]
    num_c <- length(unique(clus))
    print(c(i, num_c))
}


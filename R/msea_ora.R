rm(list=ls(all=TRUE))

#---------------------------------
#  Prepare input data for analysis
#---------------------------------

# data files
path <- "./query data"

L <- list.files(path)  # List of data files
data <- paste(path, L[1], sep = "/")  # demo

# MSEA list
#iDMET_pathwaylist_ori <- "./data/iDMET_enrimentlist.csv"
iDMET_pathwaylist_ori <- "iDMET_pathwaylist_set01.csv"

# output files
out_file <- "./out put/MSEA result/MseaResult.csv"
out_file_sub <- "./out put/MSEA result/unmatch_metabo.csv"

# Load source file
source(file = "./iDMET_msea_source.R")





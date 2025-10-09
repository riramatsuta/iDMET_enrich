# rm(list=ls(all=TRUE))

#-----------
# Load files
#------------
# Load analysis results
file <- "iDMET_result_commonmetabolitelist.csv"
x <- read.csv(file, fileEncoding = "shift-jis", header = TRUE)

# Separate by data group
data1 <- x %>% filter(data.group == "u")  # Upregulated metabolites
data2 <- x %>% filter(data.group == "d")  # Downregulated metabolites


#---------------------------------------
# Prepare lists for enrichment analysis
#---------------------------------------
d <- data1
d <- data2

d_data_u <- NULL
d_data_d <- NULL

for(i in 1:nrow(d)){
  
  index1 <- d[i,17]  # Number of significant metabolites
  index2 <- d[i,18]
  
  x1 <- rep(d[i,3], index1)  # Comparison group name
  y1 <- rep(d[i,3], index2)
  
  x2 <- unlist(strsplit(d[i,31], "\\|"))  # Metabolite names
  y2 <- unlist(strsplit(d[i,32], "\\|"))
  
  w <- as.data.frame(x2)
  z <- as.data.frame(y2)
  
  # Combine group names and metabolites
  w <- cbind(x1, w)
  z <- cbind(y1, z)
  
  d_data_u[i] <- list(w)
  d_data_d[i] <- list(z)
}

category_a <- d_data_u  
category_b <- d_data_d


#---------------
# Merge lists
#----------------
Annotation_a <- NULL
Annotation_b <- NULL

for(i in 1:length(category_a)){
  Annotation_a <- rbind(Annotation_a, category_a[[i]])
  Annotation_b <- rbind(Annotation_b, category_b[[i]])
}

#---------
# Export
#----------
write.csv(Annotation_a, "annotaionlist_up.csv", fileEncoding = "CP932")
write.csv(Annotation_b, "annotaionlist_down.csv", fileEncoding = "CP932")
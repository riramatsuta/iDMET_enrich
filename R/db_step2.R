rm(list=ls(all=TRUE))

load(file = "ALL.RData")       # Load dataset

# ------------------------------
#   Split into each comparison
# ------------------------------
k <- 1
FIN <- NULL
CompName <- NULL
FNAME <- NULL
for(i in 1:length(ALL)){
  n <- ncol(ALL[[i]])
  for(j in 3:n){
    #for(j in 6:n){
    CompName[k] <- names(ALL[[i]])[j]
    FIN[[k]] <- data.frame(name = as.character(ALL[[i]][,1]), ratio = as.numeric(as.character(ALL[[i]][,j])))
   
    k <- k + 1
  }
}

names(FIN) <- CompName            # Assign comparison names
nf <- chartr(".", " vs ", names(FIN))  # Edit comparison names

save(FIN, file = "FIN.RData")     

# ------------------------------
#   Cluster (file) numbers
# ------------------------------
# nameR <- substr(names(ALL), 1, 12)   # Select first 8 chars of PMID
nameR <- substr(names(ALL), 1, 6)      # Select first 5 chars of PID

nn <- nameR

k <- 1
FNAME <- NULL
for(i in 1:length(ALL)){
  n <- ncol(ALL[[i]])
  for(j in 3:n){
    # for(j in 6:n){
    FNAME[k] <- nn[i]
    k <- k + 1
  }
}

gclass <- as.numeric(factor(FNAME))    # Assign numbers for each study/paper



#--------------------------
# Organize ID information
#--------------------------
nameR <- names(ALL)
id <- strsplit(nameR, "\\_")

pmid_sub <- NULL
for (i in 1:length(id)){
  x <- id[[i]][1]
  pmid_sub[i] <- x
}

#-----------------------
# Aggregate by category
#-----------------------
library(dplyr)
library(magrittr)

FP <- data.frame(FNAME) %>%
  group_by(FNAME) %>%
  mutate(num = row_number()) %>%
  ungroup() %>%
  dplyr::mutate(pn2 = paste(FNAME, "-", num, sep = ""))

# ---------------
# Save results
# ---------------
save(FNAME, file = "FNAME.RData")    
save(gclass, file = "gclass.RData")  
save(FP, file = "FP_set2.RData")   

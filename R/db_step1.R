rm(list=ls(all=TRUE))

# ----------------------
#   Data preparation
# ----------------------

file <- "./metabodic_ver4.csv"
K <- read.csv(file)[,c(1,2)]  # Metabolite dictionary (850 metabolites)

K <- read.csv(file, fileEncoding = "shift-jis")

# ----------------------
#   Load new data
# ----------------------
path <- "./data/csv_set4" 
L <- list.files(path)  # List of actual data files

# ----------------------
#   Matching with metabolite dictionary
# ----------------------
ALL <- NULL
#ALL_un <- NULL
k <- 1
for(i in 1:length(L)){
  #for(i in 1){
  file <- paste(path,L[i], sep ="/")
  D <- read.csv(file, fileEncoding = "shift-jis", header=TRUE)[-1,] 
  print(i)
  
  # Metabolite names
  M <- as.character(D[,1]) # Metabolite names from each project
  
  MT <- as.numeric(matrix(NA,length(M)))
  for(k in 1:length(M)){
    
    # Matching with dictionary
    for(j in 1:nrow(K)){
      a <- as.character(K[j,2])  
      b <- strsplit(a,";")[[1]] # Split entries
      for(l in 1:length(b)){
        if (tolower(M[k])==tolower(b[l])){  # Convert to lowercase for matching
          MT[k] <- as.character(K[j,1])
        }
      }
    }
  }
  
  index <- !is.na(MT)
  all <- cbind(MT,D)[index,]
  
  ALL[i] <- list(all)
  
}

names(ALL) <- L

#ALL <- ALL[-57] # Optional removal

# ----------------------
#   Save results
# ----------------------
save(ALL, file = "ALL.RData")

# ----------------------
#   Get dimensions of each dataset
# ----------------------
m <- NULL
n <- NULL
nn <- NULL
for(i in 1:length(ALL)){
  x <-  dim(ALL[[i]])
  m[i] <- list(x)
  
  n[i] <- as.numeric(m[[i]][[1]])  
  nn[i] <- as.numeric(m[[i]][[2]])-2
}

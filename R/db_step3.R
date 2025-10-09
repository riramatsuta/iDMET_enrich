rm(list=ls(all=TRUE))

# ----------------------
#   Load data
# ----------------------
load(file = "FIN.RData")     
load(file = "gclass.RData")   
load(file = "FNAME.RData")  

# ----------------------------------------------------------
#   Collect metabolites with/without significant change
# ----------------------------------------------------------
k <- 1
D <- NULL
E2 <- NULL
for(i in 1:length(FIN)){
  
  ratio <- as.numeric(as.character(FIN[[i]][,2]))
  metabo <- as.character(FIN[[i]][,1])  # List of analyzed metabolites
  
  
  # Thresholds
  #  index1 <- which(ratio >= 1.5)   # fold-change >= 1.5
  #  index2 <- which(ratio < (1/(1.5)))   # fold-change < 1/1.5
  index1 <- which(ratio >= 1.2)        
  index2 <- which(ratio < (1/(1.2))) 
  
  index <- which(ratio == 9999)
  index3 <- which(is.na(ratio))
  index3 <- c(index, index3)
  
  A <- as.character(FIN[[i]][index1,1])  
  B <- as.character(FIN[[i]][index2,1]) 
  
  E1 <- unique(as.character(FIN[[i]][index3,1])) # Metabolites with missing/unusual values
  e2 <- length(E1) # Count
  E2[k] <- e2
  
  D[[i]] <- list(A, B, metabo, E1)
  
  k <- k + 1                 
}

#-------------------------
# Summarize counts
#-------------------------
w <- NULL
z <- NULL
for(i in 1:length(D)){
  x <- length(D[[i]][[1]])
  y <- length(D[[i]][[2]])
  w[i] <- as.numeric(as.data.frame(x))
  z[i] <- as.numeric(as.data.frame(y))
  d_data <- cbind(w,z) 
}

colnames(d_data) <- c("metabo", "NA")


# ----------------------
#   Compute odds ratios
# ----------------------
k <- 1
i_index <- NULL; j_index <- NULL; odd <- NULL; pv <- NULL;  i_pmid <- NULL;  j_pmid <- NULL
A_metabolite <- NULL; B_metabolite <- NULL; C_metabolite <- NULL; D_metabolite <- NULL
an <- NULL; bn <- NULL; cn <- NULL; dn <- NULL; M_metabolite <- NULL; mn <- NULL; i_name <- NULL; j_name <- NULL

for(i in 1:length(D)){
  for(j in 1:length(D)){
    
    print(i)
    
    # Cross-tabulation
    x1 <- unique(D[[i]][[3]]) # Upregulated metabolites in group i
    x2 <- unique(D[[i]][[4]]) # Downregulated/missing metabolites in group i
    
    y1 <- unique(D[[j]][[3]])
    y2 <- unique(D[[j]][[4]])
    
    a <- sum(!is.na(match(x1, y1))) # up/up
    b <- sum(!is.na(match(x1, y2))) # up/down
    c <- sum(!is.na(match(x2, y1))) # down/up
    d <- sum(!is.na(match(x2, y2))) # down/down
    
    # Store counts
    an[k] <- a
    bn[k] <- b
    cn[k] <- c
    dn[k] <- d
    
    # Common metabolites between two groups
    x3 <- unique(D[[i]][[5]])
    y3 <- unique(D[[j]][[5]])
    
    m <- sum(!is.na(match(x3, y3)))
    mn[k] <- m
    
    index <- match(x3, y3)
    index <- index[!is.na(index)]
    M_metabolite[k] <- paste0(y3[index], collapse = "|")
    #-----------------------------------------------
    
    # Common significant metabolites in each category
    # (a) up/up
    index <- match(x1, y1)
    index <- index[!is.na(index)]
    A_metabolite[k] <- paste0(y1[index], collapse="|")
    
    # (b) up/down
    index <- match(x1, y2)
    index <- index[!is.na(index)]
    B_metabolite[k] <- paste0(y2[index], collapse="|")
    
    # (c) down/up
    index <- match(x2, y1)
    index <- index[!is.na(index)]
    C_metabolite[k] <- paste0(y1[index], collapse="|")
    
    # (d) down/down
    index <- match(x2, y2)
    index <- index[!is.na(index)]
    D_metabolite[k] <- paste0(y2[index], collapse="|")
    
    # Compute odds ratio (add 0.5 if zero count)
    if (a*b*c*d == 0) {
      a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5
    }
    
    # Chi-squared test
    x <- matrix(c(a, c, b, d), nrow=2)
    chi <- chisq.test(x)
    
    # Store results
    i_index[k] <- i
    j_index[k] <- j
    
    i_pmid[k] <- FNAME[i]
    j_pmid[k] <- FNAME[j] 
    
    i_name[k] <- nf[i]
    j_name[k] <- nf[j]
    
    odd[k] <- log2(a*d/(b*c))
    pv[k] <- chi$p.value
    
    k <- k + 1
  }
}

# ---------------------------
#   Organize results
# ---------------------------
R <- data.frame(
  i = i_index, j = j_index, odds = odd, p_value = pv, match_count = mn,
  a_count = an, b_count = bn, c_count = cn, d_count = dn, m = M_metabolite,
  a = A_metabolite, b = B_metabolite, c = C_metabolite, d = D_metabolite,
  ip = i_pmid, jp = j_pmid, iname = i_name, jname = j_name
)

# Remove duplicate combinations
Rsub <- subset(R, i < j)

# Export results
write.csv(Rsub, file = "Rsub.csv")

# Generate OR and P matrices
OR <- matrix(odd, nrow = length(D), byrow = T)
P <- matrix(pv, nrow = length(D), byrow = T)

# ----------------------------------------
#   Explore pairs with high odds ratios
# -----------------------------------------
OR2 <- abs(OR)

z_all <- NULL
for (i in 1:nrow(OR)){
  or2 <- abs(OR2[i,])
  
  index <- (or2 > 4)
  index[i] <- FALSE
  
  if (length(which(index)) > 0){
    z <- cbind(i, which(index))
    z_all <- rbind(z_all, z)
  }
}

pair_odds <- cbind(z_all, OR[z_all])

# ---------------------------
# Graph adjacency matrix
# ---------------------------
library(igraph)

Q <- P
Q[P < 0.05] <- 1     # Significant
Q[P >= 0.05] <- 0    # Not significant
diag(Q) <- 0

row.names(P) <- 1:1111  

# -------------------------
#   Graph clustering
# -------------------------
g <- graph.adjacency(Q*OR, weighted = TRUE)
V(g)$name <- c(1:111) # Required names  


g4 <- decompose.graph(g) 

# ----------------
#   Network 11
# ----------------
g11 <- g4[[2]]

E(g11)$color[E(g11)$weight >= 0] <- 'black'
E(g11)$color[E(g11)$weight < 0] <- 'red'

g11_c <- g11
E(g11_c)$weight <- 1

weight <- E(g11)$weight
weight[weight < -10] <- -10
weight[weight > 10] <- 10

mm <- V(g11)$name

#----------------------
# Output network list
#----------------------
library(dplyr)
ne <- get.data.frame(g11)
networklist <- subset(ne, from < to)


# ----------------
# Save results
# ----------------
# save(g11, file = "C:/R/g11.RData")
# save(D, file = "C:/R/D.RData")
# save(P, file = "C:/R/P.RData")

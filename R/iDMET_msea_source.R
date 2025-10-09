#   MSEA Source Code

# rm(list = ls(all = TRUE))
# ----------------------
# Load input data 
# ----------------------
INDATA <- read.csv(data, fileEncoding = "CP932", na.strings = c("", "N.A.", "N.D."))  
HMT_list <- read.csv(file = "./data/conpound_list.csv", fileEncoding = "CP932")

# ----------------------
# Select metabolites
# ----------------------
INDATA2 <- INDATA %>% filter(p.value < 0.05) # threshold

SIG <- as.character(INDATA2[,1])  # Significant metabolites
ALL <- as.character(INDATA[,1])   # All metabolites

# --------------------------
# Load enrichment file
# --------------------------
Pathway_list <- read.csv(iDMET_pathwaylist_ori, fileEncoding = "CP932")  

Pathway_HMT <- as.character(Pathway_list[,4])     
Pathway_class <- as.character(Pathway_list[,3])   
Uni_pathway <- unique(Pathway_class)


# ----------
# MSEA 
# ----------
NumM <- NaN; NumAll <- NaN; P_chi <- NaN; P_fish <- NaN; N <- NaN; N_all <- NaN; Nset <- NaN; W <- NaN; N_ori <- NaN
R <- NaN; M_metabolite <- NaN  


for (i in 1:length(Uni_pathway)){
  print(i)
  
  # Match pathway names
  m <- match(Pathway_class, Uni_pathway[i], 0)
  
  meta_set <- Pathway_HMT[m != 0]
  N_ori[i] <- length(meta_set)
  
  m0 <- match(ALL, meta_set, 0)
  meta <- meta_set[m0]
  n <- length(meta)
  NumAll[i] <- n
  
  sigM <- match(SIG, meta, 0)
  NumM[i] <- sum(sigM != 0)
  
  ns <- length(ALL) - length(SIG)
  
  a1 <- NumM[i]            # Significant & in pathway
  a2 <- NumAll[i] - a1     # Non-significant & in pathway
  a3 <- length(SIG) - a1   # Significant & not in pathway
  a4 <- ns - a2            # Non-significant & not in pathway
  tab <- t(matrix(c(a1, a2, a3, a4), 2))
  
  # Fisher's exact test
  resfish <- fisher.test(tab, alternative = "greater")
  P_fish[i] <- resfish$p.value
  
  index <- sigM[sigM != 0]
  M_metabolite[i] <- paste0(meta[index], collapse = "|")
}

# Adjust p-values using BH method
Q <- p.adjust(P_fish, method = "BH")


# ----------------------
#   Output results
# ----------------------
P <- cbind(Uni_pathway, N_ori, NumAll, NumM, P_fish, Q, M_metabolite)
P <- data.frame(P)
colnames(P) <- c("Metabolic_Pathway", "Total_Metabolites", "Detected_Metabolites", "Selected_Metabolites", 
  "p_value(Fisher's exact test)", "q_value(BH method)", "Selected_Metabolites_Names"
)

write.csv(P, file = out_file, fileEncoding = "CP932", row.names = FALSE)

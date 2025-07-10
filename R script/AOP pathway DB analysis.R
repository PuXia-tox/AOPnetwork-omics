# ===========================
# Pathway Perturbation Analysis
# ===========================

# ---- Path Setup ----
db_path <- '/rds/projects/x/xiap-xia-transcriptomics/AOP network/Database/'
deg_path <- '/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results/DEG/'
output_path <- '/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results/'
perturbed_path <- file.path(output_path, 'perturbed_pathways')
pdf_path <- file.path(output_path, 'Figure SI-pathway HP.pdf')

# ---- Load Packages ----
library(msigdbr)
library(magrittr)
library(ggrepel)
library(readxl)
library(pheatmap)

# ---- Load Data ----
Ampliseq <- read.csv(file.path(db_path, 'Zebrafish 1637 Reduced Genome.csv'))
path.db <- readRDS(file.path(db_path, 'RZT_AOP_KE_Pathway.rds'))
wp <- path.db[[1]] %$% .[entrez_gene %in% Ampliseq[, 4], ]
wp_list <- unique(as.character(wp[, 1]))
l <- length(wp_list)

# ---- DEG Files ----
sample_names <- list.files(deg_path)
condition_names <- unlist(strsplit(sample_names, '.csv'))
DEG_list <- lapply(sample_names, function(f) read.csv(file.path(deg_path, f)))
names(DEG_list) <- sample_names

# ---- Molecular Weight Adjustment ----
nanoMW <- setNames(c(150.71, 81.38, 86.9368, 79.866), c('NP1','NP10','NP5','NP9'))

# ---- Initialize Matrix ----
CMap_RHT <- matrix(10^6, nrow = l, ncol = length(sample_names),
                   dimnames = list(wp_list, condition_names))

# ---- Enrichment Loop ----
for (mmm in seq_along(DEG_list)) {
  DEG <- DEG_list[[mmm]][, c(1, 14)]
  colnames(DEG)[1] <- colnames(Ampliseq)[3]
  DEG <- merge(DEG, Ampliseq[, c(3, 4)], by = 'Associated.Gene.Name')
  DEG <- DEG[DEG[, 2] > 0, ]
  DEG[, 2] <- log10(DEG[, 2])
  
  # Dose adjustment
  sample_tag <- names(nanoMW)[grepl(names(nanoMW), sample_names[mmm])]
  if (length(sample_tag)) {
    DEG[, 2] <- log10(10^DEG[, 2] * 1000 / nanoMW[sample_tag])
  }
  
  # Per-pathway calculations
  my_family <- matrix(NA, nrow = l, ncol = 6,
                      dimnames = list(wp_list, c('Num', 'POD', 'SD', '%', 'Genes', 'KE')))
  l.mat <- matrix(NA, nrow = l, ncol = 2)
  
  for (i in seq_along(wp_list)) {
    a <- which(wp[, 1] == wp_list[i])
    b <- which(DEG[, 3] %in% wp[a, 2])
    
    l.mat[i, ] <- c(length(b), length(a))
    x <- DEG[b, 2]
    
    my_family[i, 'Num'] <- length(b)
    my_family[i, 'POD'] <- 10^mean(x)
    my_family[i, 'SD']  <- sd(x) / length(b)
    my_family[i, '%']   <- round(100 * length(b) / length(a), 2)
    my_family[i, 'Genes'] <- paste(DEG[b, 1], collapse = ', ')
    my_family[i, 'KE'] <- paste(path.db[[2]][path.db[[2]][, 2] == wp_list[i], 1], collapse = '; ')
  }
  
  # Filter and save
  valid_rows <- which(l.mat[, 1] >= 3 | (l.mat[, 1] == 2 & l.mat[, 2] < 15))
  my_family <- my_family[valid_rows, , drop = FALSE]
  
  if (nrow(my_family) >= 3) {
    valid_pathways <- path.db[[2]][, 2]
    my_family <- my_family[rownames(my_family) %in% valid_pathways, , drop = FALSE]
    
    write.csv(my_family[order(as.numeric(my_family[, 'POD'])), ],
              file.path(perturbed_path, sample_names[mmm]))
    
    CMap_RHT[rownames(my_family), mmm] <- as.numeric(my_family[, 'POD'])
  }
}

# ---- Final Matrix Cleanup ----
valid_rows <- which(rowSums(CMap_RHT == 10^6) != ncol(CMap_RHT))
CMap_RHT <- log10(CMap_RHT[valid_rows, ])
colnames(CMap_RHT) <- condition_names

# ---- Visualization ----
pdf(file = pdf_path, width = 10, height = 10)
pheatmap(CMap_RHT, show_rownames = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("red", "white"))(100))
dev.off()
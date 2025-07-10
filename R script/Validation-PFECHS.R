# ===========================
# Pathway Perturbation Analysis
# ===========================

# ---- Path Setup ----
db_path <- '/rds/projects/x/xiap-xia-transcriptomics/AOP network/Database/'
deg_path <- '/rds/projects/x/xiap-xia-transcriptomics/AOP network/Cross_Validation/DEG/'
output_path <- '/rds/projects/x/xiap-xia-transcriptomics/AOP network/Cross_Validation/AOP/'
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
DEG_list <- list()
DEG_list[[1]] <- read_excel(paste0(deg_path,'PFECHS.xlsx'),  na='NA') %>% as.matrix()
sample_names <- c('PFECHS')
names(DEG_list) <- sample_names

# ---- Enrichment Loop ----
mmm <- 1
DEG <- DEG_list[[mmm]][, c(1, 14)]
DEG <- apply(DEG, 2, as.numeric)
DEG[, 2] <- log10(DEG[, 2])
if(sum(is.na(DEG[, 2]))){
  DEG <- DEG[!is.na(DEG[, 2]), ]
}

# Per-pathway calculations
my_family <- matrix(NA, nrow = l, ncol = 6,
                    dimnames = list(wp_list, c('Num', 'POD', 'SD', '%', 'Genes', 'KE')))
l.mat <- matrix(NA, nrow = l, ncol = 2)

for (i in seq_along(wp_list)) {
  
  a <- which(wp[, 1] == wp_list[i])
  b <- which(as.character(DEG[,1]) %in% wp[a, 2])
  
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
valid_rows <- which(l.mat[, 1] >= 1)
my_family <- my_family[valid_rows, , drop = FALSE]

if (nrow(my_family) >= 3) {
  valid_pathways <- path.db[[2]][, 2]
  my_family <- my_family[rownames(my_family) %in% valid_pathways, , drop = FALSE]
  
  write.csv(my_family[order(as.numeric(my_family[, 'POD'])), ],
            file.path('/rds/projects/x/xiap-xia-transcriptomics/AOP network/Cross_Validation/Pathway/', sample_names[mmm]))
}

# nrow(my_family) == 0
# no matched pathways

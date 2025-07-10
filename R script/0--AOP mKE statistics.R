# ===========================================
# Pathway Perturbation Analysis - Neuro AOP
# ===========================================

## ---- Setup Paths ----
db_path        <- "/rds/projects/x/xiap-xia-transcriptomics/AOP network/Database/"
output_path    <- "/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results/"
perturbed_path <- file.path(output_path, "perturbed_pathways")
pdf_path       <- file.path(output_path, "Figure SI-pathway HP.pdf")

## ---- Load Packages ----

library(readr)
library(dplyr)
library(networkD3)
library(magrittr)

## ---- Load Input Data ----
setwd(db_path)

all_path_array <- readRDS("all MIE to AO paths-updated 2020.rds")
Ampliseq       <- read.csv("Zebrafish 1637 Reduced Genome.csv")
path.db        <- readRDS("RZT_AOP_KE_Pathway.rds")

all_path_array[which(sapply(all_path_array,length)==1)]
length(which(sapply(all_path_array,length)>=3))
length(which(sapply(all_path_array,length)==max(sapply(all_path_array,length))))
length(which(sapply(all_path_array,length)==3))
summary(sapply(all_path_array,length))
pdf(file.path(output_path, "hist_Length of LAOP.pdf"), width = 12, height = 12)
hist(sapply(all_path_array,length),breaks = 30,cex.lab = 1.5,
     xlab = 'Length of LAOP',main = '',
     ylab = 'Frequency')
dev.off()

library(nortest)
path_lengths <- sapply(all_path_array, length)
ad.test(path_lengths)

mie <- table(sapply(all_path_array,function(x)x[1]))
getAOPNodeName(steatosis_aop, names(which.max(mie)))

mie_table <- table(sapply(all_path_array,function(x)x[1]))
# Sort and extract top 5 MIEs
top_mie_ids <- names(sort(mie_table, decreasing = TRUE)[1:5])
top_mie_counts <- sort(mie_table, decreasing = TRUE)[1:5]
# Get human-readable node names (assuming getAOPNodeName returns that)
top_mie_names <- sapply(top_mie_ids, function(id) getAOPNodeName(steatosis_aop, id))
# Combine results into a data frame for clarity
top_mie_df <- data.frame(
  NodeID = top_mie_ids,
  NodeName = top_mie_names,
  Frequency = as.integer(top_mie_counts)
)
print(top_mie_df)


mie_table <- table(sapply(all_path_array,function(x)x[length(x)]))
# Sort and extract top 5 MIEs
top_mie_ids <- names(sort(mie_table, decreasing = TRUE)[1:5])
top_mie_counts <- sort(mie_table, decreasing = TRUE)[1:5]
# Get human-readable node names (assuming getAOPNodeName returns that)
top_mie_names <- sapply(top_mie_ids, function(id) getAOPNodeName(steatosis_aop, id))
# Combine results into a data frame for clarity
top_mie_df <- data.frame(
  NodeID = top_mie_ids,
  NodeName = top_mie_names,
  Frequency = as.integer(top_mie_counts)
)
print(top_mie_df)



mie <- table(sapply(all_path_array,function(x)x[1]))
getAOPNodeName(steatosis_aop, names(which.max(mie)))

wp     <- path.db[[1]] %$% .[entrez_gene %in% Ampliseq[, 4], ]
wp_list <- unique(as.character(wp[, 1]))
mke    <- path.db[[2]] %$% .[X.1 %in% wp_list, ]
mke[,1] <- toupper(gsub(",", "", mke[,1]))

## ---- Ontology & Annotation ----
AOP_ontology <- suppressMessages(read_excel("AO_wiki_downloads.xlsx", sheet = "Sheet1") %>% as.matrix())
AOP_ontology[,2] <- sapply(strsplit(AOP_ontology[,2], ":"), function(x) x[2])

AOP_keID <- suppressMessages(read_excel("AO_wiki_downloads.xlsx", sheet = "Sheet2") %>% as.matrix())
AOP_keID <- AOP_keID[!is.na(AOP_keID[,1]),]
AOP_keID[,1] <- trimws(AOP_keID[,1])

AOP <- merge(AOP_ontology, AOP_keID, by.x = "key event id", by.y = "ID", all.x = TRUE)
KE_AOP_annotation <- suppressMessages(read_excel("KE-annotation_Shanshan.xlsx") %>% as.matrix())
KE_AOP_annotation_2 <- KE_AOP_annotation[,c(2,4,6)]

AOP_2 <- AOP[,c(1,2,14)]

AOP_RZT <- merge(KE_AOP_annotation_2,AOP_2,by.x = 'AOP ID', by.y = 'key event id',all.x = TRUE)
AOP_RZT <- AOP_RZT[,c(2,5,3)]
colnames(AOP_RZT) <- c('KE ID','AOP ID','Molecular Pathways')

df_aop_rzt <- AOP_RZT


# --- 2. Handle missing values and convert data types ---
# Drop rows where 'AOP ID' or 'KE ID' are NA, as these are essential for relationships
df_cleaned <- df_aop_rzt %>%
  filter(!is.na(`AOP ID`) & !is.na(`KE ID`)) %>%
  mutate(
    `KE ID` = as.character(`KE ID`), # Convert KE ID to integer then character for cleaner labels
    `AOP ID` = as.character(`AOP ID`),           # Ensure AOP ID is character
    `Molecular Pathways` = as.character(`Molecular Pathways`) # Ensure Molecular Pathways is character
  )
df_cleaned <- df_cleaned %>%
     distinct(.keep_all = TRUE)

cat("--- Individual Column Counts ---\n")


#number of 1-2-1 GO-KE match
sum(table(df_cleaned[,3])==1) # AOP-RZT
round(sum(table(df_cleaned[,3])==1)/length(table(df_cleaned[,3]))*100,2)
pdf(file.path(output_path, "hist_AOPwiki.pdf"), width = 12, height = 12)
hist(table(df_cleaned[,3]),breaks = 200,cex.lab = 1.5,
     xlab = 'Number of covered KEs',main = '',
     ylab = 'Number of molecular pathways')
dev.off()


sum(table(AOP_ontology[,5])==1) # AOP-wiki
round(sum(table(AOP_ontology[,5])==1)/length(table(AOP_ontology[,5]))*100,2)
pdf(file.path(output_path, "hist_AOPwiki.pdf"), width = 12, height = 12)
hist(table(AOP_ontology[,5]),breaks = 200,cex.lab = 1.5,
     xlab = 'Number of covered KEs',main = '',
     ylab = 'Number of molecular pathways')
dev.off()




# Basic graph stats
## ---- Load Network ----
steatosis_aop     <- convert_cytoscape_to_aop("aopwiki-KER_2020_09_21.cyjs")
aop_graph_matrix  <- convert_aop_to_graph(steatosis_aop) %>% as("matrix")
aop_graph         <- graph.adjacency(aop_graph_matrix, mode = "directed")

num_nodes <- vcount(aop_graph)
num_edges <- ecount(aop_graph)

# Node degree centrality
deg <- igraph::degree(aop_graph)  # total degree
top_central <- names(sort(deg, decreasing = TRUE)[1:10])  # top 10 most connected mKEs
for(i in 1:5){
  print(tolower(getAOPNodeName(steatosis_aop,top_central[i])))
}

# If you have pathway info in vertex attributes, e.g. a list of pathway IDs per node:
# aop_graph <- set_vertex_attr(aop_graph, "pathways", value = list_of_pathways)
# Replace list_of_pathways with your actual attribute source

# Node with highest pathway coverage (assuming 'pathways' is a list or vector)
if ("pathways" %in% vertex_attr_names(aop_graph)) {
  pathway_counts <- sapply(V(aop_graph)$pathways, function(x) length(unique(x)))
  top_pathway_nodes <- names(sort(pathway_counts, decreasing = TRUE)[1:10])
} else {
  top_pathway_nodes <- NA  # Pathway coverage data not found
}

# Biological category mapping (assuming custom vertex attribute exists, e.g. "category")
if ("category" %in% vertex_attr_names(aop_graph)) {
  bio_coverage <- table(V(aop_graph)$category)
} else {
  bio_coverage <- NA  # No biological categories available
}

# Final report
cat(paste0("📊 The graph contains ", num_nodes, " nodes and ", num_edges, " directed edges.\n"))
cat("🏆 Most centered mKEs (by degree):\n")
print(top_central)

if (!is.na(top_pathway_nodes[1])) {
  cat("🧬 mKEs with highest pathway coverage:\n")
  print(top_pathway_nodes)
}

if (!is.na(bio_coverage[1])) {
  cat("🧠 Biological coverage of core/key categories:\n")
  print(bio_coverage)
}

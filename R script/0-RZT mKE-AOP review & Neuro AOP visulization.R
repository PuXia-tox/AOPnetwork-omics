# ===========================================
# Pathway Perturbation Analysis - Neuro AOP
# ===========================================

## ---- Setup Paths ----
db_path        <- "/rds/projects/x/xiap-xia-transcriptomics/AOP network/Database/"
output_path    <- "/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results/"
perturbed_path <- file.path(output_path, "perturbed_pathways")
pdf_path       <- file.path(output_path, "Figure SI-pathway HP.pdf")

## ---- Load Packages ----
library(msigdbr)
library(magrittr)
library(ggrepel)
library(readxl)
library(pheatmap)
library(stringdist)
library(aop)
library(igraph)
library(graph)

## ---- Load Input Data ----
setwd(db_path)

all_path_array <- readRDS("all MIE to AO paths-updated 2020.rds")
Ampliseq       <- read.csv("Zebrafish 1637 Reduced Genome.csv")
path.db        <- readRDS("RZT_AOP_KE_Pathway.rds")

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

AOP_match <- AOP[AOP[,1] %in% KE_AOP_annotation[,6], ]
x         <- unique(AOP_match[,14])
AOP_neuro_ke <- toupper(read.csv('neuro_KE.csv')[,-1])
AOP_neuro_ao <- suppressMessages(read_excel("Neurotoxicological AOPs.xlsx") %>% as.matrix())

y <- KE_AOP_annotation[which(toupper(KE_AOP_annotation[,2]) %in% AOP_neuro_ke),2] %>% toupper()
AOP_neuro_ao <- na.omit(unique(as.vector(AOP_neuro_ao[, -c(1,2)]))) %>% toupper() %>% gsub("[(),]", "", .)

y.merge <- merge(KE_AOP_annotation[which(toupper(KE_AOP_annotation[,2]) %in% AOP_neuro_ke),c(2,6)],
                 AOP[,c(1,14)],
                 by.x='AOP ID',by.y='key event id',all.x = TRUE
)

df_no_duplicates <- y.merge %>%
  distinct(.keep_all = TRUE)

setdiff(AOP_neuro_ke, toupper(KE_AOP_annotation[,2]))





## ---- Load Network ----
steatosis_aop     <- convert_cytoscape_to_aop("aopwiki-KER_2020_09_21.cyjs")
aop_graph_matrix  <- convert_aop_to_graph(steatosis_aop) %>% as("matrix")
aop_graph         <- graph.adjacency(aop_graph_matrix, mode = "directed")

nodes <- convert_aop_to_graph(steatosis_aop) %>%
  as("matrix") %>%
  colnames() %>%
  sapply(function(x) getAOPNodeName(steatosis_aop, x)) %>%
  gsub("\\(.*\\)", "", .) %>%
  gsub(",", "", .) %>%
  toupper() %>%
  cbind(matrix(.), names(.))


#set.seed(read.csv("aop network seed.csv")[, 2])

## ---- KE Node Matching ----
ke.l <- which(nodes[,2] %in% y)
ke.l <- union(ke.l, c(399,1107))
setdiff(y,nodes[,2])
diff_data <- 'NEUROINFLAMMATION'
grep(diff_data,nodes[,2])
nodes[grep(diff_data,nodes[,2]),]

ke.neru_ao <- which(nodes[,2] %in% AOP_neuro_ao) # 22 neuro, 16 can be presented in current aop network
ke.neru_ao <- union(ke.neru_ao, c(755, 306, 756, 817, 839))

ke_node           <- na.omit(rownames(aop_graph_matrix)[ke.l])
ke_neuro_ao_node  <- rownames(aop_graph_matrix)[ke.neru_ao]

## ---- Filter Valid Paths ----
valid_paths <- Filter(function(path) {
  any(path %in% ke_node) && any(path %in% ke_neuro_ao_node)
}, all_path_array)

## ---- Edge Extraction from Shortest Paths ----
highlight_edge_ids <- c()
for (i in ke_node) {
  for (j in ke_neuro_ao_node) {
    paths <- all_shortest_paths(aop_graph, from = i, to = j, mode = "out")$vpaths
    for (p in paths) {
      path_nodes <- V(aop_graph)$name[p]
      if (length(path_nodes) > 1) {
        for (k in seq_len(length(path_nodes) - 1)) {
          edge_id <- get.edge.ids(aop_graph, vp = c(path_nodes[k], path_nodes[k+1]))
          if (edge_id != 0) highlight_edge_ids <- c(highlight_edge_ids, edge_id)
        }
      }
    }
  }
}
highlight_edge_ids <- unique(highlight_edge_ids)

## ---- Node Annotations ----
Wiki_MIE  <- read_excel("Wiki-MIE.xlsx") %>% as.matrix()
non_MIE   <- read_excel("non-MIE.xlsx") %>% as.matrix()
Wiki_AO   <- read_excel("Wiki-AO.xlsx") %>% as.matrix()

mie_clean <- toupper(gsub("[(),]", "", unique(Wiki_MIE[,2][!Wiki_MIE[,2] %in% non_MIE[,4]])))
ao_clean  <- toupper(gsub("[(),]", "", unique(unlist(Wiki_AO[,4:7][!is.na(Wiki_AO[,4:7])]))))

mie.l     <- which(nodes[,1] %in% mie_clean)
aop.l     <- which(nodes[,1] %in% ao_clean)

## ---- Layout & Plotting ----
# layout_coords <- layout.fruchterman.reingold(aop_graph)
# V(aop_graph)$plotX <- layout_coords[,1]
# V(aop_graph)$plotY <- layout_coords[,2]
# layout_coords <- layout.fruchterman.reingold(aop_graph,dim=3)
# layout_coords <- layout_coords * 2
# saveRDS(layout_coords, file = file.path(output_path, "AOP_layout_coords.rds"))
layout_coords <- readRDS(file.path(output_path, "AOP_layout_coords.rds"))

# Node styles
vertex_sizes <- rep(3, gsize(aop_graph))
vertex_sizes[c(ke.l,ke.neru_ao)] <- 5

vSubCol <- rep("white", gsize(aop_graph))
vSubCol[c(ke.l,ke.neru_ao)]       <- "lightblue"
vSubCol[mie.l]      <- "green"
vSubCol[aop.l]      <- "red"
#vSubCol[ke.neru_ao] <- "yellow"

# Edge styles
e_col   <- rep("gray", ecount(aop_graph))
e_width <- rep(0.5, ecount(aop_graph))
e_col[highlight_edge_ids]   <- "purple"
e_width[highlight_edge_ids] <- 3

## ---- Save PDF Plot ----
pdf(file.path(output_path, "Neuro_AOP_Network_Overview.pdf"), width = 12, height = 12)
par(mar = c(0, 0, 0, 0))

plot(aop_graph,
     vertex.color        = vSubCol,
     vertex.frame.color  = "black",
     edge.color          = e_col,
     edge.width          = e_width,
     edge.arrow.size     = 0.2,
     vertex.size         = vertex_sizes,
     vertex.shape        = "circle",
     vertex.label        = NA,
     layout              = layout_coords,
     vertex.label.cex    = 0.4,
     vertex.label.family = "Helvetica")

legend("bottomleft",
       legend = c("KE", "MIE", "AO", "Highlighted Path"),
       col    = c("lightblue", "green", "red", "purple"),
       pch    = 19, pt.cex = 1.5)

dev.off()
library(aop)
library(igraph)
library(graph)
library(magrittr)


requireNamespace("readxl", quietly = TRUE)
requireNamespace("igraph", quietly = TRUE)
requireNamespace("aop", quietly = TRUE)
requireNamespace("graph", quietly = TRUE)

# Load AOP nodes and annotation information
setwd("/rds/projects/x/xiap-xia-transcriptomics/AOP network/Database")
steatosis_aop <- convert_cytoscape_to_aop("AOP-wiki.cyjs")
aop_graph_matrix <- convert_aop_to_graph(steatosis_aop) %>% as("matrix")

nodes <- colnames(aop_graph_matrix) %>%
  sapply(function(x) getAOPNodeName(steatosis_aop, x)) %>%
  gsub("\\(.*\\)", "", .) %>%
  gsub(",", "", .) %>%
  toupper() %>%
  cbind(matrix(.), names(.)) %>%
  .[, c(1, 2)]

aop_graph <- graph.adjacency(aop_graph_matrix, mode = "directed")
set.seed(read.csv("aop network seed.csv")[, 2])

# Load supporting data
Wiki_MIE     <- read_excel("Wiki-MIE.xlsx", na = "NA") %>% as.matrix()
non_MIE      <- read_excel("non-MIE.xlsx", na = "NA") %>% as.matrix()
Wiki_AO      <- read_excel("Wiki-AO.xlsx", na = "NA") %>% as.matrix()
match2AOP    <- read_excel("match to AOP wiki KE-Zebrafish.xlsx", na = "NA") %>% as.matrix()
AOP_ontology <- suppressMessages(read_excel("AOP_Ontology_Assignments.033017.xlsx", na = "NA") %>% as.matrix())
MIE          <- AOP_ontology[grep("GO", AOP_ontology[, 8]), ]

# Match key events
path.db <- readRDS("RZT_AOP_KE_Pathway.rds")
wp <- path.db[[2]]

# Identify KE node locations
ke.l <- which(nodes[, 2] %in% y)
ke.l <- c(ke.l,588,544,454,713,1153,506)

# reverse checking
setdiff(y, nodes[, 2])
nn <- 'RENAL'
nodes[grep(nn,nodes[, 2]),]
grep(nn,nodes[, 2])


if (length(ke.l) <= 3) {
  stop("Too few key events matched; not applicable for this analysis.")
}

# MIE nodes
mie_raw <- unique(Wiki_MIE[, 2])
mie_clean <- toupper(gsub("[(),]", "", mie_raw[!mie_raw %in% non_MIE[, 4]]))
mie.l <- which(nodes[, 1] %in% mie_clean)

# AO nodes
ao_raw <- unique(unlist(Wiki_AO[, 4:7]))
ao_clean <- toupper(gsub("[(),]", "", ao_raw[!is.na(ao_raw)]))
aop.l <- which(nodes[, 1] %in% ao_clean)

# Layout and primary AOP network plot
layout_coords <- layout.fruchterman.reingold(aop_graph)
V(aop_graph)$plotX <- layout_coords[, 1]
V(aop_graph)$plotY <- layout_coords[, 2]
plotLay <- cbind(V(aop_graph)$plotX, V(aop_graph)$plotY)


# Base size for all nodes
vertex_sizes <- rep(4, gsize(aop_graph))
# Increase size for KE nodes
vertex_sizes[ke.l] <- 8  # You can adjust this number as needed

vSubCol <- rep("white", gsize(aop_graph))
vSubCol[ke.l] <- "lightblue"
vSubCol[mie.l] <- "green"
vSubCol[aop.l] <- "red"

# Plot with modified sizes
pdf("/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results/Neuro_AOP_Network_Overview.pdf", width = 10, height = 10)
par(mar = c(0, 0, 0, 0))

plot(aop_graph,
     vertex.color = vSubCol, vertex.frame.color = "black", edge.color = "black",
     vertex.label = NA, vertex.shape = "circle",
     vertex.size = vertex_sizes,  # ← new sizing vector here
     edge.arrow.size = 0.1, edge.width = 0.5, layout = plotLay,
     vertex.label.cex = 0.5, vertex.label.family = "Helvetica")

dev.off()

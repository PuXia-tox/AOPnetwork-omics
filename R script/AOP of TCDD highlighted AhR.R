# Load packages
library(aop)
library(graph)
library(igraph)
library(magrittr)
requireNamespace("readxl", quietly = TRUE)
requireNamespace("igraph", quietly = TRUE)
requireNamespace("aop", quietly = TRUE)
requireNamespace("graph", quietly = TRUE)

# Define directories
db_dir <- "/rds/projects/x/xiap-xia-transcriptomics/AOP network/Database/"
results_dir <- "/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results"
perturbed_dir <- file.path(results_dir, "perturbed_pathways")
qaop_output_dir <- file.path(results_dir, "qAOP")

# Load scripts and data
setwd(db_dir)

message("Loading AOP path array...")
all_path_array <- readRDS("all MIE to AO paths-updated 2020.rds")  # load once
steatosis_aop <- convert_cytoscape_to_aop("aopwiki-KER_2020_09_21.cyjs")
aop_graph_matrix <- convert_aop_to_graph(steatosis_aop) %>% as("matrix")
aop_graph <- graph.adjacency(aop_graph_matrix, mode = "directed")

message("Parsing AOP network...")
steatosis_aop <- convert_cytoscape_to_aop("aopwiki-KER_2020_09_21.cyjs")
nodes <- convert_aop_to_graph(steatosis_aop) %>%
  as("matrix") %>%
  colnames() %>%
  sapply(function(x) getAOPNodeName(steatosis_aop, x)) %>%
  gsub("\\(.*\\)", "", .) %>%
  gsub(",", "", .) %>%
  toupper() %>%
  cbind(matrix(.), names(.)) %>%
  .[, c(1, 3)]

# Process perturbed pathway samples
setwd(perturbed_dir)
sample_files <- list.files()

cmKE_file <- list.files('/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results/cmKE/',
                        pattern = "\\.rds$", full.names = TRUE) # Use full.names=TRUE for complete paths
cmKE <- list()
cmKE <- sapply(cmKE_file, readRDS, simplify = FALSE) # Use simplify = FALSE to ensure a list is returned



# Adjust the range for batch processing if needed
  
  for(mmm in 37:39){
  
  sample_name <- sample_files[mmm]
  sample_path <- file.path(perturbed_dir, sample_name)
  message(paste("Processing:", sample_name))
  data <- read.csv(sample_path)[, c(1, 3)]
  
  # Load AOP nodes and annotation information
  
  nodes <- colnames(aop_graph_matrix) %>%
    sapply(function(x) getAOPNodeName(steatosis_aop, x)) %>%
    gsub("\\(.*\\)", "", .) %>%
    gsub(",", "", .) %>%
    toupper() %>%
    cbind(matrix(.), names(.))
  
  # Load supporting data
  setwd("/rds/projects/x/xiap-xia-transcriptomics/AOP network/Database")
  Wiki_MIE     <- read_excel("Wiki-MIE.xlsx", na = "NA") %>% as.matrix()
  non_MIE      <- read_excel("non-MIE.xlsx", na = "NA") %>% as.matrix()
  Wiki_AO      <- read_excel("Wiki-AO.xlsx", na = "NA") %>% as.matrix()
  match2AOP    <- read_excel("match to AOP wiki KE-Zebrafish.xlsx", na = "NA") %>% as.matrix()
  AOP_ontology <- suppressMessages(read_excel("AOP_Ontology_Assignments.033017.xlsx", na = "NA") %>% as.matrix())
  MIE          <- AOP_ontology[grep("GO", AOP_ontology[, 8]), ]
  
  # Match key events
  path.db <- readRDS("RZT_AOP_KE_Pathway.rds")
  wp <- path.db[[2]]
  data.new <- data[data[, 1] %in% wp[, 2], ]
  
  for (i in seq_len(nrow(data.new))) {
    s <- unique(wp[wp[, 2] == data.new[i, 1], 1])
    data.new[i, 1] <- s[1]
    if (length(s) > 1) {
      extras <- cbind(s[-1], rep(data.new[i, 2], length(s) - 1))
      colnames(extras) <- colnames(data.new)
      data.new <- rbind(data.new, extras)
    }
  }
  
  data.new <- data.new[order(data.new[, 2], decreasing = TRUE), ]
  
  # Align to AOP node IDs
  for (i in seq_len(nrow(data.new))) {
    s <- nodes[nodes[, 1] == data.new[i, 1], 2]
    data.new[i, 1] <- s[1]
    if (length(s) > 1) {
      extras <- cbind(s[-1], rep(data.new[i, 2], length(s) - 1))
      colnames(extras) <- colnames(data.new)
      data.new <- rbind(data.new, extras)
    }
  }
  
  # Identify KE_neuro nodes locations
  AOP_ontology <- suppressMessages(read_excel("AO_wiki_downloads.xlsx", sheet = "Sheet1") %>% as.matrix())
  AOP_ontology[,2] <- sapply(strsplit(AOP_ontology[,2], ":"), function(x) x[2])
  AOP_keID <- suppressMessages(read_excel("AO_wiki_downloads.xlsx", sheet = "Sheet2") %>% as.matrix())
  AOP_keID <- AOP_keID[!is.na(AOP_keID[,1]),]
  AOP_keID[,1] <- trimws(AOP_keID[,1])
  AOP <- merge(AOP_ontology, AOP_keID, by.x = "key event id", by.y = "ID", all.x = TRUE)
  KE_AOP_annotation <- suppressMessages(read_excel("KE-annotation_Shanshan.xlsx") %>% as.matrix())
  AOP_match <- AOP[AOP[,1] %in% KE_AOP_annotation[,6], ]
  x         <- unique(AOP_match[,14])
  
  ke_nodes <- unique(na.omit(as.character(data.new[, 1])))
  ke.l <- which(nodes[, 2] %in% ke_nodes)
  AOP_neuro <- suppressMessages(read_excel("Neurotoxicological AOPs.xlsx") %>% as.matrix())
  x <- x[as.numeric(AOP_neuro[,1])] %>% gsub("\\s*\\([^\\)]+\\)", "", .)
  AOP_match[,14] <- gsub("\\s*\\([^\\)]+\\)", "", AOP_match[,14])
  y <- KE_AOP_annotation[KE_AOP_annotation[,1] %in% AOP_match[AOP_match[,14] %in% x, 1], 2] %>% toupper()
  AOP_neuro_ao <- na.omit(unique(as.vector(AOP_neuro[, -c(1,2)]))) %>% toupper() %>% gsub("[(),]", "", .)
  ke.neru_ao <- which(nodes[,2] %in% AOP_neuro_ao)
  ke.neru_ao <- union(ke.neru_ao, c(755, 306, 756, 817, 839)) 
  ke_node           <- na.omit(rownames(aop_graph_matrix)[ke.l])
  ke_neuro_ao_node  <- rownames(aop_graph_matrix)[ke.neru_ao]
  
  AhR.l <- grep('AHR/ARNT',nodes[,2])

  ## ---- Filter Valid Paths ----
  cmke_match <- grep(paste0('s',unlist(strsplit(sample_name,'.csv'))),names(cmKE))
  qAOP_path <- sapply(cmKE[[cmke_match]],function(x)x$path)
  
  valid_paths <- Filter(function(path) {
    any(path %in% ke_node) && any(path %in% ke_neuro_ao_node)
  }, qAOP_path)
  
  highlight_edge_ids <- unique(unlist(
    lapply(valid_paths, function(path) {
      edge_ids <- c()
      path_nodes <- as.character(path)
      for (i in seq_len(length(path_nodes) - 1)) {
        eid <- get.edge.ids(aop_graph, vp = c(path_nodes[i], path_nodes[i + 1]))
        if (eid != 0) edge_ids <- c(edge_ids, eid)
      }
      edge_ids
    })
  ))
  
  highlight_edge_ids <- unique(highlight_edge_ids)

  
  if (length(ke.l) <= 3) {
    stop("Too few key events matched; not applicable for this analysis.")
  }
  
  ## ---- Node Annotations ----
  Wiki_MIE  <- read_excel("Wiki-MIE.xlsx") %>% as.matrix()
  non_MIE   <- read_excel("non-MIE.xlsx") %>% as.matrix()
  Wiki_AO   <- read_excel("Wiki-AO.xlsx") %>% as.matrix()
  
  mie_clean <- toupper(gsub("[(),]", "", unique(Wiki_MIE[,2][!Wiki_MIE[,2] %in% non_MIE[,4]])))
  ao_clean  <- toupper(gsub("[(),]", "", unique(unlist(Wiki_AO[,4:7][!is.na(Wiki_AO[,4:7])]))))
  
  mie.l     <- which(nodes[,1] %in% mie_clean)
  aop.l     <- which(nodes[,1] %in% ao_clean)
  
  ## ---- Layout & Plotting ----
  layout_coords <- readRDS("/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results/AOP_layout_coords.rds")
  
  # Node styles
  vertex_sizes <- rep(3, gsize(aop_graph))
  vertex_sizes[ke.l] <- 5
  vertex_sizes[AhR.l] <- 10
  
  vSubCol <- rep("white", gsize(aop_graph))
  vSubCol[ke.l]       <- "lightblue"
  vSubCol[mie.l]      <- "green"
  vSubCol[aop.l]      <- "red"
  vSubCol[ke.neru_ao] <- "yellow"
  vSubCol[AhR.l] <- 'orange'
  
  # Edge styles
  e_col   <- rep("gray", ecount(aop_graph))
  e_width <- rep(0.5, ecount(aop_graph))
  e_col[highlight_edge_ids]   <- "purple"
  e_width[highlight_edge_ids] <- 3
  
  # 📄 Save plots to PDF
  base_name_raw <- strsplit(sample_name, '.csv')[[1]]
  sanitized_base_name <- gsub("[^A-Za-z0-9_.-]", "_", base_name_raw)
  sanitized_base_name <- gsub("^[. ]+|[. ]+$", "", sanitized_base_name)
  if (nchar(sanitized_base_name) == 0) {
    sanitized_base_name <- "unnamed_sample" # Fallback name
  }
  
  # Construct the full PDF filename using the sanitized base name
  pdf_filename <- paste0('AOP_Network_AhR_', sample_name, '.pdf')
  pdf_filepath <- file.path('/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results', pdf_filename)
  
  pdf(pdf_filepath, width = 10, height = 10)
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
         legend = c("KE", 'AhR Activation',"Neuro_AO", "MIE", "AO", "Perturbed Path"),
         col    = c("lightblue",'orange', "yellow", "green", "red", "purple"),
         pch    = 19, pt.cex = 1.5)
  
  dev.off()
  
  gc()
}




cmKE_file <- list.files('/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results/cmKE/',
                        pattern = "\\.rds$", full.names = TRUE) # Use full.names=TRUE for complete paths
cmKE <- list()
cmKE <- sapply(cmKE_file, readRDS, simplify = FALSE) # Use simplify = FALSE to ensure a list is returned


## ---- Filter Valid Paths ----
highlight_edge_ids <- list()
qAOP_path <- list()
for(i in 1:length(cmKE)){
  print(i)
  cmke <- cmKE[[i]]
  qAOP_path[[i]] <- sapply(cmke,function(x)x$path)
}

valid_paths <- list()
for(i in 1:length(qAOP_path)){
  print(i)
  cmke <- qAOP_path[[i]]
  valid_paths[[i]] <- Filter(function(path) {
    any(path %in% ke_node) && any(path %in% ke_neuro_ao_node)
  }, cmke)
}

# Function to transform a single path (vector of node IDs) into a vector of "A->B" strings
transform_path_to_edges <- function(path_nodes) {
  # Ensure path_nodes is character type
  path_nodes <- as.character(path_nodes)
  
  # If path has less than 2 nodes, no edges can be formed
  if (length(path_nodes) < 2) {
    return(character(0)) # Return an empty character vector
  }
  
  # Create the "A->B" strings
  # This uses vectorized operations: paste combines corresponding elements
  # path_nodes[-length(path_nodes)] gets all elements except the last (A)
  # path_nodes[-1] gets all elements except the first (B)
  edges <- paste0(path_nodes[-length(path_nodes)], "->", path_nodes[-1])
  
  return(edges)
}

paths_as_edge_strings <- list()
for(i in 1:length(valid_paths)){
  print(i)
  paths_as_edge_strings[[i]] <- lapply(valid_paths[[i]], transform_path_to_edges)
}


edges_AOP <- sapply(paths_as_edge_strings,function(x) unique(unlist(x)))


# Load AOP nodes and annotation information

# Process perturbed pathway samples
setwd(perturbed_dir)
cm_samples <- paste0(unlist(strsplit(sapply(strsplit(cmKE_file,'_coords'),function(x)x[2]),'.rds')),'.csv')
data_samples <- list()
data_samples <- lapply(cm_samples, read.csv)

ke_nodes <- list()

setwd("/rds/projects/x/xiap-xia-transcriptomics/AOP network/Database")
steatosis_aop <- convert_cytoscape_to_aop("aopwiki-KER_2020_09_21.cyjs")
aop_graph_matrix <- convert_aop_to_graph(steatosis_aop) %>% as("matrix")

nodes <- colnames(aop_graph_matrix) %>%
  sapply(function(x) getAOPNodeName(steatosis_aop, x)) %>%
  gsub("\\(.*\\)", "", .) %>%
  gsub(",", "", .) %>%
  toupper() %>%
  cbind(matrix(.), names(.))

aop_graph <- graph.adjacency(aop_graph_matrix, mode = "directed")

# Match key events
path.db <- readRDS("RZT_AOP_KE_Pathway.rds")
wp <- path.db[[2]]


for( k in 1:length(data_samples)){
  data <- data_samples[[k]][, c(1, 3)]
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
  ke_nodes[[k]] <- unique(na.omit(as.character(data.new[, 1])))
}

network_names <- unlist(strsplit(sapply(strsplit(cmKE_file,'_coords'),function(x)x[2]),'.rds'))
highlight_edge_ids_list <- edges_AOP
names(highlight_edge_ids_list) <- network_names
enlarged_nodes_list <- ke_nodes
names(enlarged_nodes_list) <- network_names

# Initialize network_attributes list
num_networks <- length(network_names)
network_attributes <- vector("list", num_networks)
names(network_attributes) <- names(highlight_edge_ids_list) # Use names from one of the lists

# Loop to combine the attributes for each network
for (i in 1:num_networks) {
  network_name <- names(highlight_edge_ids_list)[i] # Get the name of the current network
  network_attributes[[network_name]] <- list(
    highlighted_edges = highlight_edge_ids_list[[i]],
    enlarged_nodes = enlarged_nodes_list[[i]]
  )
}


chemical_mapping_short <- data.frame(
  Number_or_Code = c("27#", "3#", "60#", "NP1", "NP5", "NP9", "NP10", "70#", "7#", "47#",
                     "40#", "4#", "34#", "21#", "17#", "16#", "BPA", "25#", "P8", "38#",
                     "41#", "TCDD", "BPAF", "BPE", "BPS"),
  Short_Name = c("CLP", "DEHP", "MFA", "SnO2", "MnO2", "TiO2", "ZnO", "CPD", "DFC",
                 "TCS", "BbF", "DZN", "BaP", "DiR", "TPP", "PPC", "BPA", "CPS",
                 "CBP", "PFOA", "TDCPP", "TCDD", "BPAF", "BPE", "BPS"),
  stringsAsFactors = FALSE
)
network_attributes <- network_attributes[-which(sapply(strsplit(names(network_attributes),'-'),function(x)x[2])!='1')]
num_networks <- length(network_attributes)
names(network_attributes) <- sapply(strsplit(names(network_attributes),'-'),function(x)x[1])

current_network_names <- names(network_attributes)
new_short_names <- character(length(current_network_names))
names(new_short_names) <- current_network_names # Keep original names as keys for direct lookup


for (i in seq_along(current_network_names)) {
  current_code <- current_network_names[i]
  
  # Find the matching row in our mapping table using 'Number_or_Code'
  match_row <- chemical_mapping_short[chemical_mapping_short$Number_or_Code == current_code, ]
  
  if (nrow(match_row) > 0) {
    new_short_names[i] <- match_row$Short_Name[1] # Take the first match if duplicates exist (shouldn't be in this case)
  } else {
    # Handle cases where a code might not be found in your mapping
    # This ensures no NA names are introduced. You might want to log these.
    new_short_names[i] <- current_code # Keep the original code if no short name is found
    warning(paste("No short name found for code:", current_code))
  }
}

# Apply the new short names to your network_attributes list
names(network_attributes) <- new_short_names



# Calculate Pairwise Similarities ---
similarity_matrix <- matrix(NA, nrow = num_networks, ncol = num_networks,
                            dimnames = list(names(network_attributes), names(network_attributes)))

for (i in 1:num_networks) {
  for (j in i:num_networks) {
    attrs1 <- network_attributes[[i]]
    attrs2 <- network_attributes[[j]]
    
    # Jaccard Similarity for Highlighted Edges
    common_edges <- intersect(attrs1$highlighted_edges, attrs2$highlighted_edges)
    union_edges <- union(attrs1$highlighted_edges, attrs2$highlighted_edges)
    jaccard_edge_sim <- ifelse(length(union_edges) == 0, 1,
                               length(common_edges) / length(union_edges))
    
    # Jaccard Similarity for Enlarged Nodes
    common_nodes <- intersect(attrs1$enlarged_nodes, attrs2$enlarged_nodes)
    union_nodes <- union(attrs1$enlarged_nodes, attrs2$enlarged_nodes)
    jaccard_node_sim <- ifelse(length(union_nodes) == 0, 1,
                               length(common_nodes) / length(union_nodes))
    
    avg_jaccard_sim <- (jaccard_edge_sim + jaccard_node_sim) / 2
    
    similarity_matrix[i, j] <- avg_jaccard_sim
    similarity_matrix[j, i] <- avg_jaccard_sim
  }
}


cat("\n--- Pairwise Similarity Matrix (Average Jaccard) ---\n")
print(round(similarity_matrix, 3))


# --- 4. Clustering Analysis ---

distance_matrix <- as.dist(1 - similarity_matrix)

cat("\n--- Distance Matrix (1 - Average Jaccard) ---\n")
print(round(distance_matrix, 3))

# Perform Hierarchical Clustering
hc_result <- hclust(distance_matrix, method = "complete")

cat("\n--- Hierarchical Clustering Result ---\n")
print(hc_result)

# --- DEBUGGING fviz_dend ISSUE ---
# Try plotting the basic dendrogram first to see its structure
pdf(file.path(output_path, "AOP_Network Similarity.pdf"), width = 12, height = 12)
par(mar = c(0, 0, 0, 0))

plot(hc_result,
     main = "Basic Dendrogram: Hierarchical Clustering of Network Visual Attributes",
     xlab = "Networks",
     ylab = "Distance")
dev.off()

# Then try fviz_dend with different 'k' values if 3 doesn't work.
# The error "3, 0" strongly suggests it couldn't find 3 distinct clusters to draw rectangles around.
# This often happens if the distances are all very similar or some networks are identical.


# You can also inspect the heights to see where natural cuts might occur
cat("\n--- Dendrogram Heights ---\n")
print(hc_result$height)

# Get the cluster assignments for each network
clusters_k3 <- cutree(hc_result, k = 3)
cat("\n--- Cluster Assignments (k=3) ---\n")
print(clusters_k3)

# If the issue persists, consider generating random numbers in a way that
# guarantees more diversity in your attribute sets to force clearer clusters.
# For example, by ensuring a certain percentage of overlap/difference.
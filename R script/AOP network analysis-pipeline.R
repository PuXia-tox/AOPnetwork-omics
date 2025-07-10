# Load packages
library(aop)
library(graph)
library(magrittr)

# Define directories
db_dir <- "/rds/projects/x/xiap-xia-transcriptomics/AOP network/Database/"
results_dir <- "/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results"
perturbed_dir <- file.path(results_dir, "perturbed_pathways")
qaop_output_dir <- file.path(results_dir, "qAOP")

# Load scripts and data
setwd(db_dir)
source("AOP network analysis.R")
source("qAOPanalysis.R")
source("qAOPstatistic.R")

message("Loading AOP path array...")
all_path_array <- readRDS("all MIE to AO paths-updated 2020.rds")  # load once

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


# Adjust the range for batch processing if needed
for (mmm in 29:32) {
  sample_name <- sample_files[mmm]
  sample_path <- file.path(perturbed_dir, sample_name)
  
  message(paste("Processing:", sample_name))
  data <- read.csv(sample_path)[, c(1, 3)]
  
  # Run qAOP analysis
  cmke <- qAOPanalysis(data = data, all_path_array = all_path_array)
  saveRDS(cmke, file = file.path('/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results/cmKE/',
                                 paste0("AOP_layout_coords",strsplit(sample_name, '.csv')[[1]],'.rds')))
  #AopNetworkPlot(data, cmke, sample_name)
  all_aop <- qAOPstatistic(cmke)
  
  #Filter qMIE and qAO results
  qmie <- all_aop[all_aop[, 4] == "yes" & all_aop[, 5] == "yes", ]
  qao  <- all_aop[all_aop[, 4] == "yes" & all_aop[, 6] == "yes", ]
  
  # Save output
  write.csv(all_aop, file.path(qaop_output_dir, paste0("qaop-", sample_name)), row.names = FALSE)
  write.csv(qmie,   file.path(qaop_output_dir, paste0("qmie-", sample_name)), row.names = FALSE)
  write.csv(qao,    file.path(qaop_output_dir, paste0("qao-", sample_name)), row.names = FALSE)
  
  # Clean up memory
  gc()
}
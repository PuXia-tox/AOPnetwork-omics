# Load packages
library(aop)
library(graph)
library(magrittr)

# Define directories
db_dir <- "/rds/projects/x/xiap-xia-transcriptomics/AOP network/Database/"
output_path <- '/rds/projects/x/xiap-xia-transcriptomics/AOP network/Cross_Validation/Results/'
perturbed_dir <- '/rds/projects/x/xiap-xia-transcriptomics/AOP network/Cross_Validation/Pathway/'

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
  sample_name <- "TBT"
  sample_path <- file.path(perturbed_dir, sample_name)
  # then run through the '0-Figure-updated' Script from begin to end
  
  message(paste("Processing:", sample_name))
  
  # Run qAOP analysis
  # cmke <- qAOPanalysis(data = data, all_path_array = all_path_array)
  #saveRDS(cmke, file = file.path('/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results/cmKE/',
  #                               paste0("AOP_layout_coords",strsplit(sample_name, '.csv')[[1]],'.rds')))
  #AopNetworkPlot(data, cmke, sample_name)
  # gc()

  data <- read.csv(sample_path)[, c(1, 3)]
  data[,1] <- path.db[[2]][match(data[,1],path.db[[2]][,2]),1]
  data <- data[data[,1] %in% colnames(x),]
  my_mat <- matrix(ncol=ncol(x),nrow=1,6)
  colnames(my_mat) <- colnames(x)
  for(i in 1:nrow(data)){
    my_mat[,which(colnames(my_mat) == data[i,1])] <- log10(data[i,2])-3
  }

  predictions <- predict(plsMods[[1]], my_mat)
  nComps <- plsMods[[1]]$ncomp.selected
  predicted_value <- 10^(predictions$y.pred[1, nComps, 1])*1000
  print(predicted_value)

  
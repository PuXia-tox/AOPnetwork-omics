# load packages
library(aop)
library(graph)
library(magrittr)

setwd('C:/Users/xiapu/OneDrive/desktop/AOP network files')
source("AOP network analysis.R")
source("qAOPanalysis.R")
source("qAOPstatistic.R")

# linear aop database
all_path_array<-readRDS('all MIE to AO paths-updated 2020.rds') # take ~10s; just load this once

# load aop nodes and annotation information
steatosis_aop <- convert_cytoscape_to_aop('aopwiki-KER_2020_09_21.cyjs') 
# nodes in aop
nodes <-  convert_aop_to_graph(steatosis_aop) %>%
  as(., "matrix") %>%
  colnames(.) %>%
  sapply(., function(x){getAOPNodeName(steatosis_aop, x)}) %>%
  gsub("[(.*)]", "", .) %>%
  gsub(",", "", .) %>%
  cbind(matrix(.),names(.)) %>%
  toupper  %>%
  .[,c(1,3)]

# load data
setwd('C:/Users/xiapu/OneDrive/desktop/WPP/pathDB')
sample_names <- list.files(getwd())

for(mmm in 30){
  
  setwd('C:/Users/xiapu/OneDrive/desktop/WPP/pathDB')
  data <- sample_names[mmm] %>% read.csv(.) %>% .[,c(1,3)]

  # run qAOP analysis functions
  AopNetworkPlot(data)
  ppp <- qAOPanalysis(data = data,all_path_array = all_path_array) # this step may take 5~10 mins
  all_aop <- qAOPstatistic(ppp)
  
  # output results
  qmie <- all_aop[all_aop[,4]=='yes'&all_aop[,5]=='yes',]
  qao <- all_aop[all_aop[,4]=='yes'&all_aop[,6]=='yes',]
  
  write.csv(all_aop, paste0('C:/Users/xiapu/OneDrive/desktop/WPP/qAOP/','qaop-',sample_names[mmm],'.csv'),row.names = F)
  write.csv(qmie, paste0('C:/Users/xiapu/OneDrive/desktop/WPP/qAOP/','qmie-',sample_names[mmm],'.csv'),row.names = F)
  write.csv(qao, paste0('C:/Users/xiapu/OneDrive/desktop/WPP/qAOP/','qao-',sample_names[mmm],'.csv'),row.names = F)
  
  # release memory
  gc()
  
}




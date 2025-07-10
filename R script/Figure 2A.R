# load packages
library(aop)
library(graph)
library(magrittr)

setwd('C:/Users/xiapu/OneDrive/desktop/AOP network files')
source("AOP network analysis.R")
source("qAOPanalysis.R")
source("qAOPstatistic.R")

# linear aop database
#all_path_array<-readRDS('all MIE to AO paths-updated 2020.rds') # take ~10s; just load this once

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

mmm <- 39
  
setwd('C:/Users/xiapu/OneDrive/desktop/WPP/pathDB')
data <- sample_names[mmm] %>% read.csv(.) %>% .[,c(1,3)]
  
# run qAOP analysis functions

require(readxl)
require(igraph)
require(aop)
require(graph)

### loading AOP database
#steatosis_json_file <- system.file("extdata", "AOP wiki-2020 version.cyjs", package = "aop")
# AOP wiki-2020 version contain 1173 nodes
# load aop nodes and annotation information
setwd('C:/Users/xiapu/OneDrive/desktop/AOP network files')
steatosis_aop <- convert_cytoscape_to_aop('AOP-wiki.cyjs') 
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

aop_graph_u <- convert_aop_to_graph(steatosis_aop) %>%
  as(., "matrix")  %>%
  graph.adjacency(., mode="directed")

aop.seed <- read.csv('aop network seed.csv')[,2]
.Random.seed <- aop.seed

#################################################################################
### Step 1: Load in AOP database
#################################################################################
Wiki_MIE<-as.matrix(read_excel('Wiki-MIE.xlsx', sheet=1, na='NA'))
non_MIE<-as.matrix(read_excel('non-MIE.xlsx', sheet=1, na='NA'))
Wiki_AO<-as.matrix(read_excel('Wiki-AO.xlsx', sheet=1, na='NA'))
match2AOPwiki<-as.matrix(read_excel('match to AOP wiki KE-zebrafish.xlsx', sheet=1, na='NA'))
AOP<-suppressMessages(as.matrix(read_excel('AOP_Ontology_Assignments.033017.xlsx', sheet=1, na='NA')))
MIE<-AOP[grep('GO',AOP[,8]),]

#################################################################################
### Step 3: Visualization of selected AOP networks
#################################################################################

# loci of KE
path.db <- readRDS('RZT_AOP_KE_Pathway.rds')
wp <- path.db[[2]]

data.new <- data[data[,1]%in%wp[,2],]
for(i in 1:nrow(data.new)){
  s <- unique(wp[wp[,2]%in%data.new[i,1],1])
  if(length(s)>1){
    data.new[i,1] <- s[1]
    s2 <- cbind(s[-1],rep(data.new[i,2],length(s)-1))
    colnames(s2) <- colnames(data.new)
    data.new <- rbind(data.new,s2)
  } else {
    data.new[i,1] <- unique(wp[wp[,2]%in%data.new[i,1],1])
  }
}
data.new <- data.new[order(data.new[,2])[nrow(data.new):1],]


for(i in 1:nrow(data.new)){
  s <- nodes[which(nodes[,1]%in%data.new[i,1]),2]
  if(length(s)>1){
    data.new[i,1] <- s[1]
    s2 <- cbind(s[-1],rep(data.new[i,2],length(s)-1))
    colnames(s2) <- colnames(data.new)
    data.new <- rbind(data.new,s2)
    #print(s)
  }
  data.new[i,1] <- s[1]
}

ke_nodes <- as.character(data.new[,1]) # matched ke_nodes
names(ke_nodes) <- data.new[,2] # matched ke_nodes' pod values
ke.l <-which(nodes[,2]%in%ke_nodes) # ke.l = key events locations in AOP network

if(length(ke.l)<=3){
  stop("too few key events matched by your omics data; not applicable for this analysis")
}
ke.nodes <-nodes[ke.l,2] # nodes of matched KEs for further permutation analysis

# loci of MIE (this is the MIEs originally in AOP network)
s1 <- unique(as.matrix(Wiki_MIE[,2]))
s1 <- s1[!(s1%in%non_MIE[,4])] # remove non-MIE
mie <- gsub("[(.*)]", "", s1)
mie <- gsub(",", "", mie)
mie <- toupper(mie)

mie.l <-which(nodes[,1]%in%mie)
mie.nodes <-nodes[mie.l,2] # nodes of MIEs for further permutation analysis

# loci of AO (this is the AOs originally in AOP network)
s2 <- unique(as.matrix(Wiki_AO[,4:7]))
s2 <- unlist(apply(s2,2,function(x) x[-which(is.na(x))]))
s2 <- unique(unlist(s2))
aop <- gsub("[(.*)]", "",s2)
aop <- gsub(",", "", aop)
aop <- toupper(aop)

aop.l <-which(nodes[,1]%in%aop)
aop.nodes <-nodes[aop.l,2] # nodes of MIEs for further permutation analysis


# find and set the nodes to be enlarged
xx<-c()
for(i in 1:length(unclass(V(toPlot)))){
  xx[i]<-getAOPNodeName(steatosis_aop, names(unclass(V(toPlot))[i]))
}
yy <- cbind(names(unclass(V(toPlot))),xx)
ahr <- grep('AHR/ARNT',yy[,2])
node.size<-setNames(rep(7,length(names(unclass(V(toPlot))))),names(unclass(V(toPlot))))
node.size[ahr] <- 20

## Plotting primary overview of AOP network
set.seed(52)
toPlot <- aop_graph_u
toPlot.layout<-layout.fruchterman.reingold(toPlot)
V(toPlot)$plotX<-toPlot.layout[,1]
V(toPlot)$plotY<-toPlot.layout[,2]
plotLay<-cbind(V(toPlot)$plotX,V(toPlot)$plotY)

vSubCol <- rep(rgb(1,1,1,alpha=0),vcount(toPlot))
vSubCol[ke.l]<-'lightblue'
vSubCol[ahr]<-'orange'

plot(toPlot,
     vertex.label.dist=0,vertex.label.cex=0.5,
     vertex.color=vSubCol,
     vertex.frame.color=vSubCol,
     vertex.label.family="Helvetica",
     vertex.shape='circle',
     edge.arrow.size=0,edge.width=0,
     vertex.label=NA,
     vertex.size=node.size,
     #vertex.label=V(toPlot)$vertex.label,
     vertex.label.color='black',
     layout=plotLay)

vSubCol <- rep('white',vcount(toPlot))
vSubCol[mie.l]<-'green'
vSubCol[aop.l]<-'red'


node.size<-setNames(rep(4,length(names(unclass(V(toPlot))))),names(unclass(V(toPlot))))
node.size[ahr] <- 15

par(mar=c(0,0,0,0))
plot(toPlot,
     vertex.label.dist=0,vertex.label.cex=0.5,
     vertex.color=vSubCol,vertex.frame.color='black',edge.color='black',
     vertex.label.family="Helvetica",vertex.shape='circle',
     edge.arrow.size=0.1,edge.width=0.5,
     vertex.label=NA,
     vertex.size=node.size,
     vertex.label.color='black',layout=plotLay,add=T)








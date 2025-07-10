# load packges
library(aop)
library(graph)
# load aop nodes and annotation information
steatosis_json_file <- system.file("extdata", "aopwiki-KER_2020_09_21.cyjs", package = "aop")
steatosis_aop <- convert_cytoscape_to_aop(steatosis_json_file)
steatosis_graph <- convert_aop_to_graph(steatosis_aop)
aop_adjm <- as(steatosis_graph, "matrix")
node.name <- function(x){getAOPNodeName(steatosis_aop, x)}
nodes <-sapply(colnames(aop_adjm), node.name)
nodes <- gsub("[(.*)]", "", nodes)
nodes <- gsub(",", "", nodes)
nodes <- cbind(matrix(nodes),names(nodes))
nodes[,1] <- toupper(nodes[,1])

#################################################################################

setwd('/Users/sharp/Desktop/WPP/qAOP/')
files <- list.files(getwd())
files <- files[grep('qaop',files)]
sample_names <- sapply(strsplit(files,'.csv'),function(x)x[1])
sample_names <- sapply(strsplit(sample_names,'qaop-'),function(x)x[2])

DEG_list<-c()
for(i in 1:length(files)){
  DEG_list[[i]] <-read.csv(files[i])
  dp<-duplicated(DEG_list[[i]][,2])
  if(sum(dp)){
    DEG_list[[i]] <- DEG_list[[i]][-which(dp==T),]
  }
  DEG_list[[i]] <- DEG_list[[i]][order(DEG_list[[i]][,1]),]
}
names(DEG_list)<-sample_names

all.ke <- sort(unique(unlist(sapply(DEG_list,function(x) x[,1]))))
CMap <- matrix(nrow=length(all.ke),ncol=length(sample_names))
rownames(CMap) <- all.ke
colnames(CMap) <- sample_names

for(i in 1:length(DEG_list)){
  l <- which(rownames(CMap)%in%DEG_list[[i]][,1])
  CMap[l,i] <- DEG_list[[i]][,7]
}

aopData <- CMap
chem.name <- c('25#','60#','TCDD','BPA','P8')
mmm <- 5
######## Identify time-dependent laop
CMap <- aopData[,grep(chem.name[mmm],colnames(aopData))]
# load linear aop database
setwd('/Users/sharp/Desktop/Genotox')
all_path_array<-readRDS('all MIE to AO paths-updated 2020.rds') # take ~10s; just load this once


KEmatch<-function(x,CMap){
  y1<-which(x[1]%in%rownames(CMap)[!is.na(CMap[,1])])
  y2<-which(x%in%rownames(CMap)[!is.na(CMap[,2])])
  y3<-which(x%in%rownames(CMap)[!is.na(CMap[,3])])
  y4<-which(x%in%rownames(CMap)[!is.na(CMap[,4])])
  if(length(y1)&length(y2)&length(y3)&length(y4)){
    if(min(y1)<max(y4)){
      y2<-y2[which(y2>min(y1))]
      y3<-y3[which(y3<max(y4))]
      if(min(y2)<max(y3)){
        return('True')
      } else {return('False')}
    } else {return('False')}
  } else {return('False')}
}

library(pbapply) 
op <- pboptions(type = "timer")
pboptions(type = "txt")
mieChem<-pbsapply(all_path_array,function(x){
  KE<-return(KEmatch(x,CMap))
  return(KE)
})

qPath<-all_path_array[which(mieChem=='True')]
mie<-names(table(sapply(qPath,function(x)x[1])))
mie<-nodes[which(nodes[,2]%in%mie),]
mie<-cbind(mie,rep('mie',nrow(mie)))

ao<-names(table(sapply(qPath,function(x)x[length(x)])))
ao<-nodes[which(nodes[,2]%in%ao),]
ao<-cbind(ao,rep('ao',nrow(ao)))

if(nrow(mie)&nrow(ao)){
setwd('/Users/sharp/Desktop/WPP/tAOP/')
write.csv(rbind(mie,ao),paste0('tAOP-MIE-',chem.name[mmm],'.csv'))
}

gc()

###

library(pbapply) 
op <- pboptions(type = "timer")
pboptions(type = "txt")
KEmatch<-function(x,CMap){
  y1<-which(x%in%rownames(CMap)[!is.na(CMap[,1])])
  y2<-which(x%in%rownames(CMap)[!is.na(CMap[,2])])
  y3<-which(x%in%rownames(CMap)[!is.na(CMap[,3])])
  y4<-which(x%in%rownames(CMap)[!is.na(CMap[,4])])
  if(length(y1)&length(y2)&length(y3)&length(y4)){
    if(min(y1)<max(y4)){
      y2<-y2[which(y2>min(y1))]
      y3<-y3[which(y3<max(y4))]
      if(min(y2)<max(y3)){
        return('True')
      } else {return('False')}
    } else {return('False')}
  } else {return('False')}
}

mieChem<-pbsapply(all_path_array,function(x){
  KE<-return(KEmatch(x,CMap))
  return(KE)
})

qPath<-all_path_array[which(mieChem=='True')]
mie<-names(table(sapply(qPath,function(x)x[1])))
mie<-nodes[which(nodes[,2]%in%mie),]
mie<-cbind(mie,rep('mie',nrow(mie)))

ao<-names(table(sapply(qPath,function(x)x[length(x)])))
ao<-nodes[which(nodes[,2]%in%ao),]
ao<-cbind(ao,rep('ao',nrow(ao)))

setwd('/Users/sharp/Desktop/WPP/tAOP/')
write.csv(rbind(mie,ao),paste0('tAOP-ALL-',chem.name[mmm],'.csv'))

gc()



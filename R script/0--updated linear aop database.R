### loading packages
library(readxl)
library(igraph)
library(aop)
library(graph)
library(rTensor)
library(NetworkDistance)
library(ComplexHeatmap)
library(circlize)
library(graphkernels)

set.seed(1)
### load aop ker (network file)
steatosis_json_file <- system.file("extdata", "aopwiki-KER_2020_09_21.cyjs", package = "aop")
steatosis_aop <- convert_cytoscape_to_aop(steatosis_json_file)
steatosis_graph <- convert_aop_to_graph(steatosis_aop)
aop_adjm <- as(steatosis_graph, "matrix")
aop_graph <- as_adjacency_matrix(graph.adjacency(aop_adjm, mode="directed"))
aop_graph <- graph.adjacency(aop_adjm, mode="directed")
ke_id <- cbind(matrix(rownames(aop_adjm)),
               apply(matrix(rownames(aop_adjm)),1,function(x) getAOPNodeName(steatosis_aop, x)))

### load aop ke file
setwd('/Users/sharp/Desktop/toSciome/AOP network')
ke <- as.matrix(read_excel('Wiki-KE-2020-9-21.xlsx', sheet=1, na='NA'))
ke <- ke[-which(duplicated(ke[,2])==1),]
ke <- ke[which(ke[,4]%in%ke_id[,2]),] # only select ke present in aop network

colnames(ke_id)[2] <- colnames(ke)[4]
ke_id <- merge(ke_id,ke[,2:4],by=colnames(ke)[4])
colnames(ke_id) <- c('description','ke_nodes','ke_id','ke')

mie <- ke_id[which(ke_id[,4]=='MolecularInitiatingEvent'),]
ao <- ke_id[which(ke_id[,4]=='AdverseOutcome'),]

match2AOPwiki<-as.matrix(read_excel('match to AOP wiki KE.xlsx', sheet=1, na='NA'))
aop_go<-as.matrix(read_excel('AOP_Ontology_Assignments.033017.xlsx', sheet=1, na='NA'))


### get ao position (pos) and nodes id (nodes) on aop network
ao_pos <- apply(matrix(ao[,2]),1,
                function(x) which(names(components(aop_graph)$membership)==x))
ao_nodes <- as.character(ao[,2])


### get mie position (pos) and nodes id (nodes) on aop network
mie_pos <- apply(matrix(mie[,2]),1,
                function(x) which(names(components(aop_graph)$membership)==x))
mie_nodes <- as.character(mie[,2])


### get all mie to ao paths
path <- array(list(),length(ao_nodes))
for(k in 1:length(ao_nodes)){
  for(i in 1:length(mie_nodes)){
    a<-all_simple_paths(aop_graph, mie_nodes[i], ao_nodes[k],mode='out')
    if(length(a)){
      path[[k]]<-c(path[[k]],sapply(a,function(x) rownames(as.matrix(x))))
    }
  }
  print(k)
} 

all_path_array <- list()
for(i in 1:length(path)){
  all_path_array<-c(all_path_array,path[[i]])
  print(i)
}

setwd('/Users/sharp/Desktop/Genotox')
saveRDS(all_path_array,'all MIE to AO paths-updated 2020.rds')
all_path_array<-readRDS('all MIE to AO paths-updated 2020.rds')


### match omics to aop_go
setwd('/Users/sharp/Desktop/Genotox/Results/DEG')
files <- sort(list.files(getwd()))
sample_names <- sapply(strsplit(files,'-pod'),function(x)x[1])

omics <- list()
for (mmm in 1:length(sample_names)){
  
  chem.name <- sample_names[mmm]
  setwd('/Users/sharp/Desktop/Genotox/Results/GO')
  my_family <- read.csv(paste0(chem.name,'-go.csv'))
  rownames(my_family) <- my_family[,1]
  colnames(my_family)[1] <- 'GO'
  my_family <- merge(my_family[,c(1,3)],
                     match2AOPwiki[which(match2AOPwiki[,4]%in%my_family[,1]),c(1,4)],by='GO')
  if(length(which(is.na(my_family[,3])))){
    my_family <- my_family[-which(is.na(my_family[,3])),]}
  colnames(my_family)[3] <- colnames(ke_id)[1]
  my_family <- merge(my_family,ke_id[which(ke_id[,1]%in%my_family[,3]),])
  mke_nodes <- as.character(my_family[,4]) # matched ke_nodes
  names(mke_nodes) <- my_family[,3] # matched ke_nodes' pod values
  
  # all paths matched by > 1 omics go terms 
  pp<-sapply(all_path_array,function(x) if(length(which(x%in%mke_nodes))>1){x})
  pp<-pp[which(sapply(pp,function(x) class(x)!='NULL')==1)]
  
  # filter out sorted paths
  is.sorted = Negate(is.unsorted)
  ppp<-sapply(pp,function(x) {
    a <- mke_nodes[which(mke_nodes%in%x)]
    b <- x[which(x%in%mke_nodes)]
    if(is.sorted(as.numeric(names(a[apply(as.matrix(a),1,function(x) grep(paste0('^',x,'$'),b))])))){
      x
    }
  })
  ppp<-ppp[which(sapply(ppp,function(x) class(x)!='NULL')==1)]
  
  # show the number of kes meeting the sorting criteria
  enriched.nodes <- unique(unlist(ppp))
  print(paste0(mmm,', ',chem.name,', ',length(enriched.nodes)))
  
  aop_adjm <- as(steatosis_graph, "matrix")
  l<-which(rownames(aop_adjm)%in%enriched.nodes)
  
  if(length(l)){
    aop_graph_omics <- graph.adjacency(aop_adjm[l,l], mode="directed")
    plot(aop_graph_omics,layout=layout_with_kk,
         edge.arrow.size=0.3,vertex.size=6,vertex.label=NA,
         vertex.label.dist=0,vertex.label.cex=0.8,vertex.label.color='black',
         main=chem.name)
    aop_graph_omics <- as_adjacency_matrix(aop_graph_omics)
    omics[[mmm]] <- as.matrix(aop_graph_omics)
  } else {
    omics[[mmm]] <- 'NULL'
  }
  
}


setwd('/Users/sharp/Desktop/Genotox')
saveRDS(omics,'79 chemicals-aop_omics_matrix.rds')
omics<-readRDS('79 chemicals-aop_omics_matrix.rds')


na.sample <- which(sapply(omics,nrow)=='NULL')
if(length(na.sample)){
  na.sample <- which(sapply(omics,nrow)=='NULL')
  B<-omics[-na.sample]
} else {
  B<-omics  
}

xp.data <- list()
for(i in 1:length(B)){
  xp.data[[i]] <- graph.adjacency(B[[i]], mode="directed")
  V(xp.data[[i]])$name <- V(xp.data[[i]])
}



KEH <- CalculateWLKernel(xp.data, 20)
KEH <- CalculateKStepRandomWalkKernel(xp.data,rep(1,2))
KEH <- CalculateShortestPathKernel(xp.data)

###ComplexHeatmap
mydata<-KEH
if(length(na.sample)){
  rownames(mydata)<-sample_names[-na.sample]
  colnames(mydata)<-sample_names[-na.sample]
} else {
  rownames(mydata)<-sample_names
  colnames(mydata)<-sample_names
}

setwd('/Users/sharp/Desktop/Genotox')
MOA<-read.csv('chem_set.csv')
MOA<-MOA[order(as.character(MOA[,1])),]
MOA<-MOA[which(MOA[,1]%in%colnames(cor(mydata))),]
n<-c()
for(i in 1:nrow(MOA)){
  n<-c(n,which(MOA[,1]%in%colnames(cor(mydata))[i]))
}
MOA<-MOA[n,]

Hepato<-as.character(MOA[,8])
Hepato[which(Hepato==1)]<-'Yes'
Hepato[which(Hepato==0)]<-'No'
Hepato[which(Hepato==(-1))]<-'Potential'
Hepato<-as.factor(Hepato)

annotation_col <- data.frame(
  Hepatotoxicity = factor(Hepato),
  MOA = factor(MOA[,3]),
  Cell = factor(MOA[,2])
)

rownames(annotation_col) <- colnames(mydata)


ann_colors = list(Hepatotoxicity = c('Yes'='black',
                                     'No'='white',
                                     'Potential'='lightgrey',
                                     'na'='white'),
                  MOA = c('Metabolic'='orange',
                          'Genotox'='#66B2FF',
                          'EDC'='green',
                          'AHR'='red',
                          'SCCP'='white',
                          'Neuroactive'='purple',
                          'Unspecific'='black',
                          'na'='lightgrey'),
                  Cell = c(HepG2='red',
                           MCF7='white',
                           Zebrafish=4)
)

ha_col<-HeatmapAnnotation(Hepatotoxicity=as.character(annotation_col[,1]),
                          MOA=as.character(annotation_col[,2]),
                          Cell=as.character(annotation_col[,3]),
                          col=ann_colors[1:3],border = TRUE)

color_fun<- colorRamp2(c(-1,-0.33,-0.5,0.33,0.66,1), 
                       c('darkblue',"#377EB8", 'lightblue','pink',"#E41A1C",'darkred'))

dcols <- as.dist(1-cor(cor(mydata), method = "pearson"))

ht2<-Heatmap(cor(mydata),name = "HTT -log(pod)",
             top_annotation = ha_col,
             show_row_names=T,
             row_names_side = c("left"),
             column_names_side = c("bottom"),
             row_dend_width = unit(50, "mm"),
             column_dend_height = unit(50, "mm"),
             show_row_dend=T,
             clustering_distance_rows = dcols,
             clustering_method_rows = 'complete',
             #clustering_distance_columns = dcols,
             #clustering_method_columns = 'complete',
             column_names_gp = gpar(fontsize = 5),row_names_gp = gpar(fontsize = 5),
             col=color_fun)
ht2



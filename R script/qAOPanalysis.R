qAOPanalysis<-function(data,all_path_array){
  
  ### loading packages
  library(readxl)
  library(igraph)
  library(aop)
  library(graph)
  library(pbapply)
  
  ### load all excel files
  setwd('C:/Users/xiapu/OneDrive/desktop/AOP network files')
  ke <- as.matrix(read_excel('Wiki-KE-2020-9-21.xlsx', sheet=1, na='NA'))
  aop_go<-as.matrix(read_excel('match to AOP wiki KE-zebrafish.xlsx', sheet=1, na='NA'))
  match2AOPwiki<-suppressMessages(as.matrix(read_excel('match to AOP wiki KE.xlsx', sheet=1, na='NA')))
  match2AOPwiki[,1] <- gsub("[(.*)]", "", match2AOPwiki[,1])
  
  
  ### load aop ker (network file)
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
  
  aop_graph_u <- convert_aop_to_graph(steatosis_aop) %>%
    as(., "matrix")  %>%
    graph.adjacency(., mode="directed")
  
### match omics to aop_go
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
  
  mke_nodes <- as.character(data.new[,1]) # matched ke_nodes
  mke_nodes <- as.data.frame(cbind(mke_nodes,data.new[,2])) # matched ke_nodes' pod values
  mke_nodes[,2] <- as.numeric(mke_nodes[,2])
  
  # all paths matched by >= 3 omics-based go terms 
  op <- pboptions(type = "timer")
  pboptions(type = "txt")
  print('filter paths matched by omics > 3 ke')
  pp<-pbsapply(all_path_array,function(x) if(length(which(x%in%mke_nodes[,1]))>=3){x})
  pp<-pp[which(sapply(pp,function(x) class(x)!='NULL')==1)]
  if(length(pp)<=0){
    stop("no matched laop by your omics data")
  }
  
  # filter out sorted paths
  is.sorted = Negate(is.unsorted)
  
  QaopCal <- function(qaopDF){
    min.ke <- formatC(min(qaopDF[,2]), format = "e", digits = 2)
    mean.ke <- formatC(mean(qaopDF[,2]), format = "e", digits = 2)
    median.ke <- formatC(median(qaopDF[,2]), format = "e", digits = 2)
    stat<- paste0(min.ke,' (',length(qaopDF),') ', mean.ke,' ', median.ke)
    return(stat)
  }
  
  print('evaluation of dose-dependent paths')
  ppp<-pbsapply(pp,function(x) {
    a<-mke_nodes[mke_nodes[,1]%in%x,]
    b<-x[x%in%mke_nodes[,1]]
    y<-a[order(match(a[,1],b)),]
    if(is.sorted(y[,2])){
      stat<-QaopCal(a)
      result<-list('stat'=stat,'path'=x,'pod'=a)
      return(result)
    }
  })
  ppp<-ppp[which(sapply(ppp,function(x) class(x)!='NULL')==1)]
  return(ppp)
}


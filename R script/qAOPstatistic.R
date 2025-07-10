qAOPstatistic<-function(ppp){
  
  laop.pod <- sapply(ppp,function(x)x[[1]])
  laop.mie <- matrix(sapply(sapply(ppp,function(x) x[2]),function(x)x[1]))
  laop.pod.num <- sapply(strsplit(laop.pod,' '),function(x) x[2])
  
  pod.min <- sapply(strsplit(laop.pod,' '),function(x) x[1])
  pod.mean <- sapply(strsplit(laop.pod,' '),function(x) x[3])
  pod.median <- sapply(strsplit(laop.pod,' '),function(x) x[4])
  s <- cbind(matrix(laop.mie),pod.min,pod.mean,pod.median)
  min.min <- tapply(as.numeric(s[,2]),s[,1],min)
  mean.min <- tapply(as.numeric(s[,3]),s[,1],min)
  median.min <- tapply(as.numeric(s[,4]),s[,1],min)
  mie.pod <- cbind(names(table(s[,1])),matrix(min.min),matrix(mean.min),matrix(median.min))
  colnames(mie.pod) <- c('KE','min.aop','mean.aop','median.aop')
  
  s1 <- table(sapply(sapply(ppp,function(x) x[2]),function(x)x[1]))
  s1 <- cbind(names(s1),s1)
  s2 <- sort(table(sapply(sapply(ppp,function(x) x[2]),function(x)x[1])),decreasing = T)
  s2 <- nodes[nodes[,2]%in%names(s2),]
  colnames(s1) <- c('KE','num'); colnames(s2) <- c('name','KE')
  s <- merge(s1,s2,by='KE')
  s <- merge(s,mie.pod,by='KE')
  s <- s[order(as.numeric(s[,5])),]
  dim(s)
  head(s)
  #write.csv(s, paste0('/Users/sharp/Desktop/WPP/qAOP/',chem.name),row.names = F)
  
  
  ### all events/neurotoxicity-related KEs
  all_aop <- as.matrix(round(sort(table(unlist(sapply(ppp,function(x) x[[2]]))),decreasing = T)/length(ppp)*100,2))
  name <- c()
  for(i in 1:nrow(all_aop)){
    name[i] <- nodes[nodes[,2]%in%rownames(all_aop)[i],1]
  }
  all_aop <- cbind(name,all_aop)
  all_aop <- cbind(rownames(all_aop),as.matrix(all_aop))
  all_aop <- as.data.frame(all_aop)
  colnames(all_aop)[1:3] <- c('nodes','name','%ppp')
  
  all_aop <- cbind(all_aop,rep('n',nrow(all_aop)))
  colnames(all_aop)[4] <- c('matched')
  n <- which(all_aop[,1]%in%unique(unlist(sapply(ppp,function(x)x[[3]]))))
  if(length(n)){
    all_aop[n,4] <- 'yes'
  }
  
  mie <- unique(sapply(sapply(ppp,function(x)x[[2]]),function(x)x[1]))
  ao <- unique(sapply(sapply(ppp,function(x)x[[2]]),function(x)x[length(x)]))
  all_aop <- cbind(cbind(all_aop,rep('n',nrow(all_aop))),rep('n',nrow(all_aop)))
  n <- which(all_aop[,1]%in%mie)
  if(length(n)){
    all_aop[n,5] <- 'yes'
  }
  n <- which(all_aop[,1]%in%ao)
  if(length(n)){
    all_aop[n,6] <- 'yes'
  }
  colnames(all_aop)[5:6] <- c('mie','ao')
  
  pod<-sapply(strsplit(sapply(ppp,function(z)z[[1]]),' '),function(x)x[3])
  aop_pod <- apply(all_aop,1,function(x) {
    y<-x[1]
    k<-which(sapply(sapply(ppp,function(z)z[[3]])[seq(1,2*length(ppp),2)],function(z) sum(z%in%y))>0)
    if(length(k)){
      return(min(pod[k]))
    } else {
      return(NA)
    }
  })
  all_aop <- cbind(all_aop, matrix(aop_pod))
  colnames(all_aop)[7] <- c('pod')
  return(all_aop)
}
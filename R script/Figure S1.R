library(msigdbr)
setwd('C:/Users/xiapu/OneDrive/desktop/AOP network files')
Ampliseq<-read.csv('Zebrafish 1637 Reduced Genome.csv')
path.db <- readRDS('RZT_AOP_KE_Pathway.rds')
wp <- path.db[[1]]
wp <- wp[which(wp[,2]%in%Ampliseq[,4]),]

# pathway enrichment

sample_names <- list.files('C:/Users/xiapu/OneDrive/desktop/WPP/nano/DEG')

wp_list <- as.character(unique(wp[,1]))
l <- length(wp_list)
CMap <- matrix(nrow=l,ncol=length(sample_names))
rownames(CMap) <- wp_list;colnames(CMap)<-unlist(strsplit(sample_names,'.csv'))
CMap_RHT<-CMap

DEG_list<-c()
for(i in 1:length(sample_names)){
  DEG_list[[i]] <-read.csv(paste0('C:/Users/xiapu/OneDrive/desktop/WPP/nano/DEG/',sample_names[i]))
}
names(DEG_list)<-sample_names
sapply(DEG_list, function(x) nrow(x))
round(sapply(DEG_list, function(x) log10(median(x[,14]))),2)


for(mmm in 1:length(DEG_list)){
  
  DEG <- DEG_list[[mmm]][,c(1,14)]
  colnames(DEG)[1] <- colnames(Ampliseq)[3]
  DEG <- merge(DEG,Ampliseq[,c(3,4)],by='Associated.Gene.Name')
  if(sum(DEG[,2]<0)){DEG<-DEG[-which(DEG[,2]<0),]}
  DEG[,2] <- log10(DEG[,2])
  
  
  my_family <- matrix(nrow = l, ncol = 6)
  rownames(my_family) <- wp_list
  colnames(my_family) <- c('Num','POD','SD','%','Genes','KE')
  l.mat <- matrix(nrow = l, ncol = 2)
  
  for(i in 1:length(wp_list)){
    a <- which(wp[,1] == wp_list[i])
    b <- which(DEG[,3] %in% wp[a,2])
    l.mat[i,1] <- length(b)
    l.mat[i,2] <- length(a)
    x <- DEG[b,2] 
    my_family[i,1] <- length(b)
    my_family[i,2] <- 10^mean(x)
    my_family[i,3] <- sd(x)/length(b)
    my_family[i,4] <- round(length(b)/length(wp[a,2])*100,2)
    my_family[i,5] <- paste(DEG[b,1], collapse=', ')
    my_family[i,6] <- paste0(path.db[[2]][path.db[[2]][,2]%in%wp_list[i],1],collapse='; ')
    rm(b)
  }
  
  my_family <- my_family[c(which(l.mat[,1]>=3),which(l.mat[,1]==2&l.mat[,2]<15)),]
  
  if(nrow(my_family)>=3){
    setwd('C:/Users/xiapu/OneDrive/desktop/WPP/nano/pathDB')
    my_family<-my_family[rownames(my_family)%in%path.db[[2]][,2],]
    write.csv(my_family[order(as.numeric(my_family[,2])),],sample_names[mmm])
    CMap_RHT[which(rownames(CMap_RHT)%in%rownames(my_family)),mmm]<-my_family[,2]
  }
}




CMap_RHT[is.na(CMap_RHT)]<-10^6
if(length(which(rowSums(CMap_RHT==10^6)==ncol(CMap_RHT)))){
  CMap_RHT<-CMap_RHT[-which(rowSums(CMap_RHT==10^6)==ncol(CMap_RHT)),]
}
CMap_RHT<-t(apply(CMap_RHT,1,function(x)as.numeric(x)))
CMap_RHT<-log10(CMap_RHT)
colnames(CMap_RHT) <- unlist(strsplit(sample_names,'.csv'))

library(pheatmap)
pheatmap(CMap_RHT,show_rownames=F,cluster_cols=F,
         color = colorRampPalette(c("red", "white"))(100))

# PCA
pcaData <- as.data.frame(cbind(t(CMap_RHT),
                               colnames(CMap_RHT)))

pcaData[,-ncol(pcaData)] <- data.frame(apply(pcaData[,-ncol(pcaData)], 
                                             2, function(x) as.numeric(as.character(x))))

colnames(pcaData)[ncol(pcaData)] <- 'Chemicals'

library(ggfortify)
res.pca <- prcomp(pcaData[,-ncol(pcaData)], scale = TRUE, center = T)

PCA <- res.pca
top5_pc1<-names(PCA$rotation[order(abs(PCA$rotation[,1]), decreasing=TRUE),1])[1:5]
tolower(unique(path.db[[2]][path.db[[2]][,2]%in%top5_pc1,1]))
top5_pc2<-names(PCA$rotation[order(abs(PCA$rotation[,2]), decreasing=TRUE),2])[1:5]
tolower(unique(path.db[[2]][path.db[[2]][,2]%in%top5_pc2,1]))

top5_pc1<-names(PCA$rotation[order(abs(PCA$rotation[,1]), decreasing=TRUE),1])[1:5]
top5_pc2<-names(PCA$rotation[order(abs(PCA$rotation[,2]), decreasing=TRUE),2])[1:5]
top5<-unique(c(top5_pc1, top5_pc2))
tolower(unique(path.db[[2]][path.db[[2]][,2]%in%top5,1]))

PCA2<-PCA
PCA2$rotation<-PCA2$rotation[top5,]
PCA2$center<-PCA2$center[top5]
PCA2$scale<-PCA2$scale[top5]

autoplot(PCA2,data=pcaData,
         #frame = TRUE, frame.colour = 'Chemicals',
         shape = 16, size = 4,colour = 'Chemicals',
         loadings = TRUE, loadings.label = TRUE,
         loadings.label.size  = 2,
         loadings.colour = 'darkgrey',
         loadings.label.colour = 'darkgrey') + 
  xlim(c(-1,1)) +
  geom_text(aes(label = rownames(res.pca$x)), size = 2,
            nudge_y = rep(-0.01,length(rownames(res.pca$x)))) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA)) + 
  geom_hline(aes(yintercept=0), colour="black", linetype="dashed") +  
  geom_vline(aes(xintercept=0), colour="black", linetype="dashed")



### pathway boxplot
setwd('C:/Users/xiapu/OneDrive/desktop/WPP/nano/pathDB')
sample_names <- list.files('C:/Users/xiapu/OneDrive/desktop/WPP/nano/DEG')
chem_names <- sapply(strsplit(sample_names,'[.]'),function(x)x[1])
col <- sapply(strsplit(sample_names,'.csv'),function(x)x[1])

DEG_list<-c()
for(i in 1:length(sample_names)){
  DEG_list <-read.csv(sample_names[i])[,c(1,3)]
  DEG_list[,2] <- log10(DEG_list[,2])
  DEG_list[,1] <- rep(chem_names[i],nrow(DEG_list))
  if(i%%4==0){
    DEG_list <- cbind(DEG_list,rep(col[i],nrow(DEG_list)),rep(paste0('d',4),nrow(DEG_list))) 
  } else{
    DEG_list <- cbind(DEG_list,rep(col[i],nrow(DEG_list)),rep(paste0('d',i%%4),nrow(DEG_list)))
  }
  
  colnames(DEG_list)[1] <- 'samples'
  colnames(DEG_list)[3] <- 'col'
  colnames(DEG_list)[4] <- 'time'
  
  DEG_list <- cbind(DEG_list,rep(col[i],nrow(DEG_list)))
  
  if(i == 1){
    df <- DEG_list
  } else {
    df <- rbind(df, DEG_list)
  }
}

colnames(df)[1] <- 'samples'
colnames(df)[3] <- 'col'
colnames(df)[4] <- 'time'
df <- as.data.frame(df)

library(ggplot2)
give.n <- function(x){
  return(c(y = median(x)+0.13, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}


odd_numbers <- c(5:8,13:16)
p <- ggplot(df, aes(x=samples, y=POD),coef = 6, outlier.size=0.1)  + 
  geom_rect(data = df[odd_numbers, ], xmin = odd_numbers - 0.5, xmax = odd_numbers + 
              0.5, ymin = -Inf, ymax = Inf, fill = 'grey', alpha = 0.5) + 
  geom_boxplot() +
  stat_summary(fun.data = give.n, geom = "text", fun = median, 
               position = position_dodge(width = 0.75)) + 
  theme_classic() + coord_flip() +
  ylab("log10 bmd") + xlab("samples") + 
  scale_x_discrete(limits = sort(unique(df$samples),decreasing = T))
p


###
x<-as.data.frame(df[,c(1,2)])
y<-as.matrix(10^tapply(x$POD,x$samples,mean))
y<-as.matrix(y[sapply(strsplit(rownames(y),'-'),function(x)x[1])%in%c('NP1','NP5','NP6','NP9','NP10'),])
y[1,1]<-y[1,1]*1000000/150.71
y[2,1]<-y[2,1]*1000000/86.9368
y[3,1]<-y[3,1]*1000000/157.8743
y[4,1]<-y[4,1]*1000000/79.866
y[5,1]<-y[5,1]*1000000/81.38

poda<-c(2*1000000/150.71,
        0.1*1000000/86.9368,
        0.01*1000000/157.8743,
        0.04*1000000/79.866)
podt<-as.matrix(as.numeric(y))
pod_nano<-cbind(log10(podt),log10(poda))
rownames(pod_nano)<-c('SnO2','MnO2','Mn2O3','TiO2','ZnO')


### pathway boxplot
setwd('C:/Users/xiapu/OneDrive/desktop/WPP/pathDB')
cheminfo <- read_excel('C:/Users/xiapu/OneDrive/desktop/AOP network files/chemical info.xlsx') %>% as.data.frame()
sample_names <- list.files('C:/Users/xiapu/OneDrive/desktop/WPP/DEG')
chem_names <- sapply(strsplit(sample_names,'[.]'),function(x)x[1])
col <- sapply(strsplit(sample_names,'-'),function(x)x[1])

DEG_list<-c()
for(i in 1:length(sample_names)){
  DEG_list <-read.csv(sample_names[i])[,c(1,3)]
  DEG_list[,2] <- log10(DEG_list[,2])
  x <- unlist( strsplit(chem_names[i],'-'))
  x[1] <- cheminfo[which(cheminfo[,1]==x[1]),4]
  chem_names[i] <- paste(x[1],x[2],sep='-')
  DEG_list[,1] <- rep(chem_names[i],nrow(DEG_list))
  if(length(which(c('TCDD','BPA','CBP','MFA','CPS')%in%x[1]))){
  if(x[2]==1){
    DEG_list <- cbind(DEG_list,rep(col[i],nrow(DEG_list)),rep('32 hpf',nrow(DEG_list))) 
  } else{
    if(x[2]==2){
      DEG_list <- cbind(DEG_list,rep(col[i],nrow(DEG_list)),rep('48 hpf',nrow(DEG_list)))
    }
    else{
      if(x[2]==3){
        DEG_list <- cbind(DEG_list,rep(col[i],nrow(DEG_list)),rep('72 hpf',nrow(DEG_list)))
      }
      else{
          DEG_list <- cbind(DEG_list,rep(col[i],nrow(DEG_list)),rep('96 hpf',nrow(DEG_list)))
      }
    }
  }
  } else {
    DEG_list <- cbind(DEG_list,rep(col[i],nrow(DEG_list)),rep('32 hpf',nrow(DEG_list))) 
 }
  
  colnames(DEG_list)[1] <- 'samples'
  colnames(DEG_list)[3] <- 'col'
  colnames(DEG_list)[4] <- 'time'
  
  DEG_list <- cbind(DEG_list,rep(col[i],nrow(DEG_list)))
  
  if(i == 1){
    df <- DEG_list
  } else {
    df <- rbind(df, DEG_list)
  }
}

colnames(df)[1] <- 'samples'
colnames(df)[3] <- 'col'
colnames(df)[4] <- 'time'
df <- as.data.frame(df)
df[,1]<-factor(df[,1],levels=unique(df[,1])[c(1:2,7:14,19:20,25:27,15:18,3:6,21:24,28:35)])

library(ggplot2)
give.n <- function(x){
  return(c(y = median(x)+0.3, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}


odd_numbers <- c(1:4,9:12,17:20)
p <- ggplot(df, aes(x=samples, y=POD, fill=time),coef = 6, outlier.size=0.1)  + 
  geom_rect(data = df[odd_numbers, ], xmin = odd_numbers - 0.5, xmax = odd_numbers + 
              0.5, ymin = -Inf, ymax = Inf, fill = 'grey', alpha = 0.5) + 
  geom_boxplot() +
  stat_summary(fun.data = give.n, geom = "text", fun = median, 
               position = position_dodge(width = 0.75)) + 
  theme_classic() + coord_flip() +
  ylab("log10 pathway-level POD") + xlab("Chemicals") + 
  scale_fill_brewer(palette="Blues") + 
  scale_x_discrete(limits = sort(unique(df$samples),decreasing = T))
p



### regression analysis between PODa and PODt
x<-as.data.frame(df[,c(1,2)])
y<-as.matrix(10^tapply(x$POD,x$samples,mean))
y<-as.matrix(y[sapply(strsplit(rownames(y),'-'),function(x)x[1])%in%
                 c('25#','P8','3#','38#','41#','60#','BPA','TCDD',
                   '4#','7#','40#','34#','16#','47#','34#','21#',
                   '27#','92#','17#','BPAF','BPE','BPS'),])

poda<-c(0.1,
        10,
        2,
        10,
        10,
        0.022,
        80,
        1,
        0.73,
        0.032,
        10,
        2.5,
        1,
        0.01,
        0.001,
        0.01,
        0.01,
        200,
        0.001)

i<-1;pod<-list()
podt<-as.matrix(y[as.numeric(sapply(strsplit(rownames(y),'-'),function(x)x[2]))%in%i,])
pod[[i]]<-cbind(log10(podt),log10(poda))


pod_all<-rbind(pod[[1]],pod_nano)
col <- rep('red',nrow(pod_all))
col[c(4,7,12)] <- 'blue'
plot(pod_all,
     xlab='PODt',
     ylab='PODa',
     xlim=c(-5,6),
     ylim=c(-5,6),
     col=col, 
     pch=16,
     main=''
)
text(pod_all[,1], pod_all[,2]+0.3, labels=rownames(pod_all),cex=0.6)
tempLM<-lm(pod_all[,2]~pod_all[,1])
Intercept<-summary(tempLM)$coefficients[1,1]
Slope<-summary(tempLM)$coefficients[2,1]

#abline(tempLM)
tempLM<-lm(pod_all[,2]~pod_all[,1])
abline(tempLM, col="blue", lwd=2.5)
abline(0,1)
abline(-1,1, lty=3)
abline(1,1, lty=3)
abline(2,1, lty=2)
abline(-2,1, lty=2)

legend('topleft',bty='n',cex=0.7,
       c(paste0("p-value=", round(summary(tempLM)$coefficients[,4][2],2)),
         paste0("r2=", round(as.numeric(summary(tempLM)[8]),2))))





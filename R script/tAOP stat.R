setwd('/Users/sharp/Desktop/WPP/tAOP/')
tAOP <- list.files(getwd())
tAOP <- matrix(tAOP[grep('MIE',tAOP)])
tAOP <- apply(tAOP,1,function(x) read.csv(x))

setwd('/Users/sharp/Desktop/WPP/qAOP/')
qAOP <- list.files(getwd())
qAOP <- matrix(qAOP[c(grep('qmie-25#-1',qAOP),
                      grep('qmie-BPA-1',qAOP),
                      grep('qmie-P8-1',qAOP),
                      grep('qmie-TCDD-1',qAOP))])
qAOP <- apply(qAOP,1,function(x) read.csv(x))

stat <- matrix(nrow=5,ncol=3,0)
colnames(stat) <- c('tAOP','qAOP','overlap')
rownames(stat) <- c('25#-1','BPA-1','P8-1','TCDD-1','60#-1')
for(i in 1:4){
  x<-tAOP[[i]]
  y<-qAOP[[i]]
  stat[i,1] <- sum(x[,4]=='mie')
  stat[i,2] <- sum(y[,5]=='yes')
  stat[i,3] <- length(intersect(x[x[,4]=='mie',2],y[y[,5]=='yes',2]))
}
stat

#### Clustering analysis of all chemicals
setwd('/Users/sharp/Desktop/WPP/qAOP/')
qAOP <- list.files(getwd())
qAOP <- matrix(qAOP[c(grep('qmie-',qAOP))])
qAOP <- matrix(qAOP[c(grep('-1',qAOP))])
qAOP_name <- sapply(strsplit(sapply(strsplit(qAOP,'qmie-'),function(x)x[2]),'.csv'),function(x)x[1])
qAOP_data <- apply(qAOP,1,function(x) read.csv(x))
names(qAOP_data) <- qAOP_name

sapply(qAOP_data,nrow)



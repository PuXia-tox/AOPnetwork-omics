# load packages
library(magrittr)
library(ggrepel)
library(ggpubr)
library(readxl)
library(grid)
library(circlize)
library(ComplexHeatmap)
library(Hmisc)
library(tools)
library(mixOmics)
library(mdatools)
library(msigdbr)

setwd('/rds/projects/x/xiap-xia-transcriptomics/AOP network/Database/')
source("annotation compass.R")

dir_results <- ('/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results/')

############################  Figure 3  ############################

cheminfo <- read_excel('chemical info.xlsx') %>% as.data.frame()
y <- read.csv('PODt_reference chemicals.csv')
y <- as.matrix(y[sapply(strsplit(y[,1],'-'),function(x)x[1])%in%
                   c('25#','P8','60#','BPA','TCDD'),])
poda <- c(2,2.5,0.01,200,0.001)
time<-c('32hpf', '48hpf', '72hpf', '96hpf')
col <- rep('Neurotoxicants',5)
col[2] <- 'Non-neurotoxicants'
myplot <- list()

for(i in 1:4){
  
  podt <- as.matrix(y[as.numeric(sapply(strsplit(y[,1],'-'),function(x)x[2]))%in%i,])
  pod <- cbind(log10(as.numeric(podt[,2])),log10(poda),col,podt[,1]) %>% as.data.frame
  colnames(pod) <- c('PODt','PODa','Chemical Class','chemicals')
  pod[,4] <- pod[,4] %>% strsplit(.,'-') %>% sapply(.,function(x)x[1])
  pod[,4] <- cheminfo[match(pod[,4], cheminfo[,1]),4]
  pod[,1] <- as.numeric(pod[,1])
  pod[,2] <- as.numeric(pod[,2])
  tempLM<-lm(pod[,2]~pod[,1])
  p <- paste0('= ', round(summary(tempLM)$coefficients[,4][2],3),', ')
  r <- paste0(' = ', round(as.numeric(summary(tempLM)[8]),3))
  stat <- as.expression(bquote(italic('P')~.(p)~"R"^2~.(r)))
  class <- 'Chemical Class'
  
  myplot[[i]] <- ggplot(data=pod,aes(PODt, PODa, label=chemicals, color=.data[[class]])) +
    geom_abline(intercept = 0, slope = 1, lwd = 0.7,col = 'darkgrey') +
    geom_abline(intercept = 1, slope = 1, lwd = 0.7, col = 'darkgrey', lty = 2) + 
    geom_abline(intercept = -1, slope = 1, lwd = 0.7, col = 'darkgrey', lty = 2) + 
    geom_abline(intercept =  2, slope = 1, lwd = 0.7, col = 'black', lty = 2) + 
    geom_abline(intercept = -2, slope = 1, lwd = 0.7, col = 'black', lty = 2) +
    geom_abline(slope = coef(tempLM)[2], 
                intercept = coef(tempLM)[1], 
                lwd = 1.5, col = 'blue', lty = 1) +
    geom_point(aes(color = .data[[class]], fill = .data[[class]]),
               size = 5, alpha = 0.5,shape = 21,
               show.legend = FALSE) + 
    geom_text_repel(
      mapping=aes(x=PODt, y=PODa,label=chemicals),
      size=5,
      min.segment.length = 0,
      max.overlaps = 1000,
      point.padding = 0.1, 
      box.padding = unit(0.5, "lines"),
      show.legend = FALSE
    ) +
    xlim(-5, 3) + 
    ylim(-5, 3) +
    ggtitle(time[i]) +
    annotation_compass(stat,'NW',gp=gpar(fontsize=16)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),text = element_text(size = 20),
          axis.line = element_line(linewidth = 1, colour = "black"),
          axis.text = element_text(color="black"),
          axis.ticks = element_line(linewidth = 1, color="black"),
          axis.ticks.length.y = unit(.25, "cm"),
          axis.ticks.length.x = unit(.25, "cm"),
          plot.title = element_text(size=22,hjust = 0.5,face = "bold"),
          plot.margin = margin(1,1,1,1, "cm"))
}

p <- ggarrange(myplot[[1]],myplot[[2]],myplot[[3]],myplot[[4]],
               labels = c("(A)", "(B)", "(C)", "(D)"),
               font.label = list(size = 26),
               ncol = 2, nrow = 2)

ggsave(paste0(dir_results,"Figure 3.pdf"),p,width = 12, height = 12)


############################  Figure 4  ############################ 
# load dataSet
data.nano <- read.csv('nano_PoD.csv')[,-1]
data.nano[,3] <- cheminfo[match(data.nano[,4],cheminfo[,4]),5]
data.nano <- cbind(data.nano,cheminfo[,6])
data.nano[,3] <- factor(data.nano[,3], levels = c("Unknown","Yes", "No"))
data.nano <- data.nano[-11,]
colnames(data.nano)[3] <- 'Neurotoxic'
colnames(data.nano)[5] <- 'Chemical Class'
Neuro <- 'Neurotoxic'
Class <- 'Chemical Class'
tempLM<-lm(data.nano[,2]~data.nano[,1])
p <- paste0('= ', round(summary(tempLM)$coefficients[,4][2],5),', ')
r <- paste0(' = ', round(as.numeric(summary(tempLM)[8]),3))
stat <- as.expression(bquote(italic('P')~.(p)~"R"^2~.(r)))

p <- ggplot(data=data.nano,aes(PODt, PODa, label=chemicals, color=.data[[Neuro]])) +
  geom_abline(intercept = 0, slope = 1, lwd = 0.7,col = 'darkgrey') +
  geom_abline(intercept = 1, slope = 1, lwd = 0.7, col = 'darkgrey', lty = 2) + 
  geom_abline(intercept = -1, slope = 1, lwd = 0.7, col = 'darkgrey', lty = 2) + 
  geom_abline(intercept =  2, slope = 1, lwd = 0.7, col = 'black', lty = 2) + 
  geom_abline(intercept = -2, slope = 1, lwd = 0.7, col = 'black', lty = 2) + 
  geom_smooth(method='lm', lwd = 1.5, col = 'blue', lty = 1) +
  # geom_abline(slope = coef(tempLM)[2],intercept = coef(tempLM)[1],lwd = 1.5, col = 'blue', lty = 1) +
  geom_point(aes(color = .data[[Neuro]], fill = .data[[Neuro]]),
             size = 10, alpha = 0.5,shape = 21) +
  geom_text_repel(
    mapping=aes(x=PODt, y=PODa,label=chemicals),
    size=10,
    min.segment.length = 0.1,
    segment.size = 1,
    max.overlaps = 1000,
    point.padding = 0.1, 
    box.padding = unit(1, "lines"),
    show.legend = FALSE
  ) +
  xlim(-5, 3) + 
  ylim(-5, 3) +
  annotation_compass(stat,'NW',gp=gpar(fontsize=35)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),text = element_text(size = 40),
        axis.line = element_line(linewidth = 1, colour = "black"),
        axis.text = element_text(color="black"),
        axis.ticks = element_line(linewidth = 1, color="black"),
        axis.ticks.length.y = unit(.25, "cm"),
        axis.ticks.length.x = unit(.25, "cm"),
        plot.title = element_text(size=22,hjust = 0.5,face = "bold"))

ggsave(paste0(dir_results,"Figure 4.pdf"),p,width = 15, height = 12)


############################  Figure 2  ############################ 
path.db <- readRDS('RZT_AOP_KE_Pathway.rds')
hdata <- read.csv('pathDB.csv',row.names = 'X')
hdata <- hdata[,grep('.1',colnames(hdata))]*(-1)
hdata <- hdata[-which(rowSums(hdata==(-6))==ncol(hdata)),]
#hdata <- hdata[,-grep('NP',colnames(hdata))] # here we remove NP

ann_colors <- list('Neurotoxic' = c('Yes'='#EBA133',
                                    'No'='blue',
                                    'Unknown'='white'),
                   'Chemical Class'=c('Biocide'='#4DAF4A',
                                      'PAH'='black',
                                      'Drug'='white',
                                      'FR'='#E41A1C',
                                      'Industrial'='yellow',
                                      'NP'='#377EB8'))


colnames(hdata) <- gsub('X','',colnames(hdata))
ha_name <- cheminfo[match(sapply(strsplit(colnames(hdata),'[.]'),function(x)x[1]), gsub('#','',cheminfo[,1])),c(5,6)]
ha_name <- data.frame(ha_name)
colnames(ha_name) <- c('Neurotoxic','Chemical Class')

ha_col_bot <- HeatmapAnnotation(df=ha_name,
                                col = ann_colors[1:2],
                                annotation_legend_param = list(
                                  title_gp = gpar(fontsize = 14,fontface='bold'),
                                  labels_gp = gpar(fontsize = 14)
                                ),
                                border=T)

ha_col_top <- HeatmapAnnotation(PODt_boxplot = anno_boxplot(apply(hdata,2,function(x) x[x!=(-6)]),
                                                            height = unit(2,'cm'),
                                                            box_width = 0.3,
                                                            outline = F))

color_fun<-colorRamp2(c(-6,0,3,max(hdata)), c("white", "pink", "red", "purple"))

hdata.annotation <- path.db[[2]][match(rownames(hdata),path.db[[2]][,2]),1]
at <- grep('NEURO',hdata.annotation)
at <- at[-which(at %in% grep('NON-NEURONAL',hdata.annotation))]

label <- hdata.annotation[at] %>%
  gsub('N/A ', '', .) %>%
  tolower %>%
  capitalize

colnames(hdata) <- cheminfo[match(sapply(strsplit(colnames(hdata),'[.]'),function(x)x[1]), gsub('#','',cheminfo[,1])),c(4)]

pdf(file=paste0(dir_results,'Figure 2.pdf'),
    width = 10, height = 10)
draw(Heatmap(hdata,
             bottom_annotation = ha_col_bot,
             top_annotation = ha_col_top,
             col = color_fun,
             border_gp = gpar(col = "black"),
             clustering_distance_columns = 'pearson',
             column_dend_height = unit(2.5, "cm"),
             heatmap_legend_param = list(
               title = as.expression(bquote('-log10 POD'[t])),
               legend_height = unit(2, "cm"),
               title_gp = gpar(fontsize = 14,fontface='bold'),
               labels_gp = gpar(fontsize = 14),
               border = "black"
             ),
             show_row_names = F,
             show_row_dend = F,
             height = unit(10, "cm"),
             row_names_gp = gpar(fontsize = 2)) +
       rowAnnotation(label = anno_mark(at = at, labels = label)),
     heatmap_legend_side = "left",
     annotation_legend_side= 'bottom'
)
dev.off()


############################  Figure 4  ############################ 
x <- t(hdata)*(-1)
rownames(x) <- sapply(strsplit(rownames(x),'-'),function(x)x[1])
y <- as.matrix(data.nano[,c(1,2)])
rownames(y) <- data.nano[,4]
y <- y[rownames(y)%in%rownames(x),]
y <- y[order(rownames(y)),]
x <- x[rownames(x)%in%rownames(y),]
x <- x[order(rownames(x)),]
y <- cheminfo[match(rownames(x),cheminfo[,4]),5]
xc <- x
xc <- x[y!='Unknown',]
y <- y[y!='Unknown']

MyResult.splsda <- splsda(xc, y, keepX = c(20,20))
#plotIndiv(MyResult.splsda)
#plotVar(MyResult.plsda)
pdf(file=paste0(dir_results,'Figure S5.pdf'),
    width = 10, height = 10)
plotIndiv(MyResult.splsda, ind.names = T, legend= T,
          ellipse = T, star = T, title = 'PCA',
          X.label = 'PCA 1', Y.label = 'PCA 2')
dev.off()

#auc.plsda <- auroc(MyResult.splsda)

############################  PLS-DA modelling analysis  ############################
### all KEs-based modeling
x <- t(hdata)*(-1)
colnames(x) <- path.db[[2]][match(colnames(x),path.db[[2]][,2]),1]

y <- data.nano[,c(2:4)]
y <- y[y[,3]%in%rownames(x),]
y <- y[y[,2]!='No',]
x <- x[rownames(x)%in%y[,3],]
y <- y[order(y[,3]),1]
x <- x[order(rownames(x)),]

plsMods<-list()   # list to hold all models
plsMods[[1]]<-pls(x, y, center=TRUE, scale=TRUE, cv=1)
nComps<-plsMods[[1]]$ncomp.selected
#plot(plsMods[[1]])

bestMod <- 1
i <- 1
pdf(file=paste0(dir_results,'Figure 5A.pdf'),
    width = 10, height = 10)
plot(
  plsMods[[bestMod]]$cvres$y.pred[,nComps,i]~plsMods[[bestMod]]$cvres$y.ref[,i],
  xlab='Measured PODa',
  ylab='Predicted PODa',
  main=dimnames(plsMods[[bestMod]]$calres$y.pred)[[3]][i],
  col='red', pch=16, cex=1)
tempLM<-lm(plsMods[[bestMod]]$cvres$y.pred[,nComps,i]~plsMods[[bestMod]]$cvres$y.ref[,i])
abline(tempLM, col="blue", lwd=2)
abline(0,1, lty=2)
text( plsMods[[bestMod]]$cvres$y.ref[,i],
      plsMods[[bestMod]]$cvres$y.pred[,3,i]+0.15,
      names(plsMods[[bestMod]]$cvres$y.pred[,3,i]))
coef(tempLM)[2] # slope

nComps<-plsMods[[1]]$ncomp.selected
plsMods[[1]]$calres$slope[,nComps]
plsMods[[1]]$cvres$rmse[,nComps]

slope <- paste0('slope = ', round(as.numeric(coef(tempLM)[2]),3))
p <- paste0('= ', formatC(summary(tempLM)$coefficients[,4][2], format = "e", digits = 2))
p <- as.expression(bquote(italic('P')~.(p)))
r <- paste0(' = ', round(as.numeric(summary(tempLM)[8]),3))
r <- as.expression(bquote("R"^2~.(r)))
rmse <-  paste0('rmse = ',round(plsMods[[1]]$cvres$rmse[,nComps],3))
legend('bottomright',c(slope,p,r,rmse))
dev.off()


### Neuro mKEs-based modeling
# neuro-only pathways were too few to conduct pls analysis
x <- t(hdata)*(-1)
colnames(x) <- path.db[[2]][match(colnames(x),path.db[[2]][,2]),1]
at <- c(grep('NEURO',colnames(x)),
        grep('LEARN',colnames(x)),
        grep('AHR',colnames(x)))
x <- x[,at]
y <- data.nano[,c(2:4)]
y <- y[y[,3]%in%rownames(x),]
y <- y[y[,2]!='No',]
x <- x[rownames(x)%in%y[,3],]
y <- y[order(y[,3]),1]
x <- x[order(rownames(x)),]

plsMods<-list()   # list to hold all models
plsMods[[1]]<-try(pls(x, y, center=TRUE, scale=TRUE, cv=1),silent = T)

### Figure 5C -tAOP + qAOP
tAOP <- list.files('/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results/tAOP') %>% .[grep('MIE',.)] %>%
  matrix(.)  %>%
  apply(.,1,function(x) read.csv(paste0('/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results/tAOP/',x)))  %>%
  sapply(.,function(x) {
    y <- as.matrix(x)
    y[y[,4]=='mie',2]
  }) %>%
  unlist %>%
  unique

x <- t(hdata)*(-1)
colnames(x) <- path.db[[2]][match(colnames(x),path.db[[2]][,2]),1]
at.tMIE <- match(tAOP,colnames(x))
tAOP[which(is.na(at.tMIE))]
at <- at.tMIE[!is.na(at.tMIE)]
x <- x[,at]
y <- data.nano[,c(2:4)]
y <- y[y[,3]%in%rownames(x),]
y <- y[y[,2]!='No',]
x <- x[rownames(x)%in%y[,3],]
y <- y[order(y[,3]),1]
x <- x[order(rownames(x)),]

plsMods<-list()   # list to hold all models
plsMods[[1]]<-pls(x, y, center=TRUE, scale=TRUE, cv=1)
nComps<-plsMods[[1]]$ncomp.selected
#plot(plsMods[[1]])

plsMods[[1]]$calres$slope[,nComps]
plsMods[[1]]$cvres$rmse[,nComps]


bestMod <- 1
i <- 1
pdf(file=paste0(dir_results,'Figure 5B.pdf'),
    width = 10, height = 10)
plot(
  plsMods[[bestMod]]$cvres$y.pred[,3,i]~plsMods[[bestMod]]$cvres$y.ref[,i],
  xlab='Measured PODa',
  ylab='Predicted PODa',
  main=dimnames(plsMods[[bestMod]]$calres$y.pred)[[3]][i],
  col='red', pch=16, cex=1)
tempLM<-lm(plsMods[[bestMod]]$cvres$y.pred[,3,i]~plsMods[[bestMod]]$cvres$y.ref[,i])
abline(tempLM, col="blue", lwd=2)
abline(0,1, lty=2)
text( plsMods[[bestMod]]$cvres$y.ref[,i],
      plsMods[[bestMod]]$cvres$y.pred[,3,i]+0.15,
      names(plsMods[[bestMod]]$cvres$y.pred[,3,i]))

slope <- paste0('slope = ', round(as.numeric(coef(tempLM)[2]),3))
p <- paste0('= ', formatC(summary(tempLM)$coefficients[,4][2], format = "e", digits = 2))
p <- as.expression(bquote(italic('P')~.(p)))
r <- paste0(' = ', round(as.numeric(summary(tempLM)[8]),3))
r <- as.expression(bquote("R"^2~.(r)))
rmse <-  paste0('rmse = ',round(plsMods[[1]]$cvres$rmse[,nComps],3))
legend('bottomright',c(slope,p,r, rmse))
dev.off()

### Figure 5D
# load tAOP terms
tAOP <- list.files('/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results/tAOP') %>% .[grep('MIE',.)] %>%
  matrix(.)  %>%
  apply(.,1,function(x) read.csv(paste0('/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results/tAOP/',x)))  %>%
  sapply(.,function(x) {
    y <- as.matrix(x)
    y[y[,4]=='mie',2]
  }) %>%
  unlist %>%
  unique

x <- t(hdata)*(-1)
colnames(x) <- path.db[[2]][match(colnames(x),path.db[[2]][,2]),1]
at.tMIE <- match(tAOP,colnames(x))
tAOP[which(is.na(at.tMIE))]
at <- c(grep('NEURO',colnames(x)),
        grep('LEARN',colnames(x)),
        grep('AHR',colnames(x)))
at <- c(at,at.tMIE[!is.na(at.tMIE)])
x <- x[,at]
y <- data.nano[,c(2:4)]
y <- y[y[,3]%in%rownames(x),]
y <- y[y[,2]!='No',]
x <- x[rownames(x)%in%y[,3],]
y <- y[order(y[,3]),1]
x <- x[order(rownames(x)),]

plsMods<-list()   # list to hold all models
plsMods[[1]]<-pls(x, y, center=TRUE, scale=TRUE, cv=1)
nComps<-plsMods[[1]]$ncomp.selected
#plot(plsMods[[1]])

plsMods[[1]]$calres$slope[,nComps]
plsMods[[1]]$cvres$rmse[,nComps]

bestMod <- 1
i <- 1
pdf(file=paste0(dir_results,'Figure 5C.pdf'),
    width = 10, height = 10)
plot(
  plsMods[[bestMod]]$cvres$y.pred[,3,i]~plsMods[[bestMod]]$cvres$y.ref[,i],
  xlab='Measured PODa',
  ylab='Predicted PODa',
  main=dimnames(plsMods[[bestMod]]$calres$y.pred)[[3]][i],
  col='red', pch=16, cex=1)
tempLM<-lm(plsMods[[bestMod]]$cvres$y.pred[,3,i]~plsMods[[bestMod]]$cvres$y.ref[,i])
abline(tempLM, col="blue", lwd=2)
abline(0,1, lty=2)
text( plsMods[[bestMod]]$cvres$y.ref[,i],
      plsMods[[bestMod]]$cvres$y.pred[,3,i]+0.15,
      names(plsMods[[bestMod]]$cvres$y.pred[,3,i]))

slope <- paste0('slope = ', round(as.numeric(coef(tempLM)[2]),3))
p <- paste0('= ', formatC(summary(tempLM)$coefficients[,4][2], format = "e", digits = 2))
p <- as.expression(bquote(italic('P')~.(p)))
r <- paste0(' = ', round(as.numeric(summary(tempLM)[8]),3))
r <- as.expression(bquote("R"^2~.(r)))
rmse <-  paste0('rmse = ',round(plsMods[[1]]$cvres$rmse[,nComps],3))
legend('bottomright',c(slope,p,r,rmse))
dev.off()

a <- plsMods[[1]]$coeffs$values
which.max(abs(a[,,1][,1]))
sort(abs(a[,,1][,1]))
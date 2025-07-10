
input_names <- 'Androstendione'

# ===========================
# Pathway Perturbation Analysis
# ===========================

# ---- Path Setup ----
db_path <- '/rds/projects/x/xiap-xia-transcriptomics/AOP network/Database/'
deg_path <- '/rds/projects/x/xiap-xia-transcriptomics/AOP network/Cross_Validation/DEG/'
output_path <- '/rds/projects/x/xiap-xia-transcriptomics/AOP network/Cross_Validation/AOP/'
perturbed_dir <- '/rds/projects/x/xiap-xia-transcriptomics/AOP network/Cross_Validation/Pathway/'
perturbed_path <- file.path(output_path, 'perturbed_pathways')
pdf_path <- file.path(output_path, 'Figure SI-pathway HP.pdf')

# ---- Load Packages ----
library(msigdbr)
library(magrittr)
library(ggrepel)
library(readxl)
library(pheatmap)
library(org.Dr.eg.db)
library(AnnotationDbi)

# ---- Load Data ----
Ampliseq <- read.csv(file.path(db_path, 'Zebrafish 1637 Reduced Genome.csv'))
path.db <- readRDS(file.path(db_path, 'RZT_AOP_KE_Pathway.rds'))
wp <- path.db[[1]] %$% .[entrez_gene %in% Ampliseq[, 4], ]
wp_list <- unique(as.character(wp[, 1]))
l <- length(wp_list)

# ---- DEG Files ----
DEG_list <- list()
DEG_list[[1]] <- read.csv(paste0(deg_path,input_names,'.csv'))
sample_names <- input_names
names(DEG_list) <- sample_names
MW <- 286.41


# ---- Enrichment Loop ----
mmm <- 1
DEG <- DEG_list[[mmm]][, c(1, 14)]
entrez_ids <- mapIds(org.Dr.eg.db,
                     keys = DEG$gene.id,
                     column = "ENTREZID",
                     keytype = 'ENSEMBL',
                     multiVals = "first")
DEG$entrez.id <- entrez_ids
DEG <- DEG[DEG[, 2] > 0, ]
DEG[, 2] <- log10(DEG[, 2])
DEG <- DEG[,c(3,2)]
if(any(is.na(DEG[,1]))){
  DEG <- DEG[!is.na(DEG[,1]),]
}


# Per-pathway calculations
my_family <- matrix(NA, nrow = l, ncol = 6,
                    dimnames = list(wp_list, c('Num', 'POD', 'SD', '%', 'Genes', 'KE')))
l.mat <- matrix(NA, nrow = l, ncol = 2)

for (i in seq_along(wp_list)) {
  
  a <- which(wp[, 1] == wp_list[i])
  b <- which(as.character(DEG[,1]) %in% wp[a, 2])
  
  l.mat[i, ] <- c(length(b), length(a))
  x <- DEG[b, 2]
  
  my_family[i, 'Num'] <- length(b)
  my_family[i, 'POD'] <- 10^mean(x)
  my_family[i, 'SD']  <- sd(x) / length(b)
  my_family[i, '%']   <- round(100 * length(b) / length(a), 2)
  my_family[i, 'Genes'] <- paste(DEG[b, 1], collapse = ', ')
  my_family[i, 'KE'] <- paste(path.db[[2]][path.db[[2]][, 2] == wp_list[i], 1], collapse = '; ')
}

# Filter and save
valid_rows <- which(l.mat[, 1] >= 1)
my_family <- my_family[valid_rows, , drop = FALSE]

if (nrow(my_family) >= 3) {
  valid_pathways <- path.db[[2]][, 2]
  my_family <- my_family[rownames(my_family) %in% valid_pathways, , drop = FALSE]
  
  write.csv(my_family[order(as.numeric(my_family[, 'POD'])), ],
            file.path('/rds/projects/x/xiap-xia-transcriptomics/AOP network/Cross_Validation/Pathway/', sample_names[mmm]))
}




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


############################  Figure 1A  ############################

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

ggsave(paste0(dir_results,"Figure 1A.pdf"),p,width = 12, height = 12)


############################  Figure 1B  ############################ 
# load dataSet
data.nano <- read.csv('nano_PoD.csv')[,-1]
data.nano[,3] <- cheminfo[match(data.nano[,4],cheminfo[,4]),5]
data.nano <- cbind(data.nano,cheminfo[,6])
data.nano[,3] <- factor(data.nano[,3], levels = c("Unknown","Yes", "No"))
colnames(data.nano)[3] <- 'Neurotoxic'
colnames(data.nano)[5] <- 'Chemical Class'
Neuro <- 'Neurotoxic'
Class <- 'Chemical Class'
tempLM<-lm(data.nano[,2]~data.nano[,1])
p <- paste0('= ', round(summary(tempLM)$coefficients[,4][2],3),', ')
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

ggsave(paste0(dir_results,"Figure 1B-new.pdf"),p,width = 15, height = 12)


############################  Figure 2  ############################ 

files <- list.files('/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results/qAOP') %>%
  .[grep('qaop',.)] 

sample_names <- strsplit(files,'.csv') %>%
  sapply(.,function(x)x[1]) %>%
  strsplit(.,'qaop-') %>%
  sapply(.,function(x)x[2])

sample_names <- paste(cheminfo[match(sapply(strsplit(sample_names,'-'),function(x)x[1]), cheminfo[,1]),4],
                      sapply(strsplit(sample_names,'-'),function(x)x[2]),sep='-')

DEG_list<-c()
for(i in 1:length(files)){
  DEG_list[[i]] <- read.csv(paste0('/rds/projects/x/xiap-xia-transcriptomics/AOP network/Results/qAOP/',files[i]))
  dp <- duplicated(DEG_list[[i]][,2])
  if(sum(dp)){
    DEG_list[[i]] <- DEG_list[[i]][-which(dp==T),]
  }
  DEG_list[[i]] <- DEG_list[[i]][order(DEG_list[[i]][,2]),]
}
names(DEG_list)<-sample_names

all.ke <- sort(unique(unlist(sapply(DEG_list,function(x) x[,2]))))
CMap <- matrix(nrow=length(all.ke),ncol=length(sample_names))
rownames(CMap) <- all.ke
colnames(CMap) <- sample_names

for(i in 1:length(DEG_list)){
  l <- which(rownames(CMap)%in%DEG_list[[i]][,2])
  CMap[l,i] <- DEG_list[[i]][,7]
  print(length(l))
}

CMap[is.na(CMap)]<-10^6
if(length(which(rowSums(CMap==10^6)==ncol(CMap)))){
  CMap<-CMap[-which(rowSums(CMap==10^6)==ncol(CMap)),]
}


p1 <- pheatmap(log10(CMap[c(grep('NEURO',rownames(CMap)),
                            grep('LEARN',rownames(CMap))),grep('-1',colnames(CMap))]),
               show_rownames=T,fontsize_row=5,cluster_cols=T,
               #scale='column',
               color = colorRampPalette(c("red", "white"))(100))


###
hdata <- log10(CMap)[,grep('-1',colnames(CMap))]*(-1)
hdata <- hdata[-which(rowSums(hdata==(-6))==ncol(hdata)),]
#hdata <- hdata[,-(19:22)]
ann_colors <- list('Neurotoxic' = c('Yes'='#EBA133',
                                    'No'='blue',
                                    'Unknown'='white'),
                   'Chemical Class'=c('Biocide'='#4DAF4A',
                                      'PAH'='black',
                                      'Drug'='white',
                                      'FR'='#E41A1C',
                                      'Industrial'='yellow',
                                      'NP'='#377EB8'))

ha_name <- cheminfo[match(sapply(strsplit(colnames(hdata),'-'),function(x)x[1]), cheminfo[,4]),c(5,6)]
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
at <- c(grep('NEURO',rownames(hdata)),
        grep('LEARN',rownames(hdata)))
label <- rownames(hdata)[at] %>%
  gsub('N/A ', '', .) %>%
  tolower %>%
  capitalize

colnames(hdata) <- gsub('-1','',colnames(hdata))

pdf(file=paste0(dir_results,'Figure 2-new.pdf'),
    width = 10, height = 10)
draw(Heatmap(hdata,
             bottom_annotation = ha_col_bot,
             top_annotation = ha_col_top,
             col = color_fun,
             border_gp = gpar(col = "black"),
             clustering_distance_columns = 'euclidean',
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



############################  Figure 3  ############################ 
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

pdf(file=paste0(dir_results,'Figure 3-new.pdf'),
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
pdf(file=paste0(dir_results,'Figure 4.pdf'),
    width = 10, height = 10)
plotIndiv(MyResult.splsda, ind.names = T, legend= T,
          ellipse = T, star = T, title = 'PCA',
          X.label = 'PCA 1', Y.label = 'PCA 2')
dev.off()

#auc.plsda <- auroc(MyResult.splsda)



############################  Figure 5  ############################
### Figure 5A - qAOP
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

plsMods_allke <- plsMods[[1]]
x_allke<- x

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
plsMods_neuro_tmke <- plsMods[[1]]
x_neuro_tmke <- x


#################################################################################

setwd(perturbed_dir)
sample_path <- file.path(perturbed_dir, input_names)

x <- x_neuro_tmke
data <- read.csv(sample_path)[, c(1, 3)]
data[,1] <- path.db[[2]][match(data[,1],path.db[[2]][,2]),1]
data <- data[data[,1] %in% colnames(x),]
my_mat <- matrix(ncol=ncol(x),nrow=1,6)
colnames(my_mat) <- colnames(x)
for(i in 1:nrow(data)){
  my_mat[,which(colnames(my_mat) == data[i,1])] <- log10(data[i,2]/MW)
}
predictions <- predict(plsMods_neuro_tmke, my_mat)
nComps <- plsMods[[1]]$ncomp.selected
predicted_value_neuro_tmke <- 10^(predictions$y.pred[1, nComps, 1])/MW


x <- x_allke
data <- read.csv(sample_path)[, c(1, 3)]
data[,1] <- path.db[[2]][match(data[,1],path.db[[2]][,2]),1]
data <- data[data[,1] %in% colnames(x),]
my_mat <- matrix(ncol=ncol(x),nrow=1,6)
colnames(my_mat) <- colnames(x)
for(i in 1:nrow(data)){
  my_mat[,which(colnames(my_mat) == data[i,1])] <- log10(data[i,2]/MW)
}
predictions <- predict(plsMods_allke, my_mat)
nComps <- plsMods[[1]]$ncomp.selected
predicted_value_allke <- 10^(predictions$y.pred[1, nComps, 1])/MW



print(paste0('neuro+cmke model: ', round(predicted_value_neuro_tmke,4)))
print(paste0('all KE-based model: ', round(predicted_value_allke,4)))

4.34/MW

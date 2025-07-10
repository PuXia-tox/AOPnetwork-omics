# load packages
library(readxl)
library(RecordLinkage)


setwd('/rds/projects/x/xiap-xia-transcriptomics/AOP network/Database/')
# here we just analyze on the 1637 reduced zebrafish genes // The other users can also use whole geneome genes
Ampliseq <- read.csv('Zebrafish 1637 Reduced Genome.csv')

# load all zebrafish pathway database
my_list_1<-readRDS('go_list.rds')
my_list_2<-readRDS('org.Dr.eg.db compeleting msigdbr.rds')
my_list_3<-readRDS('obsolete.go.rds')
my_list_4<-readRDS('GO0023052.change.RZT.rds')
my_list<-c(my_list_1,my_list_2,my_list_3,my_list_4)
my_list<-my_list[order(names(my_list))]
l<-length(my_list)
wp_list<-names(my_list)

library(readxl)
AOP_ontology_2026 <- suppressMessages(read_excel("AO_wiki_downloads.xlsx", na = "NA", sheet ='Sheet1') %>% as.matrix())
AOP_ontology_2026[,2] <- sapply(strsplit(AOP_ontology_2026[,2], ":"), function(x) x[2])

#select DoA of well-constructed zebrafish aops
#DoA <- c('271','292','297','301','309','310','323','334','346','348','349','363','410') %>% paste0('Aop:',.)
#AOP_ontology_2026 <- AOP_ontology_2026[which(AOP_ontology_2026[,1] %in%DoA),]


AOP_keID_2026 <- suppressMessages(read_excel("AO_wiki_downloads.xlsx", na = "NA", sheet ='Sheet2') %>% as.matrix())
AOP_keID_2026 <- AOP_keID_2026[!is.na(AOP_keID_2026[,1]),1:2]
AOP_keID_2026[,1] <- gsub("^\\s+|\\s+$", "", AOP_keID_2026[,1])
merged_AOP <- merge(AOP_ontology_2026, AOP_keID_2026, by.x = "key event id", by.y = "ID", all.x = TRUE)

# review the aop-wiki-curated mKE annotations
head(merged_AOP[merged_AOP[,'process ontology id'] %in% names(tail(sort(table(merged_AOP[,'process ontology id'])),20)),])
# The results showed that the biological description of the most abundant GO-KE annotations were way too general!!!
# There, we decide to create our own mKE annotations by semantic similarity
AOP <- merged_AOP

################## Retrieve pathway ID and their full names
# C2: Wiki+KEGG+REACTOME+BIOCARTA
m_t2g <- msigdbr(species = "Danio rerio", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)
wp<-as.data.frame(m_t2g)
wp <- wp[which(wp[,2]%in%Ampliseq[,4]),]
wp <- wp[which(wp[,1]%in%names(which(table(wp[,1])>=3))),]
wp <- wp[c(grep('WP_',wp[,1]),grep('KEGG_',wp[,1]),grep('REAC',wp[,1])),]
wp[,1] <- substr(wp[,1], 4, nchar(wp))
wp[,1] <- gsub("_", " ",wp[,1])
C2 <- wp

# Gene Ontology
m_t2g <- msigdbr(species = "Danio rerio", category = "C5") %>% 
  dplyr::select(gs_name, entrez_gene)
wp<-as.data.frame(m_t2g)
wp <- wp[which(wp[,2]%in%Ampliseq[,4]),]
wp <- wp[which(wp[,1]%in%names(which(table(wp[,1])>=3))),]
wp <- wp[grep('GOBP_',wp[,1]),]
wp[,1] <- substr(wp[,1], 6, nchar(wp))
wp[,1] <- gsub("_", " ", wp[,1])
go <- wp

wp <- c(unique(C2[,1]),unique(go[,1]))


# Step 1: Get unique, clean inputs
aop_kes <- unique(na.omit(AOP$Title))         # Biological events from AOPs
aop_kes <- toupper(gsub('[(.*)]','',aop_kes))
wp_terms <- unique(wp)                           # Pathway or phenotype descriptions


### semantic similarity calculation
# solution source: https://stackoverflow.com/questions/57092479/finding-the-cosine-similarity-of-a-sentence-with-many-others-in-r
# loop visulization source: https://medium.com/human-in-a-machine-world/displaying-progress-in-a-loop-or-function-r-664796782c24
match.aop <- matrix(nrow = length(aop_kes), ncol = 3)
colnames(match.aop) <- c("Query", "Best_Match", "Cosine_Similarity")
pb <- txtProgressBar(min = 0, max = length(aop_kes), style = 3)
for (i in seq_along(aop_kes)) {
  Sys.sleep(0.1)
  
  # Combine current AOP query with all wp terms
  sv <- c(wp, aop_kes[i])
  
  # Clean and tokenize each entry
  sv_clean <- tolower(gsub("[^a-zA-Z0-9 ]", "", sv))
  svs <- strsplit(sv_clean, "\\s+")
  
  # Filter out empty or malformed tokens
  svs <- svs[sapply(svs, function(x) is.character(x) && length(x) > 0 && all(nzchar(x)))]
  
  # Explicitly assign names so stack() has matching indices
  names(svs) <- paste0("item_", seq_along(svs))
  
  # Proceed only if stackable
  if (length(svs) > 1) {
    termf <- table(stack(svs))
    idf <- log(1 / rowMeans(termf != 0))
    tfidf <- termf * idf
    last.dp <- dim(tfidf)[2]
    
    if (last.dp > 1) {
      dp <- t(tfidf[, last.dp]) %*% tfidf[, -last.dp]
      cosim <- dp / (sqrt(colSums(tfidf[, -last.dp]^2)) * sqrt(sum(tfidf[, last.dp]^2)))
      
      best_idx <- which.max(cosim)
      if (length(best_idx) > 0 && !is.na(cosim[best_idx])) {
        match.aop[i, 1] <- aop_kes[i]
        match.aop[i, 2] <- wp[best_idx]
        match.aop[i, 3] <- round(cosim[best_idx], 3)
      }
    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)
match_df <- as.data.frame(match.aop, stringsAsFactors = FALSE)
match_df <- match_df[order(match_df[,3]),]
write.csv(match_df,'match_df.csv')

### then manual checking on Cosine_Similarity < 0.3
saveRDS(path.db.list,'RZT_AOP_KE_Pathway.rds')


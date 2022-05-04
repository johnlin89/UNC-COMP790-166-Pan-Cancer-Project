.libPaths("/nas/longleaf/home/seyoun/R/x86_64-pc-linux-gnu-library/4.1")

library(NMF)
library(Rtsne)
library(dbscan)
library(mclust)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(scales)
library(egg)
library(reshape2)
library(diceR)
library(data.table)
library(dplyr)
library(factoextra)
nmfAlgorithm()

rnaseq=as.matrix(fread("/pine/scr/s/e/seyoun/03.comp790/project_comp790/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct",header=T))
names <- do.call(rbind,strsplit(colnames(rnaseq)," "))[,1]
colnames(rnaseq) <- names
meta <- as.matrix(read.table("/pine/scr/s/e/seyoun/06.Liz/01.scRNA_data/SCP542metadataMetadata.txt",header=T,sep='\t'))
c_nm <- as.matrix(read.table("/pine/scr/s/e/seyoun/06.Liz/01.scRNA_data/cellline_nm_198.txt",header=T))
subset_191 <- rnaseq[,is.element(colnames(rnaseq),c_nm)]
dim(subset_191)
row.names(rnaseq) <- rnaseq[,1]
dim(rnaseq)
exp_no_zero_num <-matrix(as.numeric(unlist(subset_191)),nrow=nrow(subset_191))
row.names(exp_no_zero_num) <- row.names(subset_191)
colnames(exp_no_zero_num) <- colnames(subset_191)
aa <- nmf_programs(exp_no_zero_num,is.log=F, rank=6, method="snmf/r", seed=1)
dim(aa$w_basis)


w_cluster1 <-aa$w_basis[order(aa$w_basis[,1],decreasing=TRUE),1]
w_cluster2 <-aa$w_basis[order(aa$w_basis[,2],decreasing=TRUE),2]
w_cluster3 <-aa$w_basis[order(aa$w_basis[,3],decreasing=TRUE),3]
w_cluster4 <-aa$w_basis[order(aa$w_basis[,4],decreasing=TRUE),4]
w_cluster5 <-aa$w_basis[order(aa$w_basis[,5],decreasing=TRUE),5]
w_cluster6 <-aa$w_basis[order(aa$w_basis[,6],decreasing=TRUE),6]

total_genes <- unique(c(names(w_cluster1[1:50]),names(w_cluster2[1:50]),names(w_cluster3[1:50]),names(w_cluster4[1:50]),names(w_cluster5[1:50]),names(w_cluster6[1:50])))
w_ma <- aa$w_basis[is.element(row.names(aa$w_basis),total_genes),]

fi_dimention_reduct <- exp_no_zero_num[is.element(row.names(exp_no_zero_num),total_genes),]
write.table(fi_dimention_reduct,"CL191_265genes.csv", quote=F, sep=",")

t_265genes <- t(fi_dimention_reduct)

hclustfunc <- function(x, method = "complete", dmeth = "euclidean") {    
  hclust(dist(x, method = dmeth), method = method)
}



#clustering
hclust_avg <- hclust(fi_dimention_reduct)
hclust_avg <- hclustfunc(t_265genes)
plot(hclust_avg)

cut_avg <- cutree(hclust_avg, k = 6)
plot(hclust_avg)
rect.hclust(hclust_avg , k = 6)
abline(h = 6, col = 'red')

#kras genes
ras_gene <- as.matrix(fread("/pine/scr/s/e/seyoun/06.Liz/ras_pathway_genes.txt",header=T))
kras_exp <-rnaseq[is.element(rnaseq[,"Description"],ras_gene),]
colnames(kras_exp) <- names
kras_191 <- kras_exp[,is.element(colnames(kras_exp),c_nm)]
row.names(kras_191) <- kras_exp[,1]
kras_num <-matrix(as.numeric(unlist(kras_191)),nrow=nrow(kras_191))
row.names(kras_num) <- row.names(kras_191)
colnames(kras_num) <- colnames(kras_191)
write.table(kras_num,"kras_227genes.csv", quote=F, sep=",")


#clustering for kras

t_kras <- t(kras_num)
#hclust_avg <- hclust(fi_dimention_reduct)
hclust_avg <- hclustfunc(t_kras)
plot(hclust_avg)

cut_avg <- cutree(hclust_avg, k = 6)
plot(hclust_avg)
rect.hclust(hclust_avg , k = 6)
abline(h = 6, col = 'red')


#kmean
distance <- get_dist(fi_dimention_reduct)
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

k2 <- kmeans(t_kras, centers = 3, nstart = 50)
fviz_cluster(k2, data = t_kras)

nmf_programs <- function(cpm, is.log=F, rank, method="snmf/r", seed=1) {
  
  if(is.log==F) CP100K_log <- log2((cpm/10) + 1) else CP100K_log <- cpm
  CP100K_log <- CP100K_log[apply(CP100K_log, 1, function(x) length(which(x > 3.5)) > ncol(CP100K_log)*0.02),]
  CP100K_log <- CP100K_log - rowMeans(CP100K_log)
  CP100K_log[CP100K_log < 0] <- 0
  
  nmf_programs <- nmf(CP100K_log, rank=rank, method=method, seed=seed)
  
  nmf_programs_scores <- list(w_basis=basis(nmf_programs), h_coef=t(coef(nmf_programs)))
  
  return(nmf_programs_scores)
}


#ras_gene <- fread("/pine/scr/s/e/seyoun/06.Liz/ras_pathway_genes.txt",header=T)

row.names(rnaseq) <- rnaseq[,2]
rnaseq_df <-rnaseq[,-c('Name', 'Description')]
#res <- nmf(rnaseq_df,6)



test <- cbind(rnaseq, total = rowSums(rnaseq_df))
#nonzero_test <- test[which(test[,1157] != '0'),]
#dim(nonzero_test)
df_new <- test[order(test[,1159],decreasing=TRUE),]
test_fi <- df_new[1:2000,1:1159]
test_fi[1517,2] <- paste("ENSG00000268942.1")
name <- rbind(as.matrix(test_fi[,'Description']))

row.names(test_fi) <- name

test_fi <-test_fi[,-c('Name', 'Description')]

write.table(test_fi,"final_test_set_2000.csv",quote=F, sep=",")


res <- nmf(test_fi,6)
V.hat <- fitted(res) 
print(V.hat) 

#test
n <- 50; counts <- c(5, 5, 8);
V <- syntheticNMF(n, counts)

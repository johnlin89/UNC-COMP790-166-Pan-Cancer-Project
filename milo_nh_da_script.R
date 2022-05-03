# This R script is intended to run on UNC's Longleaf Cluster
# https://its.unc.edu/research-computing/longleaf-cluster/
# Note that installment of packages and differences in versions may cause some
# deviations in installation/code reproducibility


# This code leverages miloR https://github.com/MarioniLab/miloR to put
# difference cell lines into neighborhoods based on sample (cell line) and
# condition (disease type)
# Differential abundance between neighborhoods is done to determine significant
# gene signatures between neighborhoods

# installation of required packages not found on UNC longleaf R version
# BiocManager::install("miloR")
# BiocManager::install("scater", force = T)


# loading required libraries
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
library(miloR)
library(tidyverse)
library(statmod)
library(stringr)

# specify file paths
figpath <- "~/class/UNC-COMP790-166-Pan-Cancer-Project/data/milo_figures" # output for figures
outpath <- "~/class/UNC-COMP790-166-Pan-Cancer-Project/data/" # output for data
inpath <- "~/class/UNC-COMP790-166-Pan-Cancer-Project/data/" # input for data


# alternative files
# mat <- read.delim(paste(inpath, "CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct", sep=""), sep=",")
# mat2 <- mat
# row.names(mat2) <- paste(mat2$Name, mat2$Description, sep="_")
# mat2 <- mat2[ , !(names(mat2) %in% c("Name", "Description"))]

mat <- read.delim(paste(inpath, "CL191_265genes.csv", sep=""), sep=",")
mat2 <- mat

# manipulate single cell exp to cell line exp
cell.line.exp <- SingleCellExperiment(assays = list(counts = mat2))
# get log norm counts
cell.line.exp <- scuttle::logNormCounts(cell.line.exp)
# dim reduction
cell.line.exp <- runPCA(cell.line.exp, ncomponents=30)
cell.line.exp <- runUMAP(cell.line.exp)
# each point is different cell line
plotUMAP(cell.line.exp)

# establish exp design
df <- data.frame(colnames(mat2))
colnames(df) <- "cell_line"
# samples are the different cell lines
df$sample <- sapply(strsplit(df$ cell_line, "_"), "[", 1)
# conditions are the different cancer types
df$condition <- sub('^[A-Z0-9]+_', '', df$cell_line)
# add exp design metadata to cle object
colData(cell.line.exp) <- DataFrame(df)
jpeg(paste(outpath, "/milo_figures/milo-fig1-umap.jpg", sep=""))
plotUMAP(cell.line.exp, colour_by="condition")
dev.off()

# neighborhoods - can experiment with prop, k and d
k <- 50
d <- 20
prop <- 0.3
# k <- 40
# d <- 20
# prop <- 0.3
cell.line.exp <- buildGraph(cell.line.exp, k = k, d = d)
cell.line.exp <- makeNhoods(cell.line.exp, prop = prop, k = k, d=d, refined = TRUE)
# distribution of nh sizes
jpeg(paste(outpath, '/milo_figures/milo-fig3-nh-hist.jpg', sep=""))
plotNhoodSizeHist(cell.line.exp)
dev.off()


# count cells in each nh
milo_pan_traj <-countCells(cell.line.exp, meta.data = data.frame(colData(cell.line.exp)), sample="sample")
nhoodCounts(milo_pan_traj)
colnames(milo_pan_traj@nhoodCounts) <- df$sample

# scree plot
percent.var <- attr(reducedDim(milo_pan_traj), "percentVar")
jpeg(paste(outpath, '/milo_figures/milo-fig2-scree.jpg', sep=""))
plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")
dev.off()

# get design df to input into DA analysis
design <- df
design <- distinct(design[c("sample", "condition")])
rownames(design) = design$sample

# DA analysis
# get distances between nhs
milo_pan_traj <- calcNhoodDistance(milo_pan_traj, d=30)
da_results <- testNhoods(milo_pan_traj, design=~condition, design.df = design)
# troubleshooting testNhoods function
# x <- milo_pan_traj
# min.mean <- 0
# design <- ~condition
# data <- df
# model <- model.matrix(object=design, data = df)
# rownames(model) <- rownames(df)
# keep.nh <- rowMeans(nhoodCounts(x)) >= min.mean
# keep.samps <- colnames(nhoodCounts(x)[keep.nh, ]) # keep.samps is null....bec there are no colnames of nhood counts
# dge <- DGEList(counts = nhoodCounts(x)[keep.nh, keep.samps], lib.size = colSums(nhoodCounts(x)[keep.nh, keep.samps]))
# dge <- estimateDisp(dge, model)
# fit <- glmQLFit(dge, model, robust = robust)
da_results %>%
  arrange(- SpatialFDR) %>%
  head()

# nodes are nhs and edges indicate how many cell lines are shared
milo_pan_traj <- buildNhoodGraph(milo_pan_traj)
jpeg('~/class/UNC-COMP790-166-Pan-Cancer-Project/data/milo_figures/milo-fig4-umap-nhs.jpg')
plotUMAP(milo_pan_traj) + plotNhoodGraphDA(milo_pan_traj, da_results, alpha=0.05) +
  plot_layout(guides="collect")
dev.off()

# look for gene signatures
# da_results$NhoodGroup <- as.numeric(da_results$SpatialFDR < 0.1 & da_results$logFC < 0)
# group nhs into groups
da_results <- groupNhoods(milo_pan_traj, da_results, da.fdr = 0.3)
# head(da_results)
jpeg('~/class/UNC-COMP790-166-Pan-Cancer-Project/data/milo_figures/milo-fig5-umap-nh-groups.jpg')
plotNhoodGroups(milo_pan_traj, da_results) 
dev.off()
# plotNhoodGroups(milo_pan_traj, da_results, layout="umap") 
## Exclude zero counts genes
# keep.rows <- rowSums(logcounts(milo_pan_traj)) != 0
# milo_pan_traj <- milo_pan_traj[keep.rows, ]
colnames(milo_pan_traj) <- colnames(mat2)
nhood_markers <- findNhoodGroupMarkers(milo_pan_traj, da_results, 
                                       aggregate.samples = TRUE, sample_col = "sample")
nhood_markers




# rows are cell lines and columns are neighborhoods
# https://github.com/MarioniLab/miloR/issues/163
nh <- nhoods(milo_pan_traj)
nh2 <- as.data.frame(as.matrix((nh)), row.names = colnames(mat2))
colnames(nh2) <- 1:ncol(nh2)
# default col names are index cells for each neighborhood ie nhoodIndex(milo_pan_traj)
nh2$rowsum <- rowSums(nh2)
nh3 <- nh2[base::order(nh2$rowsum),]
nh4 <- nh3[ , !(names(nh3) %in% c("rowsum"))]
# write.csv(nh4, file="/pine/scr/j/c/jclin/milo_nh_output.csv")

### get cell --> nh --> nh groups --> differential genes
colnames(nh4) <- paste("nh", colnames(nh4), sep="")
for(i in 1:12) {
  coli <- paste("nh", i, sep="")
  if(i == 1) {
    first <- nh4[nh4[coli]==1,][i]
    final <- first
    colnames(final) <- c("nh")
  } else {
    temp <- nh4[nh4[coli]==1,][i]
    colnames(temp) <- c("nh")
    temp$nh <- i
    final <- base::rbind(final, temp, make.row.names=T)
  }
}

nh_to_group <- data.frame(c(1:nrow(da_results)), da_results$NhoodGroup)
colnames(nh_to_group) = c("nh", "nh.group")

gr1_markers <- as.data.frame(nhood_markers[c("logFC_1", "adj.P.Val_1")] )
gr2_markers <- as.data.frame(nhood_markers[c("logFC_2", "adj.P.Val_2")] )
# gene signature is where logFC > 0
gr1_markers <- gr1_markers[gr1_markers$logFC_1>0,]
gr2_markers <- gr2_markers[gr2_markers$logFC_2>0,]
# create data frames
gr1_markers <- data.frame(row.names(gr1_markers), 1)
gr2_markers <- data.frame(row.names(gr2_markers), 2)
colnames(gr1_markers) <- c("gene", "nh.group")
colnames(gr2_markers) <- c("gene", "nh.group")
markers <- rbind(gr1_markers, gr2_markers)
markers$nh.group <- as.character(markers$nh.group)

final2 <- merge(final, nh_to_group, by="nh")
final2$cell_line <- row.names(final)
final2 <- final2[c("cell_line", "nh", "nh.group")]
final2$nh <- as.character(final2$nh)

write.csv(final2, file="/pine/scr/j/c/jclin/milo_nh_output.csv", row.names=F)
write.csv(markers, file="/pine/scr/j/c/jclin/milo_nh_to_gene_output.csv", row.names=F)

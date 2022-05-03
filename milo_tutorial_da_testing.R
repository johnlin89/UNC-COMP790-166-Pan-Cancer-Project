library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
library(miloR)


data("sim_trajectory", package = "miloR")

## Extract SingleCellExperiment object
traj_sce <- sim_trajectory[['SCE']]

## Extract sample metadata to use for testing
traj_meta <- sim_trajectory[["meta"]]

## Add metadata to colData slot
colData(traj_sce) <- DataFrame(traj_meta)

data("sim_trajectory", package = "miloR")

## Extract SingleCellExperiment object
traj_sce <- sim_trajectory[['SCE']]

## Extract sample metadata to use for testing
traj_meta <- sim_trajectory[["meta"]]

## Add metadata to colData slot
colData(traj_sce) <- DataFrame(traj_meta)

logcounts(traj_sce) <- log(counts(traj_sce) + 1)
traj_sce <- runPCA(traj_sce, ncomponents=30)
traj_sce <- runUMAP(traj_sce)

plotUMAP(traj_sce)

traj_milo <- Milo(traj_sce)
reducedDim(traj_milo, "UMAP") <- reducedDim(traj_sce, "UMAP")

traj_milo

traj_milo <- buildGraph(traj_milo, k = 10, d = 30)

traj_milo <- makeNhoods(traj_milo, prop = 0.1, k = 10, d=30, refined = TRUE)
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample="Sample")

head(nhoodCounts(traj_milo))

traj_design <- data.frame(colData(traj_milo))[,c("Sample", "Condition")]
traj_design <- distinct(traj_design)

traj_design

rownames(traj_design) <- traj_design$Sample

da_results <- testNhoods(traj_milo, design = ~ Condition, design.df = traj_design)

da_results %>%
  arrange(- SpatialFDR) %>%
  head()

traj_milo <- buildNhoodGraph(traj_milo)


plotUMAP(traj_milo) + plotNhoodGraphDA(traj_milo, da_results, alpha=0.05) +
  plot_layout(guides="collect")

# Assign cells to clones (input long donor id)

args <- commandArgs(TRUE)
if (length(args) >= 1) {
  don <- args[1] #donor
}

library(tidyverse)
library(variancePartition)
library(SingleCellExperiment)


## Load data
params <- list()
params$callset <- "filt_lenient.cell_coverage_sites"

sce <- readRDS(file.path("data/sces",
    paste0("sce_", don, "_with_clone_assignments.", params$callset, ".rds")))
sce_unst <- sce[, sce$well_condition == "unstimulated"]
cat(paste("reading", don, ":   ", ncol(sce), "cells total,",
    ncol(sce_unst), "cells unstimulated.\n"))
rm(sce)


## Filtering
sce_unst$cdr <- (sce_unst$total_features / nrow(sce_unst))
sce_unst$plate <- as.factor(sce_unst$plate)
sce_unst$assigned <- as.factor(sce_unst$assigned)
cell_idx <- which(sce_unst$assigned != "unassigned")

var_min <- 1.0
gene_var <- rowVars(logcounts(sce_unst)[,cell_idx])
gene_idx <- which(!sce_unst@int_elementMetadata$is_spike_ERCC &
                    gene_var > var_min)
n_genes <- length(gene_idx)
gene_idx <- gene_idx[order(gene_var[gene_idx], decreasing = TRUE)][1:n_genes]
print(c(n_genes, min(gene_var[gene_idx])))

# input data
log_cnt <- logcounts(sce_unst)[gene_idx, cell_idx, drop = FALSE]
covar <- as.data.frame(colData(sce_unst)[cell_idx,])
row.names(log_cnt)[1:5]

## fitting
form <- ~ cdr + (1|plate) + (1|assigned)
varPart_mat <- matrix(NA, nrow = n_genes, ncol = 4)
for (j in seq_len(n_genes)) {
  try(varPart_mat[j,] <-
        unlist(fitExtractVarPartModel(log_cnt[j, , drop = FALSE], form, covar)))
}
varPart_tmp <- fitExtractVarPartModel(log_cnt[j, , drop = FALSE], form, covar)

colnames(varPart_mat) <- colnames(varPart_tmp)
row.names(varPart_mat) <- row.names(log_cnt)

head(varPart_mat)

## save data
dat_dir <- "data/variance_components/donorVar/"
write.table(varPart_mat, file = paste0(dat_dir, don, ".var_part.var1.csv"),
            sep = ",")

cat(paste(don, ":   ", length(cell_idx), "cells ,", length(gene_idx),
          "genes.\n"))


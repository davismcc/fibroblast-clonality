params <-
structure(list(rdata_file = "test.RData", author = "Davis McCarthy", 
    title = "Rough QC Report", to_working_dir = "./"), .Names = c("rdata_file", 
"author", "title", "to_working_dir"))

## ----load-data, results='hide', message=FALSE----------------------------
library(scater, quietly = TRUE)
library(knitr, quietly = TRUE)
library(ggplot2)
library(viridis)
library(ggthemes, quietly = TRUE)
library(cowplot)
library(magrittr)
library(destiny)
## options
message(cat(paste("Working dir: ", getwd(), "\n")))
message(cat(params$rdata_file))
# opts_knit$set(root.dir = as.character(params$to_working_dir))
#opts_knit$set(root.dir = "../../")
opts_chunk$set(fig.align = 'center', fig.width = 10, fig.height = 5, dev = 'png',
               warning = FALSE, message = FALSE)
message(cat(paste("New Working dir: ", getwd(), "\n")))
#setwd(params$working_dir)
## load data
rdata_file <- file.path(params$to_working_dir, params$rdata_file)
sce <- readRDS(rdata_file)

## ----filter-kallist-no-exprs, results = 'hide'---------------------------
nonzero_genes <- (rowSums(counts(sce)) > 0)
sum(nonzero_genes)
sce_filt <- sce[nonzero_genes,]

## ----calculate-qc-metrics------------------------------------------------
count_thresh <- 0
is_exprs(sce_filt) <- (counts(sce_filt) > count_thresh)
ercc_genes <- grepl("^ERCC", featureNames(sce_filt))
mt_genes <- grepl("^MT-", fData(sce_filt)$hgnc_symbol)
sce_filt <- calculateQCMetrics(sce_filt, 
    feature_controls = list(ERCC = ercc_genes, MT = mt_genes))

## ----summary-ercc-pct----------------------------------------------------
stats <- summary(sce_filt$pct_counts_feature_controls_ERCC)
kable(data.frame("Stat" = names(stats), "Pct Counts from ERCC" = as.vector(stats)))
stats_MT <- summary(sce_filt$pct_counts_feature_controls_MT)
kable(data.frame("Stat" = names(stats_MT), "Pct Counts from MT" = as.vector(stats_MT)))

## ----plot-pheno-1--------------------------------------------------------
p1 <- plotPhenoData(sce_filt, 
              aes(x = log10_total_counts, y = total_features, 
                  colour =  filter_on_total_features)) +
    theme(legend.position = "bottom") +
    geom_vline(xintercept = 5, linetype = 2)
p2 <- plotPhenoData(sce_filt, 
              aes(x = log10_total_counts, y = total_features, 
                  colour =  pct_counts_top_100_features)) +
    theme(legend.position = "bottom") +
    geom_vline(xintercept = 5, linetype = 2)
p3 <- plotPhenoData(sce_filt, 
              aes(x = log10_total_counts, y = total_features, 
                  colour =  pct_counts_feature_controls_ERCC)) +
    theme(legend.position = "bottom") +
    geom_vline(xintercept = 5, linetype = 2)
plot_grid(p1, p2, p3, labels = c("A", "B", "C"), 
          nrow = 1)

## ----plot-pheno-dropout, echo = FALSE------------------------------------
p1 <- plotPhenoData(sce_filt, 
              aes(x = log10(pct_counts_feature_controls_ERCC), y = pct_dropout, 
                  colour =  pct_counts_feature_controls_MT)) +
    theme(legend.position = "bottom") + 
    geom_vline(xintercept = 0, linetype = 2)
p2 <- plotPhenoData(sce_filt, 
              aes(x = log10(pct_counts_feature_controls_ERCC), y = pct_dropout,
                  colour =  pct_counts_top_100_features)) +
    theme(legend.position = "bottom") + 
    geom_vline(xintercept = 0, linetype = 2)
p3 <- plotPhenoData(sce_filt, 
              aes(x = log10(pct_counts_feature_controls_ERCC), y = pct_dropout,
                  colour =  pct_counts_feature_controls_ERCC)) +
    theme(legend.position = "bottom") + 
    geom_vline(xintercept = 0, linetype = 2)
plot_grid(p1, p2, p3, labels = c("A", "B", "C"), nrow = 1)

## ----plot-pheno-dropout-MT-genes, echo = FALSE---------------------------
p1 <- plotPhenoData(sce_filt, 
              aes(x = pct_counts_feature_controls_MT, y = pct_dropout, 
                  colour =  pct_counts_feature_controls_MT)) +
    theme(legend.position = "bottom") + 
    geom_vline(xintercept = 0, linetype = 2)
p2 <- plotPhenoData(sce_filt, 
              aes(x = pct_counts_feature_controls_MT, y = pct_dropout,
                  colour =  pct_counts_top_100_features)) +
    theme(legend.position = "bottom") + 
    geom_vline(xintercept = 0, linetype = 2)
p3 <- plotPhenoData(sce_filt, 
              aes(x = pct_counts_feature_controls_MT, y = pct_dropout,
                  colour =  pct_counts_feature_controls_ERCC)) +
    theme(legend.position = "bottom") + 
    geom_vline(xintercept = 0, linetype = 2)
plot_grid(p1, p2, p3, labels = c("A", "B", "C"), nrow = 1)

## ----compare-overall-cumul-plots, fig.height=8---------------------------
p1 <- plot(sce_filt, exprs_values = "counts",
           colour_by = "pct_exprs_feature_controls_ERCC") +
    theme(legend.position = "bottom")
p2 <- plot(sce_filt, exprs_values = "counts",
           colour_by = "pct_exprs_feature_controls_MT") +
    theme(legend.position = "bottom")
plot_grid(p1, p2, labels = c("A", "B"), nrow = 1)

## ----plot-QC-most-exprs--------------------------------------------------
plotQC(sce_filt)

## ----plot-QC-exprs-freq-vs-mean, echo=FALSE------------------------------
plotQC(sce_filt, type = "exprs") +
    theme(legend.position = "right")

## ----pca-qc-metrics, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE, fig.width=8, fig.height=5, fig.align="center"----
plotPCA(sce_filt, size_by = "total_features", 
        shape_by = "filter_on_total_features", 
        pca_data_input = "pdata", detect_outliers = TRUE)

## ----pca-qc, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE, fig.width=10, fig.height=5, fig.align="center"----
plotPCA(sce_filt, size_by = "total_features", 
        colour_by = "pct_counts_feature_controls_ERCC")

## ----plot-QC-find-pcs-total-features, fig.height=7-----------------------
plotQC(sce_filt, type = "find-pcs", 
       variable = "total_features") 

## ----find-pcs-pct-top-50, fig.height=7-----------------------------------
plotQC(sce_filt, type = "find-pcs", 
       variable = "pct_counts_top_50_features") 

## ----expl-vars-----------------------------------------------------------
sce_filt$start_time <- NULL
zero_var <- matrixStats::rowVars(exprs(sce_filt)) == 0
plotQC(sce_filt[!zero_var,], "expl",
       variables = c("pct_dropout", "total_features", 
                     "pct_counts_top_200_features", 
                     "pct_counts_feature_controls_ERCC", 
                     "pct_counts_feature_controls_MT", 
                     "n_detected_feature_controls",
                     "log10_counts_endogenous_features",
                     "log10_counts_feature_controls_ERCC",
                     "log10_counts_feature_controls_MT"))


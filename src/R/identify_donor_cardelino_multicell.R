"
Identify donor multiple cells from a filtered VCF file

Usage:
identify_donor_small_vcf_cardelino.R --input_file <DATA_IN> --output_prefix <PRFX_OUT> [--donor_vcf <DONOR_VCF> -hv]

Options:
-i --input_file <DATA_IN>      input file
-o --output_prefix <PRFX_OUT>  prefix for output files; .RData, .feather and .csv files produced
-d --donor_vcf <DONOR_VCF>     VCF file containing the genotype data for all donors [default: /nfs/research2/stegle/projects/hipsci/data/genotypes/imputed/REL-2014-11_SS/hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.20151005.genotypes.gdid.mac1.recode.WGS.maf0.1.vcf.gz]
-h --help                      show this
-v --version                   print version and stop

This program does compares called genotypes from RNA-seq reads for a cell to
genotype data for all HipSci donors to find the best match in terms of
(weighted) average genotypic correlation across called variants.

The program returns results as a .csv file.

This program requires the R packages 'VariantAnnotation' and 'snpStats'
(Bioconductor) and 'cardelino' (github.com/davismcc/cardelino).

Davis McCarthy
May 2018
" -> doc

## Script to identify HIPSCI donor for a cell's VCF file
# suppressPackageStartupMessages(library(VariantAnnotation))
# suppressPackageStartupMessages(library(snpStats))
suppressPackageStartupMessages(library(cardelino))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(pcaMethods))

main <- function(input_vcf, output_prefix, donor_vcf) {
    ## Read in VCFs
    cell_vcf <- read_vcf(input_vcf)
    donor_vcf <- read_vcf(donor_vcf)
    ## Donor ID: first pass without checking for doublets
    ids_sing <- donor_id(cell_vcf, donor_vcf, check_doublet = FALSE,
                            n_vars_threshold = 25)
    ## Donor ID: second pass checking for doublets on subset of donors
    if (all(ids_sing$assigned$donor_id == "unassigned")) {
        ids_doub <- ids_sing
    } else {
        ids_doub <- donor_id(cell_vcf, donor_vcf, check_doublet = TRUE,
                            donors = unique(ids_sing$assigned$donor_id))
        p1 <- ggplot(ids_doub$assigned, aes(n_vars, prob_doublet, colour = donor_id)) +
                geom_point(size = 3, alpha = 0.5) +
                theme_bw()
        ggsave(paste0(output_prefix, "_pdoublet.png"), plot = p1, width = 9, 
                height = 8)
    }
    ## ppca on cell genotypes
    pp <- pcaMethods::ppca(t(ids_doub$A / ids_doub$D))
    df <- ids_doub$assigned
    df$PPCA1 <- pp@scores[, 1]
    df$PPCA2 <- pp@scores[, 2]
    p2 <- ggplot(df, aes(PPCA1, PPCA2, colour = donor_id)) +
        geom_point(alpha = 0.5) +
        theme_bw()
    ggsave(paste0(output_prefix, "_ppca.png"), plot = p2, width = 9, height = 8)
    ## setup output df
    output_df <- ids_doub$assigned
    output_df[["sample_id"]] <- colnames(ids_doub$D)
    output_df[["n_total_reads"]] <- colSums(ids_doub$D, na.rm = TRUE)
    output_df[["n_alt_reads"]] <- colSums(ids_doub$A, na.rm = TRUE)
    ## write output to file
    message("Writing output to CSV\n")
    write.csv(output_df, file = paste0(output_prefix, ".csv"), 
                row.names = FALSE)
    return("Done.")
}

## Get command line options
opt <- docopt::docopt(doc, version = "version 0.0.1\n")

message("working directory: ", getwd(), "\n")
message("input vcf: ", opt$input_file, "\n")
message("output prefix: ", opt$output_prefix, "\n")
message("donor vcf: ", opt$donor_vcf, "\n")

## Run main function
main(opt$input_file, opt$output_prefix, opt$donor_vcf)

# # params for testing
# opt <-  list()
# opt[["input_file"]] <- "Data/SS2_2017/22782_5/cells_merged_donor_id_filt.vcf.gz"
# opt[["output_prefix"]] <- "Data/SS2_2017/22782_5/donor_id_all_cardelino_v2"
# opt[["donor_vcf"]] <- "Data/SS2_2017/22782_5/filtered.hipsci.overlap.vcf.gz"

# opt[["input_file"]] <- "../hipsci-fibro/Data/SS2_2017/22782_6.cells_merged_donor_id_filt.vcf.gz"
# opt[["output_prefix"]] <- "../hipsci-fibro/Data/SS2_2017/22782_6.donor_id_all_cardelino_v2"
# opt[["donor_vcf"]] <- "../hipsci-fibro/Data/SS2_2017/22782_6.filtered.hipsci.overlap.vcf.gz"

# opt <-  list()opt[["input_file"]] <- "data_raw/scrnaseq/run_24252/cells_merged_filt.vcf.gz"
# 
# opt[["output_prefix"]] <- "data_raw/scrnaseq/run_24252/tmp.donor_id_cardelino_all"
# opt[["donor_vcf"]] <- "data_raw/scrnaseq/run_24252/filtered.hipsci.overlap.vcf.gz"
# opt[["donor_lines"]] <- "aowh_2;keui_1;meue_4;naah_2;poih_4;vils_1;bokz_5;datg_2;guss_1;nudd_1;sehl_6"

# opt <-  list()
# opt[["input_file"]] <- "data_raw/scrnaseq/run_25579/cells_merged_filt.vcf.gz"
# opt[["output_prefix"]] <- "data_raw/scrnaseq/run_25579/tmp.donor_id_cardelino_all"
# opt[["donor_vcf"]] <- "data_raw/scrnaseq/run_25579/filtered.hipsci.overlap.vcf.gz"
# opt[["donor_lines"]] <- "qaqx;sojd;rayr;vass;yemz;tolg;ciwj;hajc;hecn;kuco;liqa;tert"


# read_sample_vcf <- function(input_vcf) {
#     ## Read in VCF from this sample
#     message("Reading sample VCF\n")
#     vcf_sample <- readVcf(input_vcf, "GRCh37")
#     vcf_sample_filt <- vcf_sample[isSNV(vcf_sample)]
#     if (length(vcf_sample_filt) > 0) {
#         seqlevelsStyle(vcf_sample) <- "UCSC"
#         new_snp_names <- paste0("snp_",
#                                 gsub("chr", "",
#                                      gsub(":", "_",
#                                           gsub("_[ATCG]/[ATCG]", "",
#                                                names(vcf_sample_filt)))))    
#         names(vcf_sample_filt) <- new_snp_names
#     }
#     vcf_sample_filt
# }

# read_donor_vcf <- function(donor_vcf) {
#     message("Reading donors VCF\n")
#     vcf_donor <- readVcf(donor_vcf, "GRCh37")
#     seqlevelsStyle(vcf_donor) <- "UCSC"
#     isSNV_idx <- isSNV(vcf_donor)
#     vcf_donor <- vcf_donor[, grep("HPS.*pf-[a-z]+$", colnames(vcf_donor))]
#     colnames(vcf_donor) <- gsub("HPS.*-", "", gsub("_[0-9]+", "", colnames(vcf_donor)))
#     list(vcf = vcf_donor, isSNV_idx = isSNV_idx)
# }

# get_snp_matrices <- function(vcf_sample, vcf_donor) {
#     ## get snp matrices
#     sm_donor <- genotypeToSnpMatrix(
#         geno(vcf_donor, "GT"), ref = ref(vcf_donor), alt = alt(vcf_donor))
#     ## filter sample VCF to those variants found in donor VCF
#     vcf_sample <- sortSeqlevels(vcf_sample, X.is.sexchrom = TRUE)
#     slengths_sample <- seqlengths(vcf_sample)
#     vcf_donor <- sortSeqlevels(vcf_donor, X.is.sexchrom = TRUE)
#     seqlengths(vcf_donor) <- slengths_sample[seqlevels(vcf_donor)]
#     ovlap <- findOverlaps(vcf_sample, vcf_donor)
#     if (length(ovlap) < 1L) {
#         message("No common variants overlapping in sample VCF and Donor VCF\n")
#         return(list(stop_program = TRUE))
#     } else {
#         message("Filtering sample VCF\n")
#         vcf_sample2 <- vcf_sample[queryHits(ovlap)]
#         message("Filtering donor VCF\n")
#         vcf_donor <- vcf_donor[subjectHits(ovlap)]
#         match_alleles <- unlist(ref(vcf_sample2) == ref(vcf_donor) & alt(vcf_sample2) == alt(vcf_donor))
#         if (sum(match_alleles) < 1L)
#             stop("No variants with matching alleles in sample and donor VCFs")
#         else {
#             vcf_sample2 <- vcf_sample2[match_alleles]
#             vcf_donor <- vcf_donor[match_alleles]
#         }
#         message("Extracting sample SNP matrix\n")
#         sm_sample <- genotypeToSnpMatrix(
#                 geno(vcf_sample2, "GT"), ref = ref(vcf_sample2),
#                 alt = alt(vcf_sample2))
#         sm_sample_REF <- matrix(sapply(geno(vcf_sample2, "AD"), 
#                             function(x) x[[1]]), ncol = ncol(vcf_sample2))
#         sm_sample_ALT <- matrix(sapply(geno(vcf_sample2, "AD"), 
#                             function(x) x[[2]]), ncol = ncol(vcf_sample2))
#         sm_sample_ALT[is.na(sm_sample_ALT) & !is.na(sm_sample_REF)] <- 0
#         sm_sample_DEP <- sm_sample_REF + sm_sample_ALT
#         sm_sample_DEP[sm_sample_DEP == 0] <- NA
#         sm_sample_REF[is.na(sm_sample_DEP)] <- NA
#         sm_sample_ALT[is.na(sm_sample_DEP)] <- NA
#         rownames(sm_sample_REF) <- rownames(sm_sample_ALT) <- 
#             rownames(sm_sample_DEP) <- rownames(vcf_sample2)
#         colnames(sm_sample_REF) <- colnames(sm_sample_ALT) <- 
#             colnames(sm_sample_DEP) <- colnames(vcf_sample2)
#         na_sample <- is.na(sm_sample_DEP)
#         if (sum(!na_sample) < 1L) {
#             message("No common variants with non-missing genotypes overlapping in sample VCF and Donor VCF\n")
#             return(list(stop_program = TRUE))
#         } else {
#             message("Extracting donor SNP matrix\n")
#             sm_donor <- genotypeToSnpMatrix(
#                 geno(vcf_donor, "GT"), ref = ref(vcf_donor),
#                 alt = alt(vcf_donor))
#             donor_geno_mat <- matrix(as.numeric(as(sm_donor$genotypes, "numeric") > 0), 
#                 nrow = nrow(sm_donor$genotypes))
#             rownames(donor_geno_mat) <- rownames(sm_donor$genotypes)
#             colnames(donor_geno_mat) <- colnames(sm_donor$genotypes)
#             donor_geno_mat <- t(donor_geno_mat[order(rownames(donor_geno_mat)),])
#             message("Doing donor assignment using ", length(vcf_sample2),
#                     " variants\n")
#             return(list(sm_sample_REF = sm_sample_REF, A = sm_sample_ALT, 
#                     D = sm_sample_DEP, sm_sample = sm_sample, 
#                     stop_program = FALSE, C = donor_geno_mat))
#         }
#     }
# }


# main <- function(input_vcf, output_prefix, donor_vcf, missing_thresh = 0.01,
#                     maf_thresh = 0.01) {
#     ## Read in VCF from this sample
#     vcf_sample <- read_sample_vcf(input_vcf)
#     if (length(vcf_sample) < 1) {
#         write.csv(output_df, file = paste0(output_prefix, ".csv"),
#                   row.names = FALSE)
#         message("No variants in sample VCF after filtering.\n")
#         return("Done.")
#     }
#     message("...read ", length(vcf_sample), " variants from sample VCF\n")
#     if (length(vcf_sample) > 10000) {
#         vcf_sample <- vcf_sample[sample.int(nrow(vcf_sample), 10000)]
#         message(".......using sample of 10,000 variants from sample VCF\n")
#     }
#     ## Read in Donor VCF
#     donor_data <- read_donor_vcf(donor_vcf)
#     isSNV_idx <- donor_data$isSNV_idx
#     vcf_donor <- donor_data$vcf
#     rm(donor_data)
#     message("...read ", length(vcf_donor), " variants from donor VCF\n")
#     if (!any(isSNV_idx)) {
#         write.csv(output_df, file = paste0(output_prefix, ".csv"),
#                   row.names = FALSE)
#         message("No single-nucleotide variants overlapping in sample VCF and Donor VCF\n")
#         return("Done.")
#     }
#     vcf_donor <- vcf_donor[isSNV_idx]
#     ## get snp matrices
#     snpmat_list <- get_snp_matrices(vcf_sample, vcf_donor)
#     rm(vcf_sample)
#     rm(vcf_donor)
#     if (snpmat_list$stop_program) {
#         write.csv(output_df, file = paste0(output_prefix, ".csv"),
#                   row.names = FALSE)
#         message("No common variants overlapping in sample VCF and Donor VCF\n")
#         return("Done.")
#     }
#     ## Add extra filtering step 
#     keep_snp <- (rowMeans(!is.na(snpmat_list$D)) > missing_thresh &
#         (rowSums(snpmat_list$A, na.rm = TRUE) / 
#             rowSums(snpmat_list$D, na.rm = TRUE) > maf_thresh))
#     snpmat_list$A <- snpmat_list$A[keep_snp,]
#     snpmat_list$D <- snpmat_list$D[keep_snp,]
#     snpmat_list$C <- snpmat_list$C[keep_snp,]
#     ## setup output df
#     output_df <- data.frame(
#         sample_id = colnames(snpmat_list$D),
#         donor = NA,
#         post_prob = NA,
#         second_best_donor = NA,
#         second_best_post_prob = NA,
#         nvars_used = colSums(snpmat_list$D > 0.5, na.rm = TRUE),
#         n_total_reads = colSums(snpmat_list$D, na.rm = TRUE),
#         n_alt_reads = colSums(snpmat_list$A, na.rm = TRUE),
#         stringsAsFactors = FALSE)
#     assign <- cell_assign_EM(
#         A = snpmat_list$A, D = snpmat_list$D, C = snpmat_list$C,
#         Psi = rep(1 / ncol(snpmat_list$C), ncol(snpmat_list$C)),
#         model = "Binomial")
#     probs <- assign$prob
#     donors <- colnames(probs)
#     for (i in seq_len(nrow(probs))) {
#         o <- order(probs[i,], decreasing = TRUE)
#         output_df[["donor"]][i] <- donors[o[1]]
#         output_df[["post_prob"]][i] <- probs[i, o[1]]
#         output_df[["second_best_donor"]][i] <- donors[o[2]]
#         output_df[["second_best_post_prob"]][i] <- probs[i, o[2]]
#     }
#     ## write output to file
#     message("Writing output to CSV\n")
#     write.csv(output_df, file = paste0(output_prefix, ".csv"), 
#                 row.names = FALSE)
#     return("Done.")
# }



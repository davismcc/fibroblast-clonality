"
Identify donor from a cell's filtered VCF file

Usage:
identify_donor.R --input_file <DATA_IN> --output_prefix <PRFX_OUT> --donor_lines <DONORS> [--donor_vcf <DONOR_VCF> --fasta_idx <FASTA_IDX> -hv]

Options:
-i --input_file <DATA_IN>      input file
-o --output_prefix <PRFX_OUT>  prefix for output files; .RData, .feather and .csv files produced
-l --donor_lines <DONORS>      string providing donor line IDs separated by ';' (e.g. 'babk_2;wuye_2;podx3')
-d --donor_vcf <DONOR_VCF>     VCF file containing the genotype data for all donors [default: /nfs/research2/stegle/projects/hipsci/data/genotypes/imputed/REL-2014-11_SS/hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.20151005.genotypes.gdid.mac1.recode.WGS.maf0.1.vcf.gz]
-h --help                      show this
-v --version                   print version and stop

This program does compares called genotypes from RNA-seq reads for a cell to
genotype data for all HipSci donors to find the best match in terms of
(weighted) average genotypic correlation across called variants.

The program returns results as a .csv file.

This program requires the R packages 'VariantAnnotation' and 'snpStats'
(Bioconductor).

Davis McCarthy
March 2017
" -> doc

## Script to identify HIPSCI donor for a cell's VCF file
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(snpStats))

init_output_df <- function(input_vcf, ave_allelic_corr,
                           wght_ave_allelic_corr, these_donors,
                           short_names) {
    ## define output data frame
    output_df <- data.frame(
        sample_id = gsub(".filtered.vcf", "", basename(input_vcf)),
        donor_long_id = as.character(names(ave_allelic_corr)),
        donor_short_id = short_names,
        nvars_called = 0,
        nvars_ovlap_hipsci = 0,
        nvars_used = 0,
        ave_allelic_corr = ave_allelic_corr,
        wght_ave_allelic_corr = wght_ave_allelic_corr,
        pctile_ave_corr = NA,
        pctile_wght_ave_corr = NA,
        used_in_expt = (short_names %in% these_donors),
        best_match = FALSE, stringsAsFactors = FALSE)
    output_df
}

full_output_df <- function(input_vcf, vcf_sample, vcf_hipsci, sm_sample,
                           ave_allelic_corr, wght_ave_allelic_corr,
                           these_donors, short_names) {
    ## use ecdf of correlation values to get confidence for called donor
    ## as percentile of distribution of scores across all HipSci donors
    f_ave_corr <- ecdf(ave_allelic_corr)
    f_wght_ave_corr <- ecdf(wght_ave_allelic_corr)
    ## define output data frame
    output_df <- data.frame(
        sample_id = gsub("_filtered.vcf", "", basename(input_vcf)),
        donor_long_id = as.character(names(ave_allelic_corr)),
        donor_short_id = short_names,
        nvars_called = length(vcf_sample),
        nvars_ovlap_hipsci = length(vcf_hipsci),
        nvars_used = ncol(sm_sample$genotypes),
        ave_allelic_corr = ave_allelic_corr,
        wght_ave_allelic_corr = wght_ave_allelic_corr,
        pctile_ave_corr = f_ave_corr(ave_allelic_corr),
        pctile_wght_ave_corr = f_wght_ave_corr(wght_ave_allelic_corr),
        used_in_expt = (short_names %in% these_donors),
        best_match = FALSE, stringsAsFactors = FALSE)
    output_df
}

compute_confidence_score <- function(output_df, these_donors, short_names) {
     ## compute a confidence score
    wd <- which(short_names %in% these_donors)
    conf_ave <- conf_wght_ave <- rep(NA, length(wd))
    names(conf_ave) <- names(conf_wght_ave) <- short_names[wd]
    for (i in seq_along(conf_ave)) {
        this_ave <- output_df$ave_allelic_corr[wd[i]]
        max_ave <- max(output_df$ave_allelic_corr[wd[-i]])
        conf_ave[i] <- (this_ave - max_ave) / abs(max_ave)
        this_wght_ave <- output_df$wght_ave_allelic_corr[wd[i]]
        max_wght_ave <- max(output_df$wght_ave_allelic_corr[wd[-i]])
        conf_wght_ave[i] <- (this_wght_ave - max_wght_ave) / abs(max_wght_ave)
    }
    output_df$enrich_ave_corr <- output_df$enrich_wght_ave_corr <- NA
    output_df$enrich_ave_corr[wd] <- conf_ave
    output_df$enrich_wght_ave_corr[wd] <- conf_wght_ave
    output_df
}

read_sample_vcf <- function(input_vcf) {
    ## Read in VCF from this sample
    message("Reading sample VCF\n")
    vcf_sample <- readVcf(input_vcf, "GRCh37")
    keep_seqs <- seqnames(vcf_sample) %in% paste0("chr", c(1:22))
    vcf_sample_filt <- vcf_sample[keep_seqs]
    vcf_sample_filt <- vcf_sample_filt[isSNV(vcf_sample_filt)]
    if (length(vcf_sample_filt) > 0) {
        vcf_sample_filt <- renameSeqlevels(
            vcf_sample_filt, gsub("chr", "", seqlevels(vcf_sample_filt)))
        new_snp_names <- paste0("snp_",
                                gsub("chr", "",
                                     gsub(":", "_",
                                          gsub("_[ATCG]/[ATCG]", "",
                                               names(vcf_sample_filt)))))    
        names(vcf_sample_filt) <- new_snp_names
    }
    vcf_sample_filt
}

read_donor_vcf <- function(donor_vcf) {
    message("Reading donors VCF\n")
    vcf_donor <- readVcf(donor_vcf, "GRCh37")
    isSNV_idx <- isSNV(vcf_donor)
    list(vcf = vcf_donor, isSNV_idx = isSNV_idx)
}

get_snp_matrices <- function(vcf_sample, vcf_donor, short_names) {
    ## get snp matrices
    sm_donor <- genotypeToSnpMatrix(
        geno(vcf_donor, "GT"), ref = ref(vcf_donor), alt = alt(vcf_donor))
    donor <- gsub("_[0-9]+", "", short_names)
    donor <- gsub("2:", "",
                  gsub("CTRL0914es-", "",
                       gsub("CTRL1114es-", "", donor)))
    lines_for_stats <- !duplicated(donor)
    snp_stats_donor <- col.summary(sm_donor$genotypes[lines_for_stats,])
    common_vars <- snp_stats_donor$MAF > 0.02
    common_vars[is.na(common_vars)] <- FALSE
    if (sum(common_vars, na.rm = TRUE) == 0) {
        message("No common variants overlapping in sample VCF and Donor VCF\n")
        return(list(stop_program = TRUE))
    }
    vcf_donor2 <- vcf_donor[common_vars]
    sm_donor <- genotypeToSnpMatrix(
        geno(vcf_donor2, "GT"), ref = ref(vcf_donor2),
        alt = alt(vcf_donor2))
    snp_stats_donor <- col.summary(sm_donor$genotypes[lines_for_stats,])
    ## filter sample VCF to those variants found in donor VCF
    message("Filtering sample VCF\n")
    slengths_sample <- seqlengths(vcf_sample)[1:25]
    names(slengths_sample)[25] <- "MT"
    seqlengths(vcf_donor2) <- slengths_sample[seqlevels(vcf_donor2)]
    ovlap <- findOverlaps(vcf_sample, vcf_donor2)
    if (length(ovlap) < 1L) {
        message("No common variants overlapping in sample VCF and Donor VCF\n")
        return(list(stop_program = TRUE))
    } else {
        vcf_sample2 <- vcf_sample[queryHits(ovlap)]
        sm_sample <- genotypeToSnpMatrix(geno(vcf_sample2, "GT"),
                                         ref = ref(vcf_sample2),
                                         alt = alt(vcf_sample2))
        sm_sample$map$quality <- qual(vcf_sample2)
        sm_sample$map$weight <- qual(vcf_sample2) / sum(qual(vcf_sample2))
        na_sample <- is.na(as(sm_sample$genotypes, "numeric"))
        if (sum(!na_sample) < 1L) {
            message("No common variants with non-missing genotypes overlapping in sample VCF and Donor VCF\n")
            return(list(stop_program = TRUE))
        } else {
            vcf_donor <- vcf_donor2[subjectHits(ovlap)]
            sm_donor <- genotypeToSnpMatrix(
                geno(vcf_donor, "GT"), ref = ref(vcf_donor),
                alt = alt(vcf_donor))
            if ( !identical(as.character(sm_sample$map$allele.1),
                            as.character(sm_donor$map$allele.1)))
                stop("Alleles do not match between sample and donor VCFs")
            message("Computing relatedness using ", length(vcf_sample2),
                    " variants\n")
            return(list(sm_sample = sm_sample, sm_donor = sm_donor,
                        stop_program = FALSE))
        }
    }
}

compute_relatedness_score <- function(sm_sample, sm_donor, short_names) {
    donor <- gsub("_[0-9]+", "", short_names)
    donor <- gsub("2:", "",
                  gsub("CTRL0914es-", "",
                       gsub("CTRL1114es-", "", donor)))
    snp_stats_donor <- col.summary(sm_donor$genotypes[!duplicated(donor),])
    ## standardise sample genotypes
    message("Standardise genotypes\n")
    donor_geno_mat <- as(sm_donor$genotype, "numeric")
    donor_geno_mat[is.na(donor_geno_mat)] <- 0
    ## standardise sample genotypes
    sample_geno_vec <- as(sm_sample$genotype, "numeric")[1,]
    sm_sample_stand_geno <- ( (sample_geno_vec - 2 * snp_stats_donor$RAF)
        / sqrt(2 * snp_stats_donor$RAF * (1 - snp_stats_donor$RAF)) )
    ## standardise donor genotypes
    sm_donor_stand_geno <- t( (t(donor_geno_mat) - 2 * snp_stats_donor$RAF)
    / sqrt(2 * snp_stats_donor$RAF * (1 - snp_stats_donor$RAF)))
    ## compute weighted correlation between genotypes of each donor donor and
    ## the sample
    message("Compute genotypic correlations\n")
    corr_mat <- t(t(sm_donor_stand_geno) * sm_sample_stand_geno)
    wght_corr_mat <- t(t(sm_donor_stand_geno) * sm_sample_stand_geno *
                       sm_sample$map$weight)
    bad_snp <- apply(corr_mat, 2, function(x) any(is.na(x)))
    corr_mat <- corr_mat[, !bad_snp, drop = FALSE]
    wght_corr_mat <- wght_corr_mat[, !bad_snp, drop = FALSE]
    ## sum correlations across donors
    ave_allelic_corr <- rowMeans(corr_mat, na.rm = TRUE)
    wght_ave_allelic_corr <- rowSums(wght_corr_mat, na.rm = TRUE)
    list(ave_allelic_corr = ave_allelic_corr,
         wght_ave_allelic_corr = wght_ave_allelic_corr)
}

main <- function(input_vcf, output_prefix, donor_lines, donor_vcf, fasta_idx) {
    ## define samples in donor VCF
    hdr_donor <- scanVcfHeader(donor_vcf)
    ave_allelic_corr <- wght_ave_allelic_corr <-
        rep(0, length(samples(hdr_donor)))
    names(ave_allelic_corr) <- names(wght_ave_allelic_corr) <-
        samples(hdr_donor)
    ## define donors and names
    short_names <- gsub("HPSI[0-9a-z]+-", "", names(ave_allelic_corr))
    these_donors <- strsplit(donor_lines, ";")[[1]]
    these_donors_idx <- grepl(gsub(";", "|", donor_lines),
                              names(ave_allelic_corr))
    ## define output data frame
    output_df <- init_output_df(input_vcf, ave_allelic_corr,
                                wght_ave_allelic_corr, these_donors,
                                short_names)
    ## Read in VCF from this sample
    vcf_sample <- read_sample_vcf(input_vcf)
    if (length(vcf_sample) < 1) {
        write.csv(output_df, file = paste0(output_prefix, ".csv"),
                  row.names = FALSE)
        message("No variants in sample VCF after filtering.\n")
        return("Done.")
    }
    output_df$nvars_called <- length(vcf_sample)
    message("...read ", length(vcf_sample), " variants from sample VCF\n")
    ## Read in Donor VCF
    donor_data <- read_donor_vcf(donor_vcf)
    isSNV_idx <- donor_data$isSNV_idx
    vcf_donor <- donor_data$vcf
    if (!any(isSNV_idx)) {
        write.csv(output_df, file = paste0(output_prefix, ".csv"),
                  row.names = FALSE)
        message("No single-nucleotide variants overlapping in sample VCF and Donor VCF\n")
        return("Done.")
    }
    vcf_donor <- vcf_donor[isSNV_idx]
    ## get snp matrices
    snpmat_list <- get_snp_matrices(vcf_sample, vcf_donor, short_names)
    if (snpmat_list$stop_program) {
        write.csv(output_df, file = paste0(output_prefix, ".csv"),
                  row.names = FALSE)
        message("No common variants overlapping in sample VCF and Donor VCF\n")
        return("Done.")
    }
    sm_sample <- snpmat_list$sm_sample
    sm_donor <- snpmat_list$sm_donor
    rel_scores <- compute_relatedness_score(sm_sample, sm_donor, short_names)
    output_df <- full_output_df(input_vcf, vcf_sample, vcf_donor, sm_sample,
                                rel_scores$ave_allelic_corr,
                                rel_scores$wght_ave_allelic_corr,
                                these_donors, short_names)    
    ## define best match from those donors known to have been used
    wm <- which.max(ave_allelic_corr[these_donors_idx])
    best_donor_long_id <- output_df$donor_long_id[these_donors_idx][wm]
    output_df$best_match[output_df$donor_long_id == best_donor_long_id] <- TRUE
    output_df <- compute_confidence_score(output_df, these_donors, short_names)
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
message("lines used: ", opt$donor_lines, "\n")

## Run main function
main(opt$input_file, opt$output_prefix, opt$donor_lines, opt$donor_vcf,
     opt$fasta_idx)

## # params for testing
## opt <-  list()
## opt[["input_file"]] <- "data_raw/scrnaseq/run_21999/vcf/21999_1#56.filtered.vcf.gz"
## opt[["output_prefix"]] <- "data_raw/scrnaseq/run_21999/donor_id/tmp.donor_id"
## opt[["donor_vcf"]] <- "data_raw/scrnaseq/run_21999/vcf/21999_1#56.filtered.hipsci.overlap.vcf.gz"
## opt[["donor_lines"]] <- "fasu_2;kegd_2;zerv_8;zoio_2;xojn_3;fuai_1;eevy_7;oaqd_3;paab_4;sita_1;toss_3;zoio_2;heth_1;jogf_2;pelm_3;vass_1;wibj_2;zapk_3"


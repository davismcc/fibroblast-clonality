"
Identify donor from a cell's filtered VCF file

Usage:
identify_donor.R --input_file <DATA_IN> --output_prefix <PRFX_OUT> --donor_lines <DONORS> [--donor_vcf <DONOR_VCF> --fasta_idx <FASTA_IDX> -hv]

Options:
-i --input_file <DATA_IN>      input file
-o --output_prefix <PRFX_OUT>  prefix for output files; .RData, .feather and .csv files produced
-l --donor_lines <DONORS>      string providing donor line IDs separated by ';' (e.g. 'babk_2;wuye_2;podx3')
-d --donor_vcf <DONOR_VCF>     VCF file containing the genotype data for all donors [default: /nfs/research2/stegle/projects/hipsci/data/genotypes/imputed/REL-2014-11_SS/hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.20151005.genotypes.gdid.mac1.recode.WGS.maf0.1.vcf.gz]
-f --fasta_idx <FASTA_IDX>     FASTA index file for FASTA file that was used for read alignment [default: /nfs/research2/stegle/datasets/references/human/STAR_GRCh37.75_ERCC/GRCh37.p13.genome.ERCC92.fa.fai]
-h --help                      show this
-v --version                   print version and stop

This program does compares called genotypes from RNA-seq reads for a cell to
genotype data for all HipSci donors to find the best match in terms of
(weighted) average genotypic correlation across called variants.

The program returns results as a .csv file.

This program requires the R packages 'VariantAnnotation' and 'snpStats'
(Bioconductor).

Davis McCarthy
July 2016
" -> doc

## Script to identify HIPSCI donor for a cell's VCF file
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(snpStats))

## Define main function
main <- function(input_file, output_prefix, donor_lines, donor_vcf, fasta_idx) {
    ## Define files and other global variables
    hipsci_file <- donor_vcf
    sample_file <- input_file
    seq_lengths <- read.delim(fasta_idx, header = FALSE,
                              stringsAsFactors = FALSE)

    ## define samples in donor VCF
    hdr_hipsci <- scanVcfHeader(hipsci_file)
    ave_allelic_corr <- wght_ave_allelic_corr <-
        rep(0, length(samples(hdr_hipsci)))
    names(ave_allelic_corr) <- names(wght_ave_allelic_corr) <-
        samples(hdr_hipsci)
    ## define donors and names
    short_names <- gsub("HPSI[0-9a-z]+-", "", names(ave_allelic_corr))
    these_donors <- strsplit(donor_lines, ";")[[1]]
    these_donors_idx <- grepl(gsub(";", "|", donor_lines),
                              names(ave_allelic_corr))

    ## define output data frame
    output_df <- data.frame(
        sample_id = gsub(".filtered.vcf", "", basename(sample_file)),
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

    ## Read in VCF from this sample
    cat("Reading sample VCF\n")
    vcf_sample <- readVcf(sample_file, "GRCh37")
    keep_seqs <- seqnames(vcf_sample) %in% paste0("chr", c(1:22))
    if (sum(keep_seqs) < 1) {
        write.csv(output_df, file = paste0(output_prefix, ".csv"),
                  row.names = FALSE)
        cat("No variants called on autosomes in sample VCF\n")
        return("Done.")
    }
    vcf_sample_filt <- vcf_sample[keep_seqs]
    vcf_sample_filt <- vcf_sample_filt[isSNV(vcf_sample_filt)]
    if (length(vcf_sample_filt) < 1) {
        write.csv(output_df, file = paste0(output_prefix, ".csv"),
                  row.names = FALSE)
        cat("No variants in sample VCF after filtering.\n")
        return("Done.")
    }
    vcf_sample_filt <- renameSeqlevels(
        vcf_sample_filt, gsub("chr", "", seqlevels(vcf_sample_filt)))
    new_snp_names <- paste0("snp_",
                            gsub("chr", "",
                                 gsub(":", "_",
                                      gsub("_[ATCG]/[ATCG]", "",
                                           names(vcf_sample_filt)))))
    
    names(vcf_sample_filt) <- new_snp_names
    output_df$nvars_called <- length(vcf_sample_filt)
    cat("...read ", length(vcf_sample_filt), "variants from cell VCF\n")

    ## Read in HipSci VCF
    cat("Reading HipSci VCF\n")
    rr <- GRanges(seqnames = seqnames(vcf_sample_filt),
                  ranges = IRanges(
                      start = start(vcf_sample_filt),
                      end = end(vcf_sample_filt),
                      names = names(vcf_sample_filt)))
    ## for speed, sample 9999 if there are 10k or more variants
    if (length(rr) > 9999) {
        idx <- sample.int(length(rr), size = 9999)
        idx <- sort(idx)
        rr <- rr[idx]
        vcf_sample_filt <- vcf_sample_filt[idx]
    }
    cat("...reading ", length(rr), "variants from HipSci VCF\n")
    tab_hipsci <- TabixFile(hipsci_file)
    vcf_hipsci <- readVcf(tab_hipsci, "GRCh37", param = rr)
    isSNV_idx <- isSNV(vcf_hipsci)
    if (!any(isSNV_idx)) {
        write.csv(output_df, file = paste0(output_prefix, ".csv"),
                  row.names = FALSE)
        cat("No single-nucleotide variants overlapping in sample VCF and HipSci VCF\n")
        return("Done.")
    }
    ## filter donor VCF
    cat("Filter donor VCF\n")
    vcf_hipsci_filt <- vcf_hipsci[isSNV_idx]
    sm_hipsci <- genotypeToSnpMatrix(
        geno(vcf_hipsci_filt, "GT"), ref = ref(vcf_hipsci_filt),
        alt = alt(vcf_hipsci_filt))
    donor <- gsub("_[0-9]+", "", short_names)
    donor <- gsub("2:", "",
                  gsub("CTRL0914es-", "",
                       gsub("CTRL1114es-", "", donor)))
    lines_for_stats <- !duplicated(donor)
    snp_stats_hipsci <- col.summary(sm_hipsci$genotypes[lines_for_stats,])
    common_vars <- snp_stats_hipsci$MAF > 0.05
    common_vars[is.na(common_vars)] <- FALSE
    if (sum(common_vars, na.rm = TRUE) == 0) {
        write.csv(output_df, file = paste0(output_prefix, ".csv"),
                  row.names = FALSE)
        cat("No common variants overlapping in sample VCF and HipSci VCF\n")
        return("Done.")
    }
    vcf_hipsci_filt2 <- vcf_hipsci_filt[common_vars]
    sm_hipsci <- genotypeToSnpMatrix(
        geno(vcf_hipsci_filt2, "GT"), ref = ref(vcf_hipsci_filt2),
        alt = alt(vcf_hipsci_filt2))
    snp_stats_hipsci <- col.summary(sm_hipsci$genotypes[lines_for_stats,])
    donor_stats_hipsci <- row.summary(sm_hipsci$genotypes)

    ## filter sample VCF to those variants found in donor VCF
    cat("Filtering sample VCF\n")
    slengths_sample <- seqlengths(vcf_sample_filt)[1:25]
    names(slengths_sample)[25] <- "MT"
    seqlengths(vcf_hipsci_filt2) <- slengths_sample[seqlevels(vcf_hipsci_filt2)]
    ovlap <- findOverlaps(vcf_sample_filt, vcf_hipsci_filt2)
    if (length(ovlap) < 1L) {
        write.csv(output_df, file = paste0(output_prefix, ".csv"),
                  row.names = FALSE)
        cat("No common variants overlapping in sample VCF and HipSci VCF\n")
        return("Done.")
    }
    vcf_sample_filt2 <- vcf_sample_filt[queryHits(ovlap)]
    sm_sample <- genotypeToSnpMatrix(geno(vcf_sample_filt2, "GT"),
                                     ref = ref(vcf_sample_filt2),
                                     alt = alt(vcf_sample_filt2))
    sm_sample$map$quality <- qual(vcf_sample_filt2)
    sm_sample$map$weight <- qual(vcf_sample_filt2) / sum(qual(vcf_sample_filt2))
    if ( !identical(as.character(sm_sample$map$allele.1),
              as.character(sm_hipsci$map$allele.1)))
        stop("Alleles do not match between sample and donor VCFs")
    cat("Computing relatedness using ", length(vcf_sample_filt2), "variants\n")
    
    ## standardise sample genotypes
    cat("Standardise genotypes\n")
    hipsci_geno_mat <- as(sm_hipsci$genotype, "numeric")
    hipsci_geno_mat[is.na(hipsci_geno_mat)] <- 0
    ## standardise sample genotypes
    sample_geno_vec <- as(sm_sample$genotype, "numeric")[1,]
    sm_sample_stand_geno <- ( (sample_geno_vec - 2 * snp_stats_hipsci$RAF)
        / sqrt(2 * snp_stats_hipsci$RAF * (1 - snp_stats_hipsci$RAF)) )
        ## standardise donor genotypes
    sm_hipsci_stand_geno <- t( (t(hipsci_geno_mat) - 2 * snp_stats_hipsci$RAF)
    / sqrt(2 * snp_stats_hipsci$RAF * (1 - snp_stats_hipsci$RAF)))

    ## snp_means <- colMeans(hipsci_geno_mat)
    ## snp_vars <- matrixStats::colVars(hipsci_geno_mat)
    ## sm_sample_stand_geno <- ( (as(sm_sample$genotype, "numeric") - snp_means)
    ##     / sqrt(snp_vars))
    ## sm_sample_stand_geno <- sm_sample_stand_geno[1,]
    ## sm_hipsci_stand_geno <- t( (t(hipsci_geno_mat) - snp_means) / sqrt(snp_vars))

    ## compute weighted correlation between genotypes of each hipsci donor and
    ## the sample
    cat("Compute genotypic correlations\n")
    corr_mat <- t(t(sm_hipsci_stand_geno) * sm_sample_stand_geno)
    wght_corr_mat <- t(t(sm_hipsci_stand_geno) * sm_sample_stand_geno *
                       sm_sample$map$weight)
    bad_snp <- apply(corr_mat, 2, function(x) any(is.na(x)))
    corr_mat <- corr_mat[, !bad_snp, drop = FALSE]
    wght_corr_mat <- wght_corr_mat[, !bad_snp, drop = FALSE]
    ## sum correlations across donors
    ave_allelic_corr <- rowMeans(corr_mat, na.rm = TRUE)
    wght_ave_allelic_corr <- rowSums(wght_corr_mat, na.rm = TRUE)

    ## use ecdf of correlation values to get confidence for called donor
    ## as percentile of distribution of scores across all HipSci donors
    f_ave_corr <- ecdf(ave_allelic_corr)
    f_wght_ave_corr <- ecdf(wght_ave_allelic_corr)

    ## define output data frame
    output_df <- data.frame(
        sample_id = gsub("_filtered.vcf", "", basename(sample_file)),
        donor_long_id = as.character(names(ave_allelic_corr)),
        donor_short_id = short_names,
        nvars_called = length(vcf_sample),
        nvars_ovlap_hipsci = length(vcf_hipsci),
        nvars_used = length(vcf_sample_filt2),
        ave_allelic_corr = ave_allelic_corr,
        wght_ave_allelic_corr = wght_ave_allelic_corr,
        pctile_ave_corr = f_ave_corr(ave_allelic_corr),
        pctile_wght_ave_corr = f_wght_ave_corr(wght_ave_allelic_corr),
        used_in_expt = (short_names %in% these_donors),
        best_match = FALSE, stringsAsFactors = FALSE)

    ## define best match from those donors known to have been used
    wm <- which.max(ave_allelic_corr[these_donors_idx])
    best_donor_long_id <- output_df$donor_long_id[these_donors_idx][wm]
    output_df$best_match[output_df$donor_long_id == best_donor_long_id] <- TRUE
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
    ## write output to file
    cat("Writing output to CSV\n")
    write.csv(output_df, file = paste0(output_prefix, ".csv"),
              row.names = FALSE)
    return("Done.")
}

## Get command line options
opt <- docopt::docopt(doc, version = "version 0.0.1\n")

cat("working directory: ", getwd(), "\n")
cat("input vcf: ", opt$input_file, "\n")
cat("output prefix: ", opt$output_prefix, "\n")
cat("donor vcf: ", opt$donor_vcf, "\n")
cat("fasta idx: ", opt$fasta_idx, "\n")
cat("lines used: ", opt$donor_lines, "\n")

## Run main function
main(opt$input_file, opt$output_prefix, opt$donor_lines, opt$donor_vcf,
     opt$fasta_idx)

## params for testing
## opt <-  list()
## opt[["input_file"]] <- "data_raw/scrnaseq/run_21672/vcf/21672_4#223.filtered.vcf"
## opt[["output_prefix"]] <- "data_raw/scrnaseq/run_21672/donor_id/21672_4#223.donor_id"
## opt[["donor_vcf"]] <- "/hps/nobackup/stegle/projects/hipsci/hipsci_genotypes/unreleased/hipsci.wec.gtarray.unreleased.REL-2016-09.imputed_phased.INFO_0.4_filtered.20160912.genotypes.allchr.temp_for_endodiff_20170125.vcf.gz"
## opt[["fasta_idx"]] <- "/hps/nobackup/stegle/datasets/references/human/STAR_GRCh37.75_ERCC/GRCh37.p13.genome.ERCC92.fa.fai"
## opt[["donor_lines"]] <- "fasu_2;kegd_2;zerv_8;zoio_2;xojn_3;fuai_1;eevy_7;oaqd_3;paab_4;sita_1;toss_3;zoio_2;heth_1;jogf_2;pelm_3;vass_1;wibj_2;zapk_3"


"
Filter somatic variants from a flat file prepared by Petr Danecek

Usage:
filter_variants.R --input_file <DATA_IN> --output_file <FILE_OUT> [--donor_cell_vcf <VCF_IN> --max_fdr <fdr> --min_vaf_fibro <vaf> --max_vaf_fibro <vaf> --min_nalt_fibro <num> --max_vaf_ips <vaf> --combo_max_vaf_fibro <vaf> --combo_max_vaf_ips <vaf> --min_prop_covered_cells <prop> --donor_name <short_id> -hv]

Options:
-i --input_file <DATA_IN>        input file
-o --output_file <FILE_OUT>      output files
--donor_cell_vcf <VCF_IN>        (optional) VCF file (compressed or not) with sites genotyped in cells for a donor
--max_fdr <fdr>                  upper threshold for False Discovery Rate [default: 0.1]
--min_vaf_fibro <vaf>            minimum VAF in fibroblast sample [default: 0.01]
--max_vaf_fibro <vaf>            maximum VAF in fibroblast sample [default: 0.6]
--min_nalt_fibro <num>           minimum num ALT alleles in fibro [default: 2.5]
--max_vaf_ips <vaf>              maximum VAF in iPS sample [default: 0.7]
--combo_max_vaf_fibro <vaf>      maximum VAF in fibroblast sample to be used in conjunction with combo_max_vaf_ips [default: 0.35]
--combo_max_vaf_ips <vaf>        maximum VAF in fibroblast sample to be used in conjunction with combo_max_vaf_ips [default: 0.3]
--min_prop_covered_cells <prop>  minimum proportion of cells with coverage for site to be retained [default: 0.1]
--donor_name <short_id>          donor name (short HipSci ID) [default: NULL]
-h --help                        show this
-v --version                     print version and stop

The program returns filtered variants from a flat table file produced by Petr
Danecek (Sanger) and results as a tab-delimited table file (can be compressed
automatically based on file extension provided.)

This program requires the R packages 'tidyverse'.

Petr Danecek prepared the exome dataset as follows. For each donor:

* detect all sites where we see a non-reference base in the pileup;
* require minimum per-sample coverage of 20 reads;
* require at least 3 alternate reads in either fibroblast or iPS;
* require allele frequency in ExAC and 1000GP < 5%.

Fishers Exact Test is then performed on these variants to detect
fibroblast-iPS pairs with significantly changed proportion of REF/ALT
reads.

On output a tab-delimited text file is produced which includes fields
shown below (*) and which should be further filtered, depending on the
purpose of the analysis.

For clonal inference, we apply the following filters by default:

* a Benjamini-Hochberg FDR of less than 10%;
* a minimum VAF in high-coverage fibroblast of 1%;
* a maximum VAF in high-coverage fibroblast of 60% (to avoid the rare
  possibility of a homozygous alternative mutation);
* a maximum VAF in low-coverage iPSCs of 70% (to avoid the unlikely possibility
  that a homozygous alternative mutation in the fibroblasts corresponds to what
  must be a homozygous mutation in some cells in iPSCs);
* at least 3 alternative alleles observed for the site;
* a VAF < combo_max_vaf_fibro in the fibroblast sample *or* VAF < combo_max_vaf_ips in the iPS sample
* require uniqueness of sites across donors as it is highly unlikely to observe
  the same point mutations and they are most likely artefacts of some sort.

Optionally, sites can be filtered based on an input VCF such that only sites 
with non-zero coverage in at least 10% of cells are retained.

Davis McCarthy
April 2018
" -> doc


## Script to identify HIPSCI donor for a cell's VCF file
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(VariantAnnotation))

main <- function(input_file, output_file, vcf_file = NULL, max_fdr = 0.1,
                min_vaf_fibro = 0.01, max_vaf_fibro = 0.6, min_nalt_fibro = 2.5,
                max_vaf_ips = 0.7, combo_max_vaf_fibro = 0.35,
                combo_max_vaf_ips = 0.3, min_prop_covered_cells = 0.1, 
                donor_name = NULL) {
    max_fdr <- as.numeric(max_fdr)
    message("FDR threshold: ", max_fdr)
    min_vaf_fibro <- as.numeric(min_vaf_fibro)
    message("min VAF in fibro: ", min_vaf_fibro)
    max_vaf_fibro <- as.numeric(max_vaf_fibro)
    message("max VAF in fibro: ", max_vaf_fibro)
    min_nalt_fibro <- as.numeric(min_nalt_fibro)
    message("min num ALT counts in fibro: ", min_nalt_fibro)
    max_vaf_ips <- as.numeric(max_vaf_ips)
    exome_sites_hicov <- read_tsv(
    input_file,
    col_types = "ciccdcciiiicccccccc", comment = "#",
    col_names = c("chrom", "pos", "fibro", "ips", "Fisher_Pval", "ref", "alt",
                  "nREF_fibro", "nALT_fibro", "nREF_ips", "nALT_ips",
                  "AF_1000G", "AF_ExAC", "clinsing", "SIFT_Pdel",
                  "PolyPhen_Pdel", "gene_name", "ENSID", "bcftools_csq"))
    exome_sites_hicov <- dplyr::mutate(
            exome_sites_hicov,
            propALT_fibro = (nALT_fibro / (nREF_fibro + nALT_fibro)),
            propALT_ips = (nALT_ips / (nREF_ips + nALT_ips)),
            donor_short_id = gsub("HPS.*-", "", fibro),
            fdr = p.adjust(Fisher_Pval, method = "BH"),
            var_id = paste0(chrom, ":", pos, "_", ref, "_", alt))
    ## identify non-unique sites
    nonuniq_sites <- exome_sites_hicov %>% group_by(var_id) %>%
        summarise(n_donors = length(unique(fibro)))
    exome_sites_hicov <- dplyr::mutate(exome_sites_hicov, 
        uniq_site = (var_id %in% nonuniq_sites[["var_id"]][
            nonuniq_sites[["n_donors"]] < 1.5]),
        sig_fdr = (fdr < max_fdr & propALT_fibro > min_vaf_fibro & 
            propALT_fibro < max_vaf_fibro & nALT_fibro > min_nalt_fibro & 
            propALT_ips < max_vaf_ips & uniq_site & 
            (propALT_fibro < combo_max_vaf_fibro | propALT_ips < combo_max_vaf_ips)))
    if (!is.null(vcf_file)) {
        if (is.null(donor_name))
            stop("donor_name argument must be supplied with vcf_file option.")
        vcf <- readVcf(vcf_file, "GRCh37")
        ## keep variants with non-zero depth in at least 10% of cells (default)
        suff_dp <- rowMeans(geno(vcf, "DP") > 0.5) > min_prop_covered_cells
        var_id_suff_dp <- gsub(
            "_.*", "", gsub("chr", "", rownames(vcf[suff_dp,])))
        keep_var <- (paste(exome_sites_hicov[["chrom"]], 
                    exome_sites_hicov[["pos"]], sep = ":") %in% var_id_suff_dp)
        right_donor <- exome_sites_hicov[["donor_short_id"]] %in% donor_name
        exome_sites_hicov <- exome_sites_hicov[keep_var & right_donor,]
    }
    write_tsv(dplyr::filter(exome_sites_hicov, uniq_site, sig_fdr),
        output_file, col_names = TRUE)
    return("Done.")
}

## Get command line options
opt <- docopt::docopt(doc, version = "version 0.0.1\n")

message("working directory: ", getwd(), "\n")
message("input file: ", opt$input_file, "\n")
message("output file: ", opt$output_file, "\n")


## Run main function
main(opt$input_file, opt$output_file, opt$donor_cell_vcf, opt$max_fdr,
    opt$min_vaf_fibro, opt$max_vaf_fibro, opt$min_nalt_fibro, opt$max_vaf_ips,
    opt$combo_max_vaf_fibro, opt$combo_max_vaf_ips, opt$min_prop_covered_cells,
    opt$donor_name)

# opts <- list()
# opts$input_file <- "Data/exome-point-mutations/high-vs-low-exomes.v62.ft.txt.gz"
# opts$output_file <- "Data/exome-point-mutations/high-vs-low-exomes.v62.ft.filt_fdr0.1.txt.gz"
# opts$donor_cell_vcf <- "Data/SS2_2017/mpileup/vass.mpileup.vcf.gz"

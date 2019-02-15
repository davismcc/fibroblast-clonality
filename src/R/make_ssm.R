"
Make SSM format data from VAF site data

Usage:
make_ssm.R --input_file <DATA_IN> --output_file <DATA_OUT> [-hv]

Options:
-i --input_file <DATA_IN>      input file
-o --output_file <DATA_OUT>    output file
-h --help                      show this
-v --version                   print version and stop

This program reads in VAF data prepared by Petr Danecek and produces a data file
 in SSM format required for input to phyloWGS.

The program returns results as a .txt (tab-delimited) file.

Columns of data from Petr:
[1] chrom
[2] pos
[3] fibro
[4] ips
[5] ref
[6] alt
[7] Fisher test p-value - fibro vs ips in low-coverage data
[8] Fisher test p-value - fibro low-coverage vs fibro high-coverage
[9] nALT/(nREF+nALT) in low-cov fibro
[10] nALT/(nREF+nALT) in high-cov fibro
[11] low-cov fibro: nREF
[12] low-cov fibro: nALT
[13] low-cov ips:   nREF
[14] low-cov ips:   nALT
[15] high-cov fibro: nREF
[16] high-cov fibro: nALT

Davis McCarthy
March 2017
" -> doc

## Script to identify HIPSCI donor for a cell's VCF file
suppressPackageStartupMessages(library(tidyverse))

main <- function(input_file, output_file) {
    ## write output to file
    vaf_data <- read_tsv(input_file, col_names = FALSE, col_types = cols())
    colnames(vaf_data) <- c("chrom", "pos", "fibro", "ips", "ref", "alt",
                            "pval_fibro_ips_hicov", "pval_fibro_hi_fibro_lo",
                            "prop_alt_fibro_lo", "prop_alt_fibro_hi",
                            "nREF_fibro_lo", "nALT_fibro_lo", "nREF_ips_lo",
                            "nALT_ips_lo", "nREF_fibro_hi", "nALT_fibro_hi") 
    ssm_data <- data_frame(id = paste0("s", 0:(nrow(vaf_data) - 1)),
                           gene = paste0(vaf_data[["chrom"]], "_",
                                         vaf_data[["pos"]]),
                           a = vaf_data[["nREF_fibro_hi"]],
                           d = vaf_data[["nREF_fibro_hi"]] + vaf_data[["nALT_fibro_hi"]],
                           mu_r = 0.999, mu_v = 0.499)
    keep_site <- (vaf_data[["prop_alt_fibro_hi"]] > 0.1 &
                  vaf_data[["prop_alt_fibro_hi"]] < 0.9 &
                  ssm_data[["d"]] >= 100)
    message("Writing output to TSV\n")
    ssm_data <- ssm_data[keep_site,]
    ssm_data[["id"]] <- paste0("s", 0:(nrow(ssm_data) - 1))
    write_tsv(ssm_data, output_file)
    return("Done.")
}

## Get command line options
opt <- docopt::docopt(doc, version = "version 0.0.1\n")

message("working directory: ", getwd(), "\n")
message("input file: ", opt$input_file, "\n")
message("output file: ", opt$output_prefix, "\n")

## Run main function
main(opt$input_file, opt$output_file)

## # params for testing
## opt <-  list()
## opt[["input_file"]] <- "Data/exome-point-mutations/vass/highcov.2017-05-26.vass.txt"
## opt[["output_file"]] <- "Data/exome-point-mutations/vass/ssm_highcov.2017-05-26.vass.txt"
## input_file <- "Data/exome-point-mutations/vass/highcov.2017-05-26.vass.txt"
## output_file <- "Data/exome-point-mutations/vass/ssm_highcov.2017-05-26.vass.txt"
## vaf_data <- read_tsv(input_file, col_names = FALSE, col_types = cols())
## colnames(vaf_data) <- c("chrom", "pos", "fibro", "ips", "ref", "alt",
##                         "pval_fibro_ips_hicov", "pval_fibro_hi_fibro_lo",
##                         "prop_alt_fibro_lo", "prop_alt_fibro_hi",
##                         "nREF_fibro_lo", "nALT_fibro_lo", "nREF_ips_lo",
##                         "nALT_ips_lo", "nREF_fibro_hi", "nALT_fibro_hi") 
## ssm_data <- data_frame(id = paste0("s", 0:(nrow(vaf_data) - 1)),
##                        gene = paste0(vaf_data[["chrom"]], "_", vaf_data[["pos"]]),
##                        a = vaf_data[["nREF_fibro_hi"]],
##                        d = vaf_data[["nREF_fibro_hi"]] + vaf_data[["nALT_fibro_hi"]],
##                        mu_r = 0.999, mu_v = 0.499)
## keep_site <- (vaf_data[["prop_alt_fibro_hi"]] > 0.1 & ssm_data[["d"]] >= 100)
## write_tsv(ssm_data[keep_site,], output_file)

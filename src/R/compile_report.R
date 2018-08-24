'
Compile and automated report from an .RData file with an existing template

Usage:
compile_report.R --input_file <DATA_IN> --output_file <HTML_OUT> [--results_out <RDS_OUT> --template <TEMPLATE_FILE> --title <TITLE> --author <AUTHOR> --donor <DONOR> --to_working_dir <WORK_DIR> -hv]

Options:
-i --input_file <DATA_IN>        input file
-o --output_file <HTML_OUT>      specify output HTML file
--results_out <RDS_OUT>          specify file for .rds output file [default: NA]
-t --template <TEMPLATE_FILE>    name of template file to use for report [default: ../Rmd/rough_qc_template.Rmd]
--title <TITLE>                  report title [default: Rough QC Report]
--author <AUTHOR>                report author [default: Davis McCarthy]
--donor <DONOR>                  HipSci donor cell line [default: MISSING]
-w --to_working_dir <WORK_DIR>   set the working directory for the Rmd template file [default: ./]
-h --help                        show this
-v --version                     print version and stop

This program does two distinct things compiles an html report from an RMarkdown 
template and given a dataset and other parameters.

This program requires the R packages "docopt" (CRAN) and "scater" (github.com/davismcc/scater).

Davis McCarthy
November 2017
' -> doc

## Define main function
main <- function(input_file, output_file, results_out, template, 
                 title = "Rough QC Report", author = "Davis McCarthy",
                 donor = "MISSING", to_working_dir = "./") {
    rmarkdown::render(input = template, output_file = basename(output_file),
                      output_dir = dirname(output_file),
                      params = list(input_file = input_file, 
                                    output_file = results_out,
                                    author = author, title = title,
                                    donor = donor,
                                    to_working_dir = to_working_dir), 
                      output_format = "html_document")
}

## Get command line options
opt <- docopt::docopt(doc, version = "version 0.0.1\n")

## Run main function
main(opt$input_file, opt$output_file, opt$results_out, opt$template, opt$title, 
    opt$author, opt$donor, opt$to_working_dir)

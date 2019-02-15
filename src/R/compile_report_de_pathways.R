'
Compile and automated report from an .RData file with an existing template

Usage:
compile_report_cell_assign.R --callset <CALLSET_ID> --output_file <HTML_OUT> [--template <TEMPLATE_FILE> --title <TITLE> --author <AUTHOR> --to_working_dir <WORK_DIR> -hv]

Options:
-c --callset <CAllSET_ID>      input .rds file containing a SingleCellExperiment object
-o --output_file <HTML_OUT>      specify output HTML file
-t --template <TEMPLATE_FILE>  name of template file to use for report [default: ../Rmd/DE_pathways_FTv62_callset_clones_pairwise_vs_base.Rmd]
--title <TITLE>                  report title [default: Rough QC Report]
--author <AUTHOR>                report author [default: Davis McCarthy]
-w --to_working_dir <WORK_DIR>   set the working directory for the Rmd template file [default: ./]
-h --help                        show this
-v --version                     print version and stop

This program does two distinct things compiles an html report from an RMarkdown 
template and given a dataset and other parameters.

This program requires the R packages "docopt" (CRAN) and "scater" (github.com/davismcc/scater).

Davis McCarthy
May 2018
' -> doc

## Define main function
main <- function(callset, output_file, template, 
                 title = "DE Pathway Analysis", author = "Davis McCarthy",
                 to_working_dir = "./") {
    rmarkdown::render(input = template, output_file = basename(output_file),
                      output_dir = dirname(output_file),
                      params = list(callset = callset, 
                                    author = author, 
                                    title = title,
                                    to_working_dir = to_working_dir), 
                      output_format = "html_document")
}

## Get command line options
opt <- docopt::docopt(doc, version = "version 0.0.1\n")

cat("callset: ", opt$callset, "\n")
cat("to_working_dir: ", opt$to_working_dir, "\n")
cat("working_dir: ", getwd(), "\n")

## Run main function
main(opt$callset, opt$output_file, opt$template, opt$title, opt$author, 
     opt$to_working_dir)

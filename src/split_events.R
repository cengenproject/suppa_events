#!/bin/env Rscript

# Split psi file by neuron to run diffSplice

# See this to split with a list of arbitrary conditions: https://github.com/comprna/SUPPA/blob/master/scripts/split_file.R


library(getopt)

# Read command line arguments ----
if(! interactive()){
  spec <- matrix(c(
    'input_path',  'i', 1, "character", "Path to PSI file",
    'output_path', 'o', 1, "character", "Path to directory in which to save split PSI",
    'help',        'h', 0, "logical",   "Print this help"
  ), byrow=TRUE, ncol=5)
  
  opt <- getopt(spec)
} else{
  # Options for experimenting
  opt <- list(
    input_path = "data/240301b_psiPerEvent.psi",
    output_path = "data/240301_psi_condition"
  )
}


if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}





psi <- read.delim(input_path)


samples <- names(psi)
neurons <- stringr::str_match(samples, "^([A-Zef0-9]{2,4})r[0-9]{1,4}")[,2]

stopifnot(all(! is.na(neurons) ))

samples_by_neuron <- split(samples, neurons)

purrr::iwalk(samples_by_neuron,
            ~ {
              write.table(psi[, .x],
                          file = file.path(output_path, paste0(.y, ".psi")))
              })




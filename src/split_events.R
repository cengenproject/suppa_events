#!/bin/env Rscript

# Split psi/tpm file by neuron to run diffSplice

# See this to split with a list of arbitrary conditions: https://github.com/comprna/SUPPA/blob/master/scripts/split_file.R


library(getopt)

# Read command line arguments ----
if(! interactive()){
  spec <- matrix(c(
    'input_path',  'i', 1, "character", "Path to PSI file",
    'output_path', 'o', 1, "character", "Path to directory in which to save split PSI",
    'extension',   'e', 1, "character", "File extension for the output (e.g. psi or tpm)",
    'help',        'h', 0, "logical",   "Print this help"
  ), byrow=TRUE, ncol=5)
  
  opt <- getopt(spec)
} else{
  # Options for experimenting
  opt <- list(
    input_path = "data/240301b_psiPerEvent.psi",
    output_path = "data/240301_psi_condition",
    extension = "psi"
  )
  
  # opt <- list(
  #   input_path = "data/231208_str_q_tx_TPM.tsv",
  #   output_path = "data/240301_psi_condition",
  #   extension = "tpm"
  # )
}


if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}





psi <- read.delim(opt$input_path)


samples <- names(psi)
neurons <- stringr::str_match(samples, "^([A-Zef0-9]{2,4})r[0-9]{1,4}")[,2]

stopifnot(all(! is.na(neurons) ))

samples_by_neuron <- split(samples, neurons)

samples_by_neuron <- samples_by_neuron[lengths(samples_by_neuron) > 1]


purrr::iwalk(samples_by_neuron,
            ~ {
              write.table(psi[, .x],
                          file = file.path(opt$output_path, paste0(.y, ".", opt$extension)),
                          quote = FALSE)
              })




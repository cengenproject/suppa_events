

outliers_file <- "/gpfs/gibbs/pi/hammarlund/CeNGEN/bulk/bulk_alignments/bsn12_outliers_to_ignore.txt"

str_q_input_file <- "data/231208_str_q_tx_TPM.tsv"
str_q_output_file <- "data/240913_str_q_tx_TPM_filt.tsv"


outliers <- readLines(outliers_file) |>
  setdiff("")


dat <- read.table(str_q_input_file)

dat <- dat[,-which(colnames(dat) %in% outliers)]

write.table(dat, str_q_output_file, quote = FALSE, sep = "\t")


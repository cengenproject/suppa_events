

# A look at PSI across neuron types


library(tidyverse)


library(wbData)

tx2g <- wb_load_tx2gene(289)
gids <- wb_load_gene_ids(289)



# Load ----

psi <- read.delim("data/240301b_psiPerEvent.psi") |>
  rownames_to_column("event_id") |>
  as_tibble() |>
  separate_wider_regex(event_id,
                       patterns = c(gene_id = "^WBGene[0-9]{8}", ";",
                                    event_type = "[SEA53MXRIFL]{2}", "\\:",
                                    event_coordinates = "[IXV]+\\:[0-9:\\-]+:[+-]$"),
                       cols_remove = FALSE)

psi_lg <- pivot_longer(psi,
                       -c(gene_id, event_type, event_coordinates, event_id),
                       names_to = "sample_id",
                       values_to = "PSI") |>
  mutate(neuron_id = str_match(sample_id, "^([A-Zef0-9]{2,4})r[0-9]{1,4}")[,2])



tpm <- read.delim("data/231208_str_q_tx_TPM.tsv") |>
  rownames_to_column("transcript_id") |>
  as_tibble() |>
  mutate(gene_id = wb_tx2g(transcript_id, tx2g, warn_missing = TRUE),
         .after = "transcript_id")






# Filter ----

# Older: filter genes that have enough tpm in any sample. See below, only looking within neuron types
# sum_tpm <- tibble(gene_id = tpm_gene$gene_id,
#                   sum_tpm = rowSums(tpm_gene |> select(-transcript_id, -gene_id))) |>
#   summarize(tpm = sum(sum_tpm),
#             .by = "gene_id") |>
#   mutate(in_neurons = gene_id %in% wormDatasets::genes_by_pattern$present_in_neurons)
# 
# hist(log10(sum_tpm$tpm), breaks = 100)
# 
# ggplot(sum_tpm) +
#   theme_classic() +
#   geom_density(aes(x = tpm, fill = in_neurons), alpha = .5) +
#   scale_x_log10() +
#   geom_vline(aes(xintercept = 150), linetype = 'dashed')
# 
# genes_expressed_here <- sum_tpm$gene_id[sum_tpm$tpm > 150]
# 
# length(genes_expressed_here)
# 
# 
# # Note: 14,601 protein-coding, the rest are other biotypes
# # xx <- gids |> mutate(here = gene_id %in% genes_expressed_here)
# # 
# # table(xx$biotype, xx$here)
# 
# 
# 
# filt_psi <- psi |>
#   filter(gene_id %in% genes_expressed_here)
# 
# filt_psi_lg <- psi_lg |>
#   filter(gene_id %in% genes_expressed_here)


# filter by neuron (replace with Alec's data when available)

tpm_gene_neur <- tpm |>
  pivot_longer(-c(transcript_id, gene_id),
               names_to = "sample_id",
               values_to = "TPM") |>
  mutate(neuron_id = str_match(sample_id, "^([A-Zef0-9]{2,4})r[0-9]{1,4}")[,2]) |>
  summarize(TPM = sum(TPM),
            .by = c("gene_id", "neuron_id"))
  
tpm_gene_neur |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = TPM), bins = 200) +
  scale_x_log10() +
  geom_vline(aes(xintercept = 15), linetype = 'dashed')

# a few example (see on browser)
tpm_gene_neur |>
  filter(gene_id == "WBGene00012697") |>
  arrange(desc(TPM)) |> View()

tpm_gene_neur <- tpm_gene_neur |>
  mutate(is_expressed = TPM > 15)



# double check: relatively few cases where gene judged not expressed but gene-neuron judged expressed
# xx <- tpm_gene_neur |>
#   mutate(gene_here = gene_id %in% genes_expressed_here)
# 
# table(xx$gene_here, xx$is_expressed)


psi_lg <- psi_lg |>
  left_join(tpm_gene_neur,
            by = c("gene_id", "neuron_id"))

filt_psi_lg <- psi_lg |> filter(is_expressed)




# event types

table(psi$event_type)


# single event




my_ev <- sample(filt_psi_lg$event_id, 1)

my_ev
filt_psi_lg |>
  filter(event_id == my_ev) |>
  ggplot(aes(x = neuron_id, y = PSI)) +
  theme_classic() +
  geom_violin(fill = 'grey90', color = 'grey90', width = 2) +
  # ggbeeswarm::geom_quasirandom() +
  # geom_point(alpha = .5) +
  geom_jitter(aes(x = neuron_id, y = PSI), height = 0, width = .2, alpha = .5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle(my_ev)

# ex: my_ev <- "WBGene00010867;SE:V:13370636-13370772:13370852-13370948:-"
# => need to filter individual neurons






# Dominant spliceform ----






























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


# filter by neuron
gene_expression <- read.delim("../majiq/data/2024-03-05_alec_integration/bsn12_subtracted_integrated_binarized_expression_withVDDD_FDR0.05_030424.tsv")

neurs_sequenced <- colnames(gene_expression)


filt_psi_lg <- psi_lg |>
  filter(neuron_id %in% neurs_sequenced) |>
  left_join(gene_expression |>
              as.data.frame() |>
              rownames_to_column("gene_id") |>
              pivot_longer(-gene_id,
                           names_to = "neuron_id",
                           values_to = "expressed")) |>
  filter(expressed == 1L)




# Specificity ----



psi_by_neuron <- filt_psi_lg |>
  filter(!is.na(PSI)) |>
  summarize(PSI_neuron = mean(PSI, na.rm = TRUE),
            nb_samples = n(),
            .by = c("event_id", "gene_id","event_type",
                    "neuron_id"))

hist(psi_by_neuron$PSI_neuron, breaks = 100)

ggplot(psi_by_neuron) +
  theme_classic() +
  geom_density(aes(x = PSI_neuron, fill = event_type), alpha = .3)





min_pairwise_diff <- function(x){
  
  x <- x[! is.na(x)]
  n <- length(x)
  
  if(n < 2) return(0)
  
  diffs <- numeric(length(x))
  for(neur in seq_along(x)){
    diffs[[neur]] <- DescTools::CombPairs(x[neur], x[-neur]) |>
      as.matrix() |>
      matrixStats::rowDiffs() |>
      abs() |>
      min()
  }
  
  max(diffs)
}




psi_var <- psi_by_neuron |>
  filter(!is.na(PSI_neuron),
         nb_samples > 2) |>
  summarize(specificity = min_pairwise_diff(PSI_neuron) / sd(PSI_neuron, na.rm = TRUE),
            sd = sd(PSI_neuron, na.rm = TRUE),
            nb_neurons = n(),
            .by = c("event_id", "gene_id","event_type"))

psi_var |>
  # filter(n > 10) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x = sd, y = specificity),
                               alpha = .5) +
  ylab("Specificity index") + xlab(NULL) +
  geom_hline(aes(yintercept = 1), linetype = 'dashed', color = 'grey')



psi_var |>
  filter(nb_neurons > 5, sd > .1) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x = event_type, y = specificity),
                               alpha = .5) +
  ylab("Specificity index") + xlab(NULL) +
  geom_hline(aes(yintercept = 1), linetype = 'dashed', color = 'grey')




# ggsave("variance_specificity.pdf", path = export_dir,
#        width = 21, height = 12, units = "cm")


#~ Stat tests ----

filt_psi_var <- psi_var |>
  filter(nb_neurons > 5, sd > .1)

kruskal.test(filt_psi_var$specificity,
             g = filt_psi_var$event_type)

conover.test::conover.test(filt_psi_var$specificity,
                           g = filt_psi_var$event_type,
                           method = "holm")



krusk_var <- agricolae::kruskal(filt_psi_var$specificity,
                                filt_psi_var$event_type,
                                console = TRUE)



filt_psi_var |>
  ggplot(aes(x = event_type)) +
  theme_classic() +
  geom_violin(aes(y = specificity),
              fill = "grey90"
  ) +
  stat_summary(aes(y = specificity),
               fun.min = \(.x) quantile(.x, .25),
               fun.max = \(.x) quantile(.x, .75),
               fun = median,
               color = "red3") +
  geom_text(data = krusk_var$groups |> rownames_to_column("event_type"),
            aes(label = groups),
            y = 4, nudge_x = .1) +
  ylab("Variance") + xlab(NULL)




# ggsave("var_specificity_test.pdf", path = export_dir,
#        width = 21, height = 12, units = "cm")









#~ look at most specific events ----

extremal_af <- filt_psi_var |> 
  # filter(event_type == "AF") |>
  slice_max(order_by = specificity, n = 20)

# ev <- extremal_af$event_id[[1]]



mat_extr_var <- psi_by_neuron |>
  filter(event_id %in% extremal_af$event_id) |>
  pivot_wider(id_cols = event_id,
              names_from = neuron_id,
              values_from = PSI_neuron) |>
  column_to_rownames("event_id") |>
  as.matrix()

annot_df <- psi_by_neuron |>
  filter(event_id %in% extremal_af$event_id) |>
  select(-PSI_neuron, -neuron_id, -nb_samples) |>
  distinct() |>
  left_join(filt_psi_var,
            by = c("event_id", "gene_id", "event_type")) |>
  column_to_rownames("event_id") |>
  select(specificity) |>
  arrange(desc(specificity))

pheatmap::pheatmap(mat_extr_var[rownames(annot_df),],
                   annotation_row = annot_df,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   show_rownames = FALSE,
                   # filename = file.path(export_dir, "heat_variance_af.pdf"),
                   width = 10,
                   height = 7,
                   main = "AF",
                   border_color = NA)
ev <- (filt_psi_var |> arrange(desc(specificity)) |> pull(event_id))[[3]]

filt_psi_lg |>
  filter(event_id == ev) |>
  filter(!is.na(PSI)) |>
  mutate(nb_samples = n(), .by = neuron_id) |>
  filter(nb_samples > 2) |>
  arrange(desc(PSI)) |> mutate(neuron_id = fct_inorder(neuron_id)) |>
  ggplot() +
  theme_classic() +
  geom_violin(aes(x = neuron_id, y = PSI)) +
  ggbeeswarm::geom_quasirandom(aes(x = neuron_id, y = PSI)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#










filt_psi_var |>
  mutate(is_ev = event_id == ev) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = event_id, y = specificity, color = is_ev))


psi_by_neuron |> filter(event_id %in% ev) |> 
  group_by(event_id) |>
  arrange(event_id, desc(PSI_neuron)) |>
  mutate(rank = rank(1-PSI_neuron, ties.method = "first")) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = rank, y = PSI_neuron, color = event_id)) +
  theme(legend.position = "none")












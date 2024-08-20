

# A look at PSI across neuron types


library(tidyverse)


library(wbData)

# tx2g <- wb_load_tx2gene(289)
gids <- wb_load_gene_ids(289)

export_dir <- "data/outs/240813_fig"

# Load ----

psi <- read.delim("data/240813_psi/240813_psiPerEvent.psi") |>
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


# # filter by neuron
# gene_expression <- read.delim("../majiq/data/2024-03-05_alec_integration/bsn12_subtracted_integrated_binarized_expression_withVDDD_FDR0.05_030424.tsv")
# 
# neurs_sequenced <- colnames(gene_expression)
# 
# 
# filt_psi_lg <- psi_lg |>
#   filter(neuron_id %in% neurs_sequenced) |>
#   left_join(gene_expression |>
#               as.data.frame() |>
#               rownames_to_column("gene_id") |>
#               pivot_longer(-gene_id,
#                            names_to = "neuron_id",
#                            values_to = "expressed")) |>
#   filter(expressed == 1L)


# do not filter by neuron (rely on prefiltering)
neurs_sequenced <- unique(psi_lg$neuron_id)


filt_psi_lg <- psi_lg |>
  filter(! is.na(PSI))


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





max_deltapsi <- function(psi){
  
  psi <- psi[! is.na(psi)]
  n <- length(psi)
  
  if(n < 2) return(NA_real_)
  
  diffs <- lapply(seq_len(n - 1),
                  \(neur){
                    DescTools::CombPairs(psi[neur], psi[-(1:neur)]) |>
                      as.matrix() |>
                      matrixStats::rowDiffs() |>
                      abs() |>
                      as.numeric()
                  })
  dpsis <-  unlist(diffs)
  
  max(dpsis)
}



dpsi_gini <- function(psi){
  
  psi <- psi[! is.na(psi)]
  n <- length(psi)
  
  if(n < 2) return(NA_real_)
  
  diffs <- lapply(seq_len(n - 1),
                  \(neur){
                    DescTools::CombPairs(psi[neur], psi[-(1:neur)]) |>
                      as.matrix() |>
                      matrixStats::rowDiffs() |>
                      abs() |>
                      as.numeric()
                  })
  dpsis <-  unlist(diffs)
  
  DescTools::Gini(dpsis)
}

psi_var <- psi_by_neuron |>
  filter(!is.na(PSI_neuron),
         nb_samples > 2) |>
  summarize(max_dpsi = max_deltapsi(PSI_neuron),
            gini = dpsi_gini(PSI_neuron),
            sd_psi = sd(PSI_neuron, na.rm = TRUE),
            nb_neurons = n(),
            .by = c("event_id", "gene_id","event_type"))





## Plot Gini ----

psi_var_plot <- psi_var |>
  filter(nb_neurons > 5) |>
  filter(event_type != "RI") |>
  mutate(event_type = case_match(event_type,
                                 "A3" ~ "Alt. 3' ss",
                                 "A5" ~ "Alt. 5' ss",
                                 "AF" ~ "Alt. first exon",
                                 "AL" ~ "Alt. last exon",
                                 "MX" ~ "Multiple exons",
                                 "RI" ~ "Intron retention",
                                 "SE" ~ "Cassette exon")) |>
  mutate(gene_name = i2s(gene_id, gids))

psi_var_plot |>
  ggplot() +
  theme_classic() +
  ylab(expression(Gini~index*group("(",group("|",Delta*PSI,"|"),")"))) +
  xlab(expression(max*group("(",group("|",Delta*PSI,"|"),")"))) +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  theme(legend.position = "none") +
  annotate("segment", x = .5, xend = 1, y = .5, linetype = 'dashed', color = 'grey') +
  geom_vline(aes(xintercept = .5), linetype = 'dashed', color = 'grey') +
  facet_wrap(~event_type) +
  geom_point(aes(x = max_dpsi, y = gini, color = event_type),
             alpha = .3) +
  geom_text(data = psi_var_plot |>
              summarize(low_dpsi = sum( max_dpsi <= .5),
                        low_gini = sum(max_dpsi > .5 & gini <= .5),
                        high_gini = sum(max_dpsi > .5 & gini > .5),
                        .by = event_type) |>
              pivot_longer(-event_type,
                           names_to = "class",
                           values_to = "number") |>
              mutate(prop = round(100*number / sum(number), 1),
                     .by = event_type) |>
              mutate(label = paste0(number, "\n(",prop,"%)")) |>
              left_join(tibble(x = c(0, 1, 1),
                               y = c(.1, .1, .9),
                               hjust = c(0,1,1),
                               class = c("low_dpsi","low_gini","high_gini"))),
            aes(x=x,y=y,label=label,hjust = hjust))
  # ggrepel::geom_text_repel(aes(x = max_dpsi, y = gini, label = gene_name))

# ggsave("deltapsi_gini_max.pdf", path = export_dir,
#        width = 10, height = 8, units = "cm", scale = 1.8)

# psi_var |>
#   mutate(gene_name = i2s(gene_id, gids)) |>
#   select(event_id, event_type, gene_id, gene_name, nb_neurons, max_dpsi, gini, sd_psi) |>
#   readr::write_csv(paste0(export_dir, "/psi_gini.csv"))



#~ Gini microexons ----




source("R/extract_event_coordinates.R")




exons_psi_neur <- psi_by_neuron |>
  filter(event_type == "SE") |>
  separate_wider_regex(event_id,
                       patterns = c(gene_id2 = "^WBGene[0-9]{8}", ";",
                                    event_type2 = "[SEA53MXRIFL]{2}", "\\:",
                                    event_coordinates = "[IXV]+\\:[0-9:\\-]+:[+-]$"),
                       cols_remove = FALSE) |>
  select(-event_type2, -gene_id2) |>
  nest(psi = c(neuron_id, PSI_neuron, nb_samples)) |>
  mutate(exon_length = extract_coords("SE", event_coordinates)[["exon_length"]]) |>
  unnest(psi) |>
  select(event_id, gene_id, neuron_id, nb_samples, exon_length, nb_samples, PSI_neuron)

all.equal(psi_by_neuron |> filter(event_type == "SE") |> select(event_id, gene_id,neuron_id,nb_samples,PSI_neuron),
          exons_psi_neur |> select(event_id, gene_id,neuron_id,nb_samples,PSI_neuron))



exons_psi_var <- exons_psi_neur |>
  filter(!is.na(PSI_neuron),
         nb_samples > 2) |>
  mutate(is_microexon = exon_length <= 27) |>
  summarize(max_dpsi = max_deltapsi(PSI_neuron),
            gini = dpsi_gini(PSI_neuron),
            sd_psi = sd(PSI_neuron, na.rm = TRUE),
            nb_neurons = n(),
            .by = c("event_id", "gene_id","is_microexon"))


exons_psi_var_plot <- exons_psi_var |>
  filter(nb_neurons > 5) |>
  mutate(exon_type = if_else(is_microexon,
                              "Microexon",
                              "Longer exon")) |>
  mutate(gene_name = i2s(gene_id, gids))


exons_psi_var_plot |>
  ggplot() +
  theme_classic() +
  ylab(expression(Gini~index*group("(",group("|",Delta*PSI,"|"),")"))) +
  xlab(expression(max*group("(",group("|",Delta*PSI,"|"),")"))) +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_color_manual(values = c("#6F94CD", "#6F94CD")) +
  theme(legend.position = "none") +
  annotate("segment", x = .5, xend = 1, y = .5, linetype = 'dashed', color = 'grey') +
  geom_vline(aes(xintercept = .5), linetype = 'dashed', color = 'grey') +
  facet_wrap(~exon_type) +
  geom_point(aes(x = max_dpsi, y = gini, color = exon_type),
             alpha = .3) +
  geom_text(data = exons_psi_var_plot |>
              summarize(low_dpsi = sum( max_dpsi <= .5),
                        low_gini = sum(max_dpsi > .5 & gini <= .5),
                        high_gini = sum(max_dpsi > .5 & gini > .5),
                        .by = exon_type) |>
              pivot_longer(-exon_type,
                           names_to = "class",
                           values_to = "number") |>
              mutate(prop = round(100*number / sum(number), 1),
                     .by = exon_type) |>
              mutate(label = paste0(number, "\n(",prop,"%)")) |>
              left_join(tibble(x = c(0, 1, 1),
                               y = c(.1, .1, .9),
                               hjust = c(0,1,1),
                               class = c("low_dpsi","low_gini","high_gini"))),
            aes(x=x,y=y,label=label,hjust = hjust))


# ggsave("deltapsi_gini_max_micro.pdf", path = export_dir,
#        width = 6.67, height = 4, units = "cm", scale = 1.8)









####






#~ explorations of Gini deltapsi etc ----

# previous metric, similar to Gini in a way
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


dpsi_mid <- function(psi){
  
  psi <- psi[! is.na(psi)]
  n <- length(psi)
  
  if(n < 2) return(NA_real_)
  
  diffs <- lapply(seq_len(n - 1),
                  \(neur){
                    DescTools::CombPairs(psi[neur], psi[-(1:neur)]) |>
                      as.matrix() |>
                      matrixStats::rowDiffs() |>
                      abs() |>
                      as.numeric()
                  })
  dpsis <-  unlist(diffs)
  
  mean(dpsis > .5)
}


fake_psi <- bind_rows(
  tibble(event_id = rep("a", 7),
         PSI_neuron = rep(.5, 7) + rnorm(7, sd = .01)),
  tibble(event_id = rep("b", 7),
         PSI_neuron = rep(.9, 7) + rnorm(7, sd = .01)),
  tibble(event_id = rep("c", 7),
         PSI_neuron = rep(.1, 7) + rnorm(7, sd = .01)),
  tibble(event_id = rep("d", 7),
         PSI_neuron = rep(.1, 6) |> c(.9) + rnorm(7, sd = .01)),
  tibble(event_id = rep("e", 7),
         PSI_neuron = rep(.1, 4) |> c(rep(.9, 3)) + rnorm(7, sd = .01)),
  tibble(event_id = rep("f", 7),
         PSI_neuron = seq(.1,.9, length.out = 7) + rnorm(7, sd = .01))
)

fake_psi_var <- fake_psi |>
  summarize(specificity = dpsi_mid(PSI_neuron),
            max_dpsi = max_deltapsi(PSI_neuron),
            gini = dpsi_gini(PSI_neuron),
            sd = sd(PSI_neuron, na.rm = TRUE),
            nb_neurons = n(),
            .by = c("event_id"))

fake_psi_var |>
  ggplot() +
  theme_classic() +
  geom_label(aes(x = max_dpsi, y = gini, label = event_id),
             alpha = .1)
fake_psi_var |>
  ggplot() +
  theme_classic() +
  geom_label(aes(x = max_dpsi, y = specificity, label = event_id),
             alpha = .1)

psi_var |>
  # bind_rows(fake_psi_var |> mutate(event_type = "fake")) |>
  filter(nb_neurons > 5) |>
  ggplot() +
  theme_classic() +
  facet_wrap(~event_type) +
  geom_point(aes(x = max_dpsi, y = gini, color = event_type),
             alpha = .5) +
  # ylab(expression(symbol("%")~group("|",Delta*PSI,"|")>0.5)) +
  ylab("Gini index") +
  xlab(expression(max~group("|",Delta*PSI,"|"))) +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  geom_hline(aes(yintercept = .5), linetype = 'dashed', color = 'grey') +
  geom_vline(aes(xintercept = .5), linetype = 'dashed', color = 'grey')


psi_var |>
  bind_rows(fake_psi_var |> mutate(event_type = "fake")) |>
  filter(nb_neurons > 5,
         max_dpsi > .5) |>
  ggplot() +
  theme_classic() +
  facet_wrap(~event_type) +
  geom_point(aes(x = sd, y = gini, color = event_type),
             alpha = .5)
















fake_psi_var |>
  filter(nb_neurons > 5,
         max_dpsi > .5) |>
  ggplot() +
  theme_classic() +
  geom_label(aes(x = sd, y = gini, label = event_id),
             alpha = .1)



psi_var |>
  filter(nb_neurons > 5) |>
  ggplot() +
  theme_classic() +
  facet_wrap(~event_type) +
  geom_histogram(aes(x = max_dpsi, color = event_type),
             alpha = .5) +
  xlab(expression(max~group("|",Delta~PSI,"|")))


# check most variable
ev <- psi_var |>
  filter(specificity > .5, nb_neurons > 5) |>
  slice_head(n = 1) |>
  pull(event_id)

psi <- psi_by_neuron |>
  filter(event_id == ev) |>
  pull(PSI_neuron)
  
plot(psi)

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



psi_var |>
  filter(nb_neurons > 10) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x = sd, y = specificity),
                               alpha = .5) +
  ylab("Specificity index") +
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








# Microexon specificity ----

source("R/extract_event_coordinates.R")

exons_psi_neur <- psi_by_neuron |>
  filter(event_type == "SE") |>
  separate_wider_regex(event_id,
                       patterns = c(gene_id2 = "^WBGene[0-9]{8}", ";",
                                    event_type2 = "[SEA53MXRIFL]{2}", "\\:",
                                    event_coordinates = "[IXV]+\\:[0-9:\\-]+:[+-]$"),
                       cols_remove = FALSE) |>
  select(-event_type2, -gene_id2) |>
  nest(psi = c(neuron_id, PSI_neuron, nb_samples)) |>
  mutate(exon_length = extract_coords("SE", event_coordinates)[["exon_length"]]) |>
  unnest(psi) |>
  select(event_id, gene_id, neuron_id, nb_samples, exon_length, nb_samples, PSI_neuron)




exons_specificity <- exons_psi_neur |>
  filter(!is.na(PSI_neuron),
         nb_samples > 2) |>
  summarize(specificity = min_pairwise_diff(PSI_neuron) / sd(PSI_neuron, na.rm = TRUE),
            sd = sd(PSI_neuron, na.rm = TRUE),
            var = var(PSI_neuron, na.rm = TRUE),
            nb_neurons = n(),
            .by = c("event_id", "gene_id","exon_length"))




exons_specificity |>
  filter(nb_neurons > 2, sd > .1) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x = exon_length, y = specificity),
                               alpha = .5) +
  scale_x_log10() +
  ylab("Specificity index") +
  geom_hline(aes(yintercept = 1), linetype = 'dashed', color = 'grey') +
  geom_vline(aes(xintercept = 27), linetype = 'dashed', color = 'grey')


exons_specificity |>
  filter(nb_neurons > 2, sd > .1) |>
  mutate(microexon = exon_length <= 27) |>
  summarize(specificity = mean(specificity),
            .by = microexon)


exons_specificity |>
  filter(nb_neurons > 2, sd > .1) |>
  mutate(microexon = exon_length <= 27) |>
  mutate(microexon = sample(microexon)) |>
  summarize(specificity = mean(specificity),
            .by = microexon)

exons_specificity |>
  filter(nb_neurons > 2) |>
  mutate(microexon = exon_length <= 27) |>
  # mutate(microexon = sample(microexon)) |>
  summarize(specificity = mean(sd),
            .by = microexon)



exons_specificity |>
  filter(nb_neurons > 2) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x = sd, y = specificity, color = exon_length <= 27),
                               alpha = .5) +
  # scale_x_log10() +
  ylab("Specificity index") +
  geom_hline(aes(yintercept = 1), linetype = 'dashed', color = 'grey')

exons_specificity |>
  filter(nb_neurons > 2) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x = exon_length, y = var),
                               alpha = .5) +
  scale_x_log10() +
  ylab("Variance") +
  geom_vline(aes(xintercept = 27), linetype = 'dashed', color = 'grey')


exons_specificity |>
  filter(nb_neurons > 2) |>
  mutate(microexon = exon_length <= 27) |>
  # mutate(microexon = sample(microexon)) |>
  summarize(specificity = median(var),
            .by = microexon)










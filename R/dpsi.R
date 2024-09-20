


# Inits ----

library(Biostrings)
library(GenomicFeatures)
library(genomation)

library(tidyverse)
library(wbData)


gseq <- readDNAStringSet(wb_get_genome_path(289))


tx2g <- wb_load_tx2gene(289)
gids <- wb_load_gene_ids(289)

export_dir <- "data/outs/240920_fig"






# Load ----

#~ Neurons ----

dpsi_dta <- read.table("data/240918/neurons/240918.dpsi") |>
  rownames_to_column("event_id") |>
  filter(event_id != "WBGene00010673;MX:IV:12574301-12574620:12576543-12576754:12573993-12576902:12577002-12577062:+") |>
  as_tibble()


dpsi <- dpsi_dta |>
  pivot_longer(-event_id,
               names_sep = "_",
               names_to = c("neuron_pair", ".value")) |>
  separate_wider_regex(event_id,
                       patterns = c(gene_id = "^WBGene[0-9]{8}", ";",
                                    event_type = "[SEA53MXRIFL]{2}", "\\:",
                                    event_coordinates = "[IXV]+\\:[0-9:\\-]+:[+-]$"),
                       cols_remove = FALSE) |>
  mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE), .after = "gene_id") |>
  separate_wider_delim(neuron_pair,
                       delim = ".",
                       names = c("neurA", "neurB"),
                       cols_remove = FALSE) |>
  filter(neurA != "Ref",
         neurB != "Ref") |>
  mutate(detectable = !is.na(dPSI)) |>
  mutate(is_ds = detectable & p.val < .1 & abs(dPSI) > .3)




# qs::qsave(dpsi, "intermediates/240918/240920_dpsi_neurons.qs")




#~ Tissue ----

dpsi_dta_tissue <- read.table("data/240918/koterniak/240918.dpsi") |>
  rownames_to_column("event_id") |>
  filter(event_id != "WBGene00010673;MX:IV:12574301-12574620:12576543-12576754:12573993-12576902:12577002-12577062:+") |>
  as_tibble()


dpsi_tissue <- dpsi_dta_tissue |>
  pivot_longer(-event_id,
               names_sep = "_",
               names_to = c("neuron_pair", ".value")) |>
  separate_wider_regex(event_id,
                       patterns = c(gene_id = "^WBGene[0-9]{8}", ";",
                                    event_type = "[SEA53MXRIFL]{2}", "\\:",
                                    event_coordinates = "[IXV]+\\:[0-9:\\-]+:[+-]$"),
                       cols_remove = FALSE) |>
  mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE), .after = "gene_id") |>
  separate_wider_delim(neuron_pair,
                       delim = ".",
                       names = c("neurA", "neurB"),
                       cols_remove = FALSE) |>
  mutate(detectable = !is.na(dPSI)) |>
  mutate(is_ds = detectable & p.val < .1 & abs(dPSI) > .3)

# qs::qsave(dpsi_tissue, "intermediates/240918/240920_dpsi_tissue.qs")





#~ check a few random examples ----



my_ev <- sample(dpsi$event_id[which(abs(dpsi$dPSI) > .3)], 1)
my_neurA <- sample(dpsi$neurA, 1)
my_neurB <- sample(dpsi$neurB, 1)

dpsi |>
  filter(event_id == my_ev,
         neurA %in% c(my_neurA, my_neurB),
         neurB %in% c(my_neurA, my_neurB)) |>
  as.data.frame()









#/ ======= Load ====== / ----

# First perform preprocessing above to create the qs files.
# Upon subsequent reloading, start here directly after Inits


dpsi <- qs::qread("intermediates/240918/240920_dpsi_neurons.qs")
dpsi_tissue <- qs::qread("intermediates/240918/240920_dpsi_tissue.qs")








# Number DS ----
# note: already corrected for multiple testing (option -gc)

dpsi$p.val |> hist(breaks = 70)
dpsi_tissue$p.val |> hist(breaks = 70)

table(dpsi$p.val < .05)
table(dpsi$p.val < .1)
table(dpsi_tissue$p.val < .05)
table(dpsi_tissue$p.val < .1)

dpsi |>
  filter(detectable) |>
  slice_sample(n = 1e4) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = dPSI,
                 y = -log(p.val),
                 color = is_ds),
             alpha = .5) +
  # scale_y_continuous(limits = c(1,NA)) +
  facet_grid(cols = vars(event_type)) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))




# Event types overview ----

dpsi |>
  select(event_id, event_type) |>
  distinct() |>
  pull(event_type) |>
  table()

dpsi |>
  filter(detectable) |>
  select(event_id, event_type) |>
  distinct() |>
  pull(event_type) |>
  table()




dpsi_tissue |>
  filter(detectable) |>
  select(event_id, event_type) |>
  distinct() |>
  pull(event_type) |>
  table()






#~ type vs DS ----

# minimal version
dpsi |>
  summarize(has_ds = any(is_ds),
            .by = c(event_id, event_type)) |>
  select(event_type, has_ds) |>
  ggplot() +
  theme_classic() +
  geom_bar(aes(x = event_type, fill = has_ds)) +
  scale_fill_manual(values = c("grey", "darkred")) +
  xlab(NULL) + ylab("Number of events") +
  theme(legend.position = "none")




# type vs DS vs detectable (code for fig in preprint)

dpsi |>
  summarize(gene_expressed = any(detectable),
            has_ds = any(is_ds),
            .by = c(event_id, event_type)) |>
  mutate(
    category = case_when(
      !gene_expressed ~ "not measured",
      !has_ds ~ "no dAS in neurons",
      .default = "dAS",
    ) |> fct_inorder(),
    event_type = case_match(
      event_type,
      "A3" ~ "Alt. 3' ss",
      "A5" ~ "Alt. 5' ss",
      "AF" ~ "Alt. first exon",
      "AL" ~ "Alt. last exon",
      "MX" ~ "Multiple exons",
      "RI" ~ "Intron retention",
      "SE" ~ "Cassette exon"
    )) |>
  select(event_type, category) |>
  ggplot() +
  theme_classic() +
  scale_fill_manual(values = c("grey90", "grey30", "darkred")) +
  xlab(NULL) + ylab("Number of events") +
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  scale_y_continuous(labels = scales::label_comma()) +
  geom_bar(aes(x = event_type, fill = category),
           color = 'black')



# type vs DS vs detectable with tissue-regulated cases


events_categorized <- full_join(
  dpsi |>
    summarize(measurable_neur = any(detectable),
              has_ds_neur = any(is_ds),
              .by = c(event_id, event_type)),
  dpsi_tissue |>
    summarize(measurable_tissue = any(detectable),
              has_ds_tissue = any(is_ds),
              .by = c(event_id, event_type)),
  by = c("event_id", "event_type")
) |>
  mutate(
    category = case_when(
      !measurable_neur & !measurable_tissue ~ "not measured",
      !has_ds_neur & !has_ds_tissue ~ "no DAS",
      !has_ds_neur &  has_ds_tissue ~ "tissue-regulated",
      has_ds_neur ~ "neuron-regulated",
      .default = "other",
    ) |> ordered(levels = c("not measured", "no DAS",
                            "tissue-regulated", "neuron-regulated")),
    event_type = case_match(
      event_type,
      "A3" ~ "Alt. 3' ss",
      "A5" ~ "Alt. 5' ss",
      "AF" ~ "Alt. first exon",
      "AL" ~ "Alt. last exon",
      "MX" ~ "Multiple exons",
      "RI" ~ "Intron retention",
      "SE" ~ "Cassette exon"
    ))


events_categorized |>
  select(event_type, category) |>
  ggplot() +
  theme_classic() +
  scale_fill_manual(values = c("grey90", "grey30", 'darkblue', "darkred")) +
  xlab(NULL) + ylab("Number of events") +
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  scale_y_continuous(labels = scales::label_comma()) +
  geom_bar(aes(x = event_type, fill = category),
           color = 'black')


# Save plot and corresponding Source Data

# ggsave("ds_per_type.pdf", path = export_dir,
#        width = 16, height = 9, units = "cm")
# 
# events_categorized |>
#   writexl::write_xlsx(file.path(export_dir, "ds_per_type_sourceData.xlsx"))




#~ Manual exploration compare neurons and tissue sero vs pan DAS ----

# serotonergic: keep only comparisons of a sero neuron and a non-sero, check if DAS
neuron_ds_sero_vs_all <- dpsi |>
  mutate(sero = neurA %in% c("ADF" ,"AFD","HSN","NSM") + neurB %in% c("ADF" ,"AFD","HSN","NSM")) |>
  filter(sero == 1) |>
  summarize(neuron_ds_sero = any(is_ds),
            .by = c(event_id, event_type))

tissue_ds_sero_vs_all <- dpsi_tissue |>
  mutate(sero = (neurA =="serotonergic" & neurB == "neurons") |
           (neurA =="neurons" & neurB == "serotonergic")) |>
  filter(sero == 1) |>
  summarize(tissue_ds_sero = any(is_ds),
            .by = c(event_id, event_type))

event_id_both <- union(neuron_ds_sero_vs_all$event_id,
                       tissue_ds_sero_vs_all$event_id)

table(CeNGEN = (neuron_ds_sero_vs_all$neuron_ds_sero |>
                  setNames(neuron_ds_sero_vs_all$event_id))[event_id_both],
      Koterniak = (tissue_ds_sero_vs_all$tissue_ds_sero |>
                     setNames(tissue_ds_sero_vs_all$event_id))[event_id_both])

# example fail
left_join(neuron_ds_sero_vs_all,
          tissue_ds_sero_vs_all,
          by = c("event_id", "event_type")) |>
  filter(neuron_ds_sero,  !tissue_ds_sero) |> #pull(event_type) |> table()
  filter(event_type == "SE") |>
  slice_sample(n = 3) |>
  pull(event_id) -> myev
# myev <- "WBGene00015512;SE:II:7776296-7776351:7776505-7777667:-"

# both agree that there is DAS of sero vs all
myev <- "WBGene00006438;SE:IV:9119961-9120195:9120266-9127195:+"
myev <- "WBGene00016415;SE:II:5204991-5205444:5205632-5206578:+"
myev <- "WBGene00017053;SE:IV:7227507-7227586:7227617-7227696:-"


# Koterniak DAS not CeNGEN
myev <- "WBGene00001851;SE:V:20783222-20783326:20783528-20783738:+"
myev <- "WBGene00015512;SE:II:7776296-7776351:7776505-7777667:-"
myev <- "WBGene00019466;SE:V:9343416-9343501:9343567-9343676:-"

# CeNGEN DAS not Koterniak
myev <- "WBGene00009304;SE:I:14832132-14832841:14832867-14834698:-"
myev <- "WBGene00000201;SE:III:12052714-12053486:12053506-12053586:-"
myev <- "WBGene00018335;SE:IV:8609921-8610384:8610449-8611376:+"


dpsi |>
  filter(event_id == myev,
         neurA %in% c("ADF" ,"AFD","HSN","NSM") | neurB %in% c("ADF" ,"AFD","HSN","NSM"),
         is_ds == TRUE) |> View()

dpsi_tissue |>
  filter(event_id == myev,
         is_ds == TRUE) |> View()

# psi_tissue |>
#   filter(event_id == myev) |>
#   as.data.frame()

psi_lg |>
  filter(event_id == myev,
         neuron_id %in% (unique(psi_lg$neuron_id) |> setdiff(c("ADF" ,"AFD","HSN","NSM")) |> sample(4)))
psi_lg |>
  filter(event_id == myev,
         neuron_id %in% c("ADF" ,"AFD","HSN","NSM", "OLL"))


psi_lg_tissue |>
  filter(event_id == myev,
         neuron_id %in% c("neurons","serotonergic")) |> as.data.frame()



# >> Sequence features << ----

# below this point, only compare whether difference between DAS and non-DAS


#~ Extract coordinates of event ----
source("R/extract_event_coordinates.R")

coords_all <- dpsi |>
  select(event_id, event_type, gene_id, gene_name, event_coordinates) |>
  distinct() |>
  nest(.by = event_type) |>
  mutate(coords = map2(event_type, data,
                       ~ extract_coords(.x, .y[["event_coordinates"]]))) |>
  mutate(coords = map2(data, coords, ~bind_cols(.x["event_id"],
                                                .y)))




#~ GC content ----

source("R/sequence_properties.R")

seq_gc <- map2(coords_all$event_type |> set_names(),
               coords_all$coords,
               \(ev_type, coords){
                 get_seq_properties(ev_type,
                                    coords |> select(-matches("^c[1-9]$")))
                 
               }) |>
  imap(~ {
    .x |>
      add_column(type = .y) |>
      pivot_longer(-c(event_id, type),
                   names_pattern = "^([a-z_]+)_(width|gc)$",
                   names_to = c("measurement", ".value")) |>
      mutate(percent_GC = 100 * gc/width)
  })





#~ Conservation ----


source("R/sequence_conservation.R")
bw_phastcons <- rtracklayer::import("data/UCSC_fastcons/ce11.phastCons26way.bw")

seqlevels(bw_phastcons) <- c("I", "II","III","IV", "M", "V","X")



cons_all_nested <- coords_all |>
  mutate(conservation_score = map2(event_type, coords,
                                   ~ get_cons(.x, .y, bw_phastcons) ))
cons_all <- cons_all_nested |>
  pluck("conservation_score") |>
  set_names(cons_all_nested$event_type)



#~ Assemble data ----

# Assemble the data from these different sequence features
# Note the event types have different column names; first transform the sequence features
# into long dataframes (nested), then join by event id and feature

coords_long <- coords_all |>
  select(-data) |>
  mutate(
    coords = map(coords,
                 \(coords_one_type){
                   coords_one_type |>
                     select(event_id, ends_with("_length")) |>
                     pivot_longer(cols = ends_with("_length"),
                                  names_pattern = "^([_a-z]+)_length$",
                                  names_to = "feature",
                                  values_to = "length")
                 })
  ) |>
  unnest(coords)

cons_long <- tibble(
  event_type = names(cons_all),
  values = map(cons_all, 
             \(cons_one_type){
               cons_one_type |>
                 pivot_longer(cols = ends_with("_conservation"),
                              names_pattern = "^([a-z_]+)_conservation$",
                              names_to = "feature",
                              values_to = "conservation")
             } )
) |>
  unnest(values)

gc_long <- tibble(
  event_type = names(seq_gc),
  values = map(seq_gc,
                   \(gc_one_type){
                     gc_one_type |>
                       rename(feature = measurement) |>
                       select(event_id, feature, percent_GC)
                   })
) |>
  unnest(values)

features_long <- coords_long |>
  left_join(cons_long,
            by = c("event_type", "event_id", "feature")) |>
  left_join(gc_long,
            by = c("event_type", "event_id", "feature"))


local({
  # ensure perfect match ie left or right join identical
  stopifnot(all.equal(
    features_long,
    coords_long |>
      right_join(cons_long,
                 by = c("event_type", "event_id", "feature")) |>
      right_join(gc_long,
                 by = c("event_type", "event_id", "feature"))
  ))
})




# add back DAS information


local({
  # ensure match of event names
  stopifnot(all.equal(
    features_long |>
      left_join(dpsi |>
                  summarize(has_ds = any(is_ds),
                            is_detectable = any(detectable),
                            .by = c(event_type, event_id)),
                by = c("event_type", "event_id")
      ),
    features_long |>
      full_join(dpsi |>
                  summarize(has_ds = any(is_ds),
                            is_detectable = any(detectable),
                            .by = c(event_type, event_id)),
                by = c("event_type", "event_id")
      )
  ))
})

features_long <- features_long |>
  full_join(dpsi |>
              summarize(has_ds = any(is_ds),
                        is_detectable = any(detectable),
                        .by = c(event_type, event_id)),
            by = c("event_type", "event_id")
            )




# qs::qsave(features_long,
#           "intermediates/240918/240920_features.qs")




#~ Plots ----


features_long |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x = has_ds, y = length, color = has_ds)) +
  facet_wrap(~interaction(event_type, feature), scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Length (bp)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))


features_long |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x = has_ds, y = conservation, color = has_ds)) +
  facet_wrap(~interaction(event_type, feature), scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Conservation score") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))

features_long |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x = has_ds, y = percent_GC, color = has_ds)) +
  facet_wrap(~interaction(event_type, feature), scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Percent GC") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))










# Stats ----
tests <- features_long |>
  filter(is_detectable) |>
  group_by(event_type, feature) |>
  summarize(pval_length = wilcox.test(length[has_ds],
                                      length[! has_ds])$p.value,
            pval_gc = wilcox.test(percent_GC[has_ds],
                                  percent_GC[! has_ds])$p.value,
            pval_conservation = wilcox.test(conservation[has_ds],
                                            conservation[! has_ds])$p.value,
            .groups = 'drop') |>
  pivot_longer(starts_with("pval_"),
               names_prefix = "pval_",
               names_to = "metric",
               values_to = "p_val") |>
  mutate(padj = p.adjust(p_val, method = 'fdr'))



medians <- features_long |>
  filter(is_detectable) |>
  group_by(event_type, feature) |>
  summarize(median_length_ds = median(length[has_ds]),
             median_length_nonds = median(length[! has_ds]),
             median_gc_ds = median(percent_GC[has_ds]),
             median_gc_nonds = median(percent_GC[! has_ds]),
             median_conservation_ds = median(conservation[has_ds]),
             median_conservation_nonds = median(conservation[! has_ds]),
             .groups = 'drop') |>
  pivot_longer(starts_with("median_"),
               names_pattern = "^(median_)([a-z]+)_(nonds|ds)$",
               names_to = c(".value", "metric", ".value"))





table(tests$padj < .1)

hist(tests$p_val, breaks = 40)
hist(tests$padj, breaks = 40)

tests |>
  filter(padj < 0.05) |>
  mutate(sig = cut(padj,
                   breaks =   c(.1,  .05,  .01, .001,  0),
                   labels = rev(c("#", "*", "**", "***"))))




tests |>
  left_join(medians,
            by = c("event_type", "feature", "metric")) |>
  filter(padj < 0.1) |>
  mutate(sig = cut(padj,
                   breaks =c(0, .001, .01, .05, .1),
                   labels =   c("***", "**", "*", "#"))) |>
  select(event_type, feature, metric, median_ds, median_nonds, padj, sig) |>
  as.data.frame()

# see outs/sig_table.xlsx

# tests |>
#   left_join(medians,
#             by = c("event_type", "feature", "metric")) |>
#   # filter(padj < 0.1) |>
#   mutate(sig = cut(padj,
#                    breaks =c(0, .001, .01, .05, .1),
#                    labels =   c("***", "**", "*", "#"))) |>
#   rename(`p-value` = p_val,
#          `adjusted p-value` = padj,
#          `median in dAS events` = median_ds,
#          `median in non-dAS events` = median_nonds,
#          `significance level` = sig) |>
#   writexl::write_xlsx(paste0(export_dir, "/all_comparisons.xlsx"))



#~ Focus on some of the significant changes ----


# se exon length

d_se_med <- d_se |>
  filter(is_detectable,
         feature == "exon") |>
  summarize(median = median(length),
            .by = has_ds)

d_se |>
  filter(is_detectable) |>
  filter(feature == "exon") |>
  ggplot() +
  theme_classic() +
  geom_density(aes(x = length, fill = has_ds),
               alpha = .5) +
  geom_vline(aes(xintercept = median, color = has_ds),
             data = d_se_med,
             linetype = 'dotted', linewidth = 1) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_color_manual(values = c("grey30", "darkred")) +
  scale_y_continuous(limits = c(0,1.2)) +
  theme(legend.position = "none") +
  xlab("Skipped exon length (bp)")

d_se |>
  filter(is_detectable) |>
  filter(feature == "exon") |>
  mutate(has_ds = if_else(has_ds, "dAS", "non-dAS") |>
           factor(levels = c("non-dAS", "dAS"))) |>
  ggplot() +
  theme_classic() +
  ggridges::geom_density_ridges(aes(x = length, y = has_ds, fill = has_ds),
                                scale = 3) +
  geom_vline(aes(xintercept = median, color = has_ds),
             data = d_se_med,
             linetype = 'dashed', linewidth = 1) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_color_manual(values = c("grey50", "red3")) +
  scale_y_discrete(expand = c(0.01, 0), labels = NULL) +
  # scale_y_continuous(limits = c(0,1.2)) +
  theme(legend.position = "none") +
  xlab("Skipped exon length (bp)") + ylab(NULL)


# ggsave("ridge_se_exon_length_density.pdf", path = export_dir,
#        width = 5, height = 5, units = "cm")



# AF dist exon length

d_af_de_med <- d_af |>
  filter(is_detectable,
         feature == "distal_exon") |>
  summarize(median = median(length),
            .by = has_ds)

# d_af |>
#   filter(is_detectable) |>
#   filter(feature == "distal_exon") |>
#   ggplot() +
#   theme_classic() +
#   geom_density(aes(x = length, fill = has_ds),
#                alpha = .5) +
#   scale_x_log10(labels = scales::label_comma()) +
#   scale_fill_manual(values = c("grey30", "darkred")) +
#   scale_y_continuous(limits = c(0,1.2)) +
#   theme(legend.position = "none") +
#   xlab("Alternative first exon, distal exon length (bp)")
# 
# # ggsave("af_distal_exon_length_density.pdf", path = export_dir,
# #        width = 10, height = 7, units = "cm")

d_af |>
  filter(is_detectable) |>
  filter(feature == "distal_exon") |>
  ggplot() +
  theme_classic() +
  ggridges::geom_density_ridges(aes(x = length, y = has_ds, fill = has_ds),
                                scale = 3) +
  geom_vline(aes(xintercept = median, color = has_ds),
             data = d_af_de_med,
             linetype = 'dashed', linewidth = 1) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_color_manual(values = c("grey50", "red3")) +
  scale_y_discrete(expand = c(0.01, 0), labels = NULL) +
  # scale_y_continuous(limits = c(0,1.2)) +
  theme(legend.position = "none") +
  xlab("Alternative first exon, distal exon length (bp)") + ylab(NULL)

# ggsave("ridge_af_distal_exon_length_density.pdf", path = export_dir,
#        width = 5, height = 5, units = "cm")




# af dist intron length


# d_af |>
#   filter(is_detectable) |>
#   filter(feature == "distal_intron") |>
#   ggplot() +
#   theme_classic() +
#   geom_density(aes(x = length, fill = has_ds),
#                alpha = .5) +
#   scale_x_log10(limits = c(100, 40000)) +
#   scale_fill_manual(values = c("grey30", "darkred")) +
#   scale_y_continuous(limits = c(0,1.2)) +
#   theme(legend.position = "none") +
#   xlab("Alternative first exon, distal intron length (bp)")
# 
# # ggsave("af_distal_intron_length_density.pdf", path = export_dir,
# #        width = 15, height = 8, units = "cm")

d_af_di_med <- d_af |>
  filter(is_detectable,
         feature == "distal_intron") |>
  summarize(median = median(length),
            .by = has_ds)

d_af |>
  filter(is_detectable) |>
  filter(feature == "distal_intron") |>
  ggplot() +
  theme_classic() +
  ggridges::geom_density_ridges(aes(x = length, y = has_ds, fill = has_ds),
                                scale = 3) +
  geom_vline(aes(xintercept = median, color = has_ds),
             data = d_af_di_med,
             linetype = 'dashed', linewidth = 1) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_color_manual(values = c("grey50", "red3")) +
  scale_y_discrete(expand = c(0.01, 0), labels = NULL) +
  # scale_y_continuous(limits = c(0,1.2)) +
  theme(legend.position = "none") +
  xlab("Alternative first exon, distal intron length (bp)") + ylab(NULL)

# ggsave("ridge_af_distal_intron_length_density.pdf", path = export_dir,
#        width = 5, height = 5, units = "cm")



# af prox exon length

d_af_pe_med <- d_af |>
  filter(is_detectable,
         feature == "proximal_exon") |>
  summarize(median = median(length),
            .by = has_ds)


# d_af |>
#   filter(is_detectable) |>
#   filter(feature == "proximal_exon") |>
#   ggplot() +
#   theme_classic() +
#   geom_density(aes(x = length, fill = has_ds),
#                alpha = .5) +
#   scale_x_log10() +
#   scale_fill_manual(values = c("grey30", "darkred")) +
#   scale_y_continuous(limits = c(0,1.2)) +
#   theme(legend.position = "none") +
#   xlab("Alternative first exon, proximal exon length (bp)")
# 
# # ggsave("af_proximal_exon_length_density.pdf", path = export_dir,
# #        width = 15, height = 8, units = "cm")


d_af |>
  filter(is_detectable) |>
  filter(feature == "proximal_exon") |>
  ggplot() +
  theme_classic() +
  ggridges::geom_density_ridges(aes(x = length, y = has_ds, fill = has_ds),
                                scale = 3, bandwidth = .1) +
  geom_vline(aes(xintercept = median, color = has_ds),
             data = d_af_pe_med,
             linetype = 'dashed', linewidth = 1) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_color_manual(values = c("grey50", "red3")) +
  scale_y_discrete(expand = c(0.01, 0), labels = NULL) +
  # scale_y_continuous(limits = c(0,1.2)) +
  theme(legend.position = "none") +
  xlab("Alternative first exon, proximal exon length (bp)") + ylab(NULL)

# ggsave("ridge_af_proximal_exon_length_density.pdf", path = export_dir,
#        width = 5, height = 5, units = "cm")






# AF dist intr conservation
d_af_di_cons_med <- d_af |>
  filter(is_detectable,
         feature == "distal_intron") |>
  summarize(median = median(conservation),
            .by = has_ds)



# d_af |>
#   filter(is_detectable) |>
#   filter(feature == "distal_intron") |>
#   ggplot() +
#   theme_classic() +
#   geom_density(aes(x = conservation, fill = has_ds),
#                alpha = .5, bounds = c(0,1)) +
#   scale_fill_manual(values = c("grey30", "darkred")) +
#   scale_x_continuous(limits = c(0,1)) +
#   theme(legend.position = "none") +
#   xlab("Alternative first exon, distal intron conservation score")
# 
# # ggsave("af_distal_intron_conservation_density.pdf", path = export_dir,
# #        width = 15, height = 8, units = "cm")




d_af |>
  filter(is_detectable) |>
  filter(feature == "distal_intron") |>
  ggplot() +
  theme_classic() +
  ggridges::geom_density_ridges(aes(x = conservation, y = has_ds, fill = has_ds),
                                scale = 3) +
  geom_vline(aes(xintercept = median, color = has_ds),
             data = d_af_di_cons_med,
             linetype = 'dashed', linewidth = 1) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_color_manual(values = c("grey50", "red3")) +
  scale_y_discrete(expand = c(0.01, 0), labels = NULL) +
  # scale_y_continuous(limits = c(0,1.2)) +
  theme(legend.position = "none") +
  xlab("Alternative first exon, distal intron conservation score") + ylab(NULL)

# ggsave("ridge_af_distal_intron_conservation_density.pdf", path = export_dir,
#        width = 5, height = 5, units = "cm")













# Microexons ----

# coords_all <- qs::qread("intermediates/240425_dpsi/240425_all_coords.qs")
# dpsi <- qs::qread("intermediates/240425_dpsi/filt_dpsi.qs")

skipped_exons_lengths <- coords_all |>
  filter(event_type == "SE") |>
  pull(coords) |>
  first() |>
  select(event_id, exon_length)

skipped_exons <- dpsi |>
  filter(event_type == "SE") |>
  left_join(skipped_exons_lengths)

skipped_exons_lengths |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = exon_length), bins = 50,
                 color = 'white') +
  scale_x_log10() +
  geom_vline(aes(xintercept = 27),
             color = 'red3')


skipped_exons2 <- skipped_exons |>
  filter(detectable) |>
  summarize(has_ds = any(is_ds),
            .by = c(exon_length, event_id))

lim_high <- 27
skipped_exons2 |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = exon_length,
                     fill = has_ds),
                 bins = 30,
                 color = 'white') +
  scale_x_log10(labels = scales::label_comma()) +
  geom_vline(aes(xintercept = 27),
             color = 'grey10', linetype = "dashed") +
  geom_vline(aes(xintercept = lim_high),
             color = 'grey10', linetype = "dashed") +
  scale_fill_manual(values = c("grey30", "darkred")) +
  geom_text(data = tibble(
    x = c(11, 3000),
    y = c(90,90),
    label = c(paste0("Microexons (≤ 27 bp)\n",
                     sum(skipped_exons2$has_ds[skipped_exons2$exon_length <= 27]),
                     "/",
                     sum(skipped_exons2$exon_length <= 27),
                     " (", round(100*mean(skipped_exons2$has_ds[skipped_exons2$exon_length <= 27]))," %)"),
              paste0("Longer exons (> ",lim_high," bp)\n",
                     sum(skipped_exons2$has_ds[skipped_exons2$exon_length > lim_high]),
                     "/",
                     sum(skipped_exons2$exon_length > lim_high),
                     " (", round(100*mean(skipped_exons2$has_ds[skipped_exons2$exon_length > lim_high]))," %)"))
  ),
  aes(x=x,y=y,label=label)) +
  theme(legend.position = "none") +
  xlab("Exon length (bp)") +
  ylab("Number of exons")

# ggsave("microexons.pdf", path = export_dir,
#        width = 12, height = 8, units = "cm")

# test
skp_ex_test <- skipped_exons2 |>
  mutate(is_microexon = exon_length <= 27) |>
  summarize(nb_ds = sum(has_ds),
            nb_tot = n(),
            .by = is_microexon)

prop.test(x = skp_ex_test$nb_ds,
          n = skp_ex_test$nb_tot)




#~ Microexons: nb of diff pairs ----

skipped_ex_by_pairs <- skipped_exons |>
  filter(detectable) |>
  # mutate(is_microexon = (exon_length <= 27)) |>
  summarize(nb_ds = sum(is_ds),
            nb_tot = n(),
            nb_non_ds = sum(!is_ds),
            .by = c(event_id, exon_length))


stopifnot(all(skipped_ex_by_pairs$nb_ds + skipped_ex_by_pairs$nb_non_ds == skipped_ex_by_pairs$nb_tot))

# skipped_ex_by_pairs |>
#   ggplot() +
#   theme_classic() +
#   geom_point(aes(x = exon_length, y = nb_ds/nb_tot)) +
#   scale_x_log10(labels = scales::label_comma()) +
#   geom_vline(aes(xintercept = 27),
#              color = 'grey50', linetype = "dashed")
# 
# skipped_ex_by_pairs |>
#   filter(nb_tot > 10) |>
#   mutate(is_microexon = exon_length <= 27,
#          prop_ds = nb_ds/nb_tot) |>
#   summarize(mean_prop_ds = mean(prop_ds),
#             sd_prop_ds = sd(prop_ds),
#             .by = is_microexon) |>
#   ggplot() +
#   theme_classic() +
#   geom_col(aes(x = is_microexon, y = prop_ds))
# 
# skipped_ex_by_pairs |>
#   filter(nb_tot > 10) |>
#   mutate(is_microexon = exon_length <= 27,
#          prop_ds = nb_ds/nb_tot) |>
#   ggplot() +
#   theme_classic() +
#   geom_density(aes(x = prop_ds, fill = is_microexon),
#                alpha = .2)


pos <- ggbeeswarm::position_quasirandom()
skipped_ex_by_pairs |>
  filter(nb_tot > 10) |>
  mutate(is_microexon = exon_length <= 27,
         prop_ds = nb_ds/nb_tot) |>
  mutate(gene_name = if_else((is_microexon & prop_ds > .15) |
                               (!is_microexon & prop_ds > .38),
                             str_split_i(event_id, ";",1) |>
                               i2s(gids, warn_missing = TRUE),
                             "")) |>
  mutate(type = if_else(is_microexon, "microexon", "longer exon")) |>
  ggplot(aes(x = type, y = prop_ds, label = gene_name)) +
  theme_classic() +
  geom_point(position = pos) +
  ggrepel::geom_text_repel(box.padding = 1.2, position = pos, color = 'grey40') +
  scale_y_continuous(labels = scales::label_percent()) +
  xlab(NULL) + ylab("Proportion of neuron pairs dAS")


# ggsave("microexons_ds_pairs.pdf", path = export_dir,
#        width = 10, height = 9, units = "cm")



skipped_ex_by_pairs |>
  filter(nb_tot > 10) |>
  mutate(is_microexon = exon_length <= 27,
         prop_ds = nb_ds/nb_tot) |>
  summarize(mean_prop_ds = mean(prop_ds),
            .by = is_microexon)


skipped_ex_by_pairs_filt <- skipped_ex_by_pairs |>
  filter(nb_tot > 10) |>
  mutate(is_microexon = exon_length <= 27,
         prop_ds = nb_ds/nb_tot)

wilcox.test(prop_ds ~ is_microexon, data = skipped_ex_by_pairs_filt)
boxplot(prop_ds ~ is_microexon, data = skipped_ex_by_pairs_filt)



#~ export microexons ----

# skipped_exons2 |>
#   mutate(is_microexon = exon_length <= 27) |>
#   filter(is_microexon, has_ds) |>
#   select(event_id) |>
#   write_csv("data/outs/240426_microexons.csv")






# Inits ----

library(Biostrings)
library(GenomicFeatures)
library(plyranges)
library(genomation)

library(tidyverse)
library(wbData)


gseq <- readDNAStringSet(wb_get_genome_path(289))


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

# > Fig. 6A 

# ---> source data

# events_categorized |>
#   writexl::write_xlsx(file.path(export_dir, "ds_per_type_sourceData.xlsx"))




#~ Manual exploration compare neurons and tissue sero vs pan DAS ----
# Mentioned in discussion, not shown

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
# note, results saved in features_long qs file, can skip to end of the "Assemble data" subsection

#~ Extract coordinates of event ----
source("R/extract_event_coordinates.R")

coords_all <- dpsi |>
  filter(event_type != "RI") |>
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
                  filter(event_type != "RI") |>
                  summarize(has_ds = any(is_ds),
                            is_detectable = any(detectable),
                            .by = c(event_type, event_id)),
                by = c("event_type", "event_id")
      ),
    features_long |>
      full_join(dpsi |>
                  filter(event_type != "RI") |>
                  summarize(has_ds = any(is_ds),
                            is_detectable = any(detectable),
                            .by = c(event_type, event_id)),
                by = c("event_type", "event_id")
      )
  ))
})

features_long <- features_long |>
  full_join(dpsi |>
              filter(event_type != "RI") |>
              summarize(has_ds = any(is_ds),
                        is_detectable = any(detectable),
                        .by = c(event_type, event_id)),
            by = c("event_type", "event_id")
            )




# qs::qsave(features_long,
#           "intermediates/240918/241011_features.qs")

features_long <- qs::qread("intermediates/240918/241011_features.qs")


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





table(tests$padj < .05)

hist(tests$p_val, breaks = 40)
hist(tests$padj, breaks = 40)

tests |>
  filter(padj < 0.05) |>
  mutate(sig = cut(padj,
                   breaks =   c(.1,  .05,  .01, .001,  0),
                   labels = rev(c("#", "*", "**", "***"))))




table_signif_seq_feats <- tests |>
  left_join(medians,
            by = c("event_type", "feature", "metric")) |>
  filter(padj < 0.1) |>
  mutate(sig = cut(padj,
                   breaks =c(0, .001, .01, .05, .1),
                   labels =   c("***", "**", "*", "#"))) |>
  select(event_type, feature, metric, median_ds, median_nonds, padj, sig)

table_signif_seq_feats |>
  arrange(padj)

# Save results (table 1)

# table_signif_seq_feats |>
#   arrange(padj) |>
#   writexl::write_xlsx(file.path(export_dir, "sequence_features_signif.xlsx"))


# save all seq features (supplementary table)

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
#   writexl::write_xlsx(file.path(export_dir, "sequence_features_all_tests.xlsx"))

# save all events for source data

# features_long |>
#   writexl::write_xlsx(file.path(export_dir, "sequence_feat_all_sourceData.xlsx"))


#~ Focus on some of the significant changes ----

#~~ AF prox exon length ----

d_af_pe_med <- features_long |>
  filter(event_type == "AF",
         feature == "proximal_exon") |>
  filter(is_detectable) |>
  summarize(median = median(length),
            .by = has_ds)


features_long |>
  filter(event_type == "AF",
         feature == "proximal_exon") |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggridges::geom_density_ridges(aes(x = length, y = has_ds, fill = has_ds),
                                scale = 3, bandwidth = .1) +
  geom_vline(aes(xintercept = median, color = has_ds),
             data = d_af_pe_med,
             linetype = 'dashed', linewidth = .4) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_color_manual(values = c("grey50", "red3")) +
  scale_y_discrete(expand = c(0.01, 0), labels = NULL) +
  theme(legend.position = "none") +
  xlab("Alternative first exon,\nproximal exon length (bp)") + ylab(NULL)

# ggsave("ridge_af_pe_length.pdf", path = export_dir,
#        width = 5, height = 4.5, units = "cm")





#~~ SE ex length ----
d_se_med <- features_long |>
  filter(event_type == "SE",
         feature == "exon") |>
  filter(is_detectable) |>
  summarize(median = median(length),
            .by = has_ds)


features_long |>
  filter(event_type == "SE",
         feature == "exon") |>
  filter(is_detectable) |>
  mutate(has_ds = if_else(has_ds, "dAS", "non-dAS") |>
           factor(levels = c("non-dAS", "dAS"))) |>
  ggplot() +
  theme_classic() +
  ggridges::geom_density_ridges(aes(x = length, y = has_ds, fill = has_ds),
                                scale = 3) +
  geom_vline(aes(xintercept = median, color = has_ds),
             data = d_se_med,
             linetype = 'dashed', linewidth = .4) +
  scale_x_log10(labels = scales::label_comma()) +
  # scale_x_continuous(limits = c(0, 500)) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_color_manual(values = c("grey50", "red3")) +
  scale_y_discrete(expand = c(0.01, 0), labels = NULL) +
  theme(legend.position = "none") +
  xlab("Skipped exon,\nexon length (bp)") + ylab(NULL)


# ggsave("ridge_se_exon_length.pdf", path = export_dir,
#        width = 5, height = 4.5, units = "cm")


#~~ AF dist intron length ----

d_af_di_med <- features_long |>
  filter(event_type == "AF",
         feature == "distal_intron") |>
  filter(is_detectable) |>
  summarize(median = median(length),
            .by = has_ds)



features_long |>
  filter(event_type == "AF",
         feature == "distal_intron") |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggridges::geom_density_ridges(aes(x = length, y = has_ds, fill = has_ds),
                                scale = 3) +
  geom_vline(aes(xintercept = median, color = has_ds),
             data = d_af_di_med,
             linetype = 'dashed', linewidth = .4) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_color_manual(values = c("grey50", "red3")) +
  scale_y_discrete(expand = c(0.01, 0), labels = NULL) +
  # scale_y_continuous(limits = c(0,1.2)) +
  theme(legend.position = "none") +
  xlab("Alternative first exon,\ndistal intron length (bp)") + ylab(NULL)

# ggsave("ridge_af_di_length.pdf", path = export_dir,
#        width = 5, height = 4.5, units = "cm")




#~~ A5 intron length ----

d_a5_i_med <- features_long |>
  filter(event_type == "A5",
         feature == "intron") |>
  filter(is_detectable) |>
  summarize(median = median(length),
            .by = has_ds)



features_long |>
  filter(event_type == "A5",
         feature == "intron") |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggridges::geom_density_ridges(aes(x = length, y = has_ds, fill = has_ds),
                                scale = 3) +
  geom_vline(aes(xintercept = median, color = has_ds),
             data = d_a5_i_med,
             linetype = 'dashed', linewidth = .4) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_color_manual(values = c("grey50", "red3")) +
  scale_y_discrete(expand = c(0.01, 0), labels = NULL) +
  # scale_y_continuous(limits = c(0,1.2)) +
  theme(legend.position = "none") +
  xlab("Alternative 5' splice site,\nintron length (bp)") + ylab(NULL)

# ggsave("ridge_a5_i_length.pdf", path = export_dir,
#        width = 5, height = 4.5, units = "cm")







#~~ AF prox intr conservation ----
d_af_pi_cons_med <- features_long |>
  filter(event_type == "AF",
         feature == "proximal_intron") |>
  filter(is_detectable) |>
  summarize(median = median(conservation),
            .by = has_ds)




features_long |>
  filter(event_type == "AF",
         feature == "proximal_intron") |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggridges::geom_density_ridges(aes(x = conservation, y = has_ds, fill = has_ds),
                                scale = 3) +
  geom_vline(aes(xintercept = median, color = has_ds),
             data = d_af_pi_cons_med,
             linetype = 'dashed', linewidth = .4) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_color_manual(values = c("grey50", "red3")) +
  scale_y_discrete(expand = c(0.01, 0), labels = NULL) +
  # scale_y_continuous(limits = c(0,1.2)) +
  theme(legend.position = "none") +
  xlab("Alternative first exon,\nproximal intron conservation score") + ylab(NULL)

# ggsave("ridge_af_pi_cons.pdf", path = export_dir,
#        width = 5, height = 4.5, units = "cm")



#~~ AF dist exon length ----

d_af_de_med <- features_long |>
  filter(event_type == "AF",
         feature == "distal_exon") |>
  filter(is_detectable) |>
  summarize(median = median(length),
            .by = has_ds)



features_long |>
  filter(event_type == "AF",
         feature == "distal_exon") |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggridges::geom_density_ridges(aes(x = length, y = has_ds, fill = has_ds),
                                scale = 3) +
  geom_vline(aes(xintercept = median, color = has_ds),
             data = d_af_de_med,
             linetype = 'dashed', linewidth = .4) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_color_manual(values = c("grey50", "red3")) +
  scale_y_discrete(expand = c(0.01, 0), labels = NULL) +
  # scale_y_continuous(limits = c(0,1.2)) +
  theme(legend.position = "none") +
  xlab("Alternative first exon,\ndistal exon length (bp)") + ylab(NULL)

# ggsave("ridge_af_de_length.pdf", path = export_dir,
#        width = 5, height = 4.5, units = "cm")

# > Fig. 6B


# ---> source data

# features_long |>
#   filter(
#     (event_type == "AF" & feature == "proximal_exon") |
#       (event_type == "SE" & feature == "exon") |
#       (event_type == "AF" & feature == "distal_intron") |
#       (event_type == "AF" & feature == "distal_intron") |
#       (event_type == "A5" & feature == "intron") |
#       (event_type == "AF" & feature == "proximal_intron") |
#       (event_type == "AF" & feature == "distal_exon")
#   ) |>
#   filter(is_detectable) |>
#   writexl::write_xlsx(file.path(export_dir, "sourceDate_feature_densityplots.xlsx"))








# >> Protein sequence << ----

# ~ restrict to CDS-overlapping events ----

source("R/extract_event_coordinates.R")

cds <- GenomicFeatures::cds(wbData::wb_load_TxDb(289))


coords_all <- dpsi |>
  select(event_id, event_type, gene_id, gene_name, event_coordinates) |>
  distinct() |>
  nest(.by = event_type) |>
  mutate(coords = map2(event_type, data,
                       ~ extract_coords(.x, .y[["event_coordinates"]]))) |>
  mutate(coords = map2(data, coords,
                       ~bind_cols(.x["event_id"],
                                  .y)))


# ~~ Alt 3' ss ----

# does the overhang overlap the CDS?

a3_gr <- coords_all |>
  filter(event_type == "A3") |>
  chuck("coords", 1L) |>
  select(-matches("^c[1-4]$"),
         -ends_with("_length")) |>
  pivot_longer(cols = c(intron_start, intron_end, overhang_start, overhang_end),
               names_to = c("feature", ".value"),
               names_sep = "_") |>
  rename(seqnames = chr) |>
  as_granges()


a3_overhang_gr <- a3_gr |>
  filter(feature == "overhang")

a3_overhang_in_cds <- join_overlap_inner_directed(a3_overhang_gr, cds) |>
  as.data.frame() |>
  pull(event_id) |>
  unique()



# ~~ Alt 5' ss ----

# does the overhang overlap the CDS?

a5_gr <- coords_all |>
  filter(event_type == "A5") |>
  chuck("coords", 1L) |>
  select(-matches("^c[1-4]$"),
         -ends_with("_length")) |>
  pivot_longer(cols = c(intron_start, intron_end, overhang_start, overhang_end),
               names_to = c("feature", ".value"),
               names_sep = "_") |>
  rename(seqnames = chr) |>
  as_granges()


a5_overhang_gr <- a5_gr |>
  filter(feature == "overhang")

a5_overhang_in_cds <- join_overlap_inner_directed(a5_overhang_gr, cds) |>
  as.data.frame() |>
  pull(event_id) |>
  unique()



# ~~ SE ----

# does the exon (inclusion) overlap CDS

se_gr <- coords_all |>
  filter(event_type == "SE") |>
  chuck("coords", 1L) |>
  select(event_id, chr, strand,
         intron_start = upstream_intron_start,
         intron_end = downstream_intron_end,
         exon_start, exon_end) |>
  pivot_longer(cols = c(intron_start, intron_end, exon_start, exon_end),
               names_to = c("feature", ".value"),
               names_pattern = "^(intron|exon)_(start|end)$") |>
  rename(seqnames = chr) |>
  as_granges()


se_inclusion <- se_gr |>
  filter(feature == "exon")

se_exons_in_cds <- join_overlap_inner_directed(se_inclusion, cds) |>
  as.data.frame() |>
  pull(event_id) |>
  unique()




length(a3_overhang_in_cds)
length(a3_gr)/2

length(a5_overhang_in_cds)
length(a5_gr)/2

length(se_exons_in_cds)
length(se_gr)/2




# ~ frame-shift ----

proportions_ndf <-  features_long |>
  filter((event_type == "A5" & feature == "overhang" & (event_id %in% a5_overhang_in_cds)) |
           (event_type == "A3" & feature == "overhang" & (event_id %in% a3_overhang_in_cds)) |
           (event_type == "SE" & feature == "exon" & (event_id %in% se_exons_in_cds)) ) |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  nest(.by = event_type) |>
  mutate(tab = map(data,
                   ~ {
                     .x |>
                       summarize(n_frame = sum(length_is_triple),
                                 n_not_frame = sum(!length_is_triple),
                                 .by = c(has_ds, is_detectable)) |>
                       mutate(prop_in_frame = round( 100 * n_frame / (n_frame + n_not_frame) )) |>
                       arrange(is_detectable, has_ds)
                   }
  ))

proportions_l <- proportions_ndf$tab |> set_names(proportions_ndf$event_type)




imap(proportions_l,
     \(tab, ev_t){
       tab |>
         select(-prop_in_frame) |>
         pivot_longer(-c(has_ds, is_detectable),
                      names_to = "frame",
                      values_to = "count") |>
         add_column(event_type = ev_t)
  }) |>
  bind_rows() |>
  mutate(
    category = case_when(
      ! is_detectable ~ "not measured",
      has_ds ~ "DAS",
      .default = "not DAS") |>
      factor(levels = c("not measured", "not DAS", "DAS") |> rev()),
    frame = factor(frame, levels = c("n_not_frame", "n_frame")),
    event_type = case_match(event_type,
                            "A3" ~ "Alt. 3' ss",
                            "A5" ~ "Alt. 5' ss",
                            "SE" ~ "Cassette exon")
  ) |>
  ggplot() +
  theme_classic() +
  guides(fill = guide_legend(title = NULL )) +
  theme(legend.position = "top") +
  scale_fill_brewer(type = "qual",
                    labels = c("Frame-shifting","Frame-preserving")) +
  xlab(NULL) +
  ylab("Number of events") +
  coord_flip() +
  facet_wrap(~ event_type) +
  geom_col(aes(x = category, y = count, fill = frame))


# Save figure

# ggsave("events_length_triple.pdf", path = export_dir,
#        width = 15, height = 6, units = "cm")

# > Fig 6E



# For sourceData Fig 6E

# features_long |>
#   filter((event_type == "A5" & feature == "overhang" ) |
#            (event_type == "A3" & feature == "overhang" ) |
#            (event_type == "SE" & feature == "exon" ) ) |>
#   mutate(length_is_triple = (length %% 3) == 0,
#          intersects_cds = event_id %in% c(a5_overhang_in_cds, a3_overhang_in_cds, se_exons_in_cds)) |>
#   select(event_type, event_id, is_detectable, has_ds, length_is_triple, intersects_cds) |>
#   # count(event_type, is_detectable, has_ds, length_is_triple, intersects_cds) |> View()
#   writexl::write_xlsx(file.path(export_dir, "events_length_triple_sourceData.xlsx"))
  




# Save table

prop_pfs_tab <- proportions_l |>
  imap(~ .x |> add_column(event_type = .y, .before = 1)) |>
  map(~{
    test <- .x |>
        filter(is_detectable) |>
        select(n_frame, n_not_frame) |>
        chisq.test()

    .x |>
      add_column(pvalue = c(0,0, test[["p.value"]] ))

    }) |>
  bind_rows()

prop_pfs_tab$p_adj <- NA

prop_pfs_tab$p_adj[!is.na(prop_pfs_tab$pvalue)] = p.adjust(
  prop_pfs_tab$pvalue[!is.na(prop_pfs_tab$pvalue)]
)
prop_pfs_tab

# save table 2 essentially identical to source data Fig 6E

# prop_pfs_tab |>
#   mutate(
#     category = case_when(
#       ! is_detectable ~ "not measured",
#       has_ds ~ "DAS",
#       .default = "not DAS") |>
#       factor(levels = c("not measured", "not DAS", "DAS")),
#     event_type = case_match(
#       event_type,
#       "A3" ~ "Alt. 3' ss",
#       "A5" ~ "Alt. 5' ss",
#       "SE" ~ "Cassette exon"
#     ),
#     n_total = n_frame + n_not_frame) |>
#   select(`Event type` = event_type,
#          `Quantification` = category,
#          `Frame-preserving events` = n_frame,
#          `Total events` = n_total,
#          `% in frame` = prop_in_frame,
#          `p-value` = pvalue,
#          `p adjusted` = p_adj) |>
#   writexl::write_xlsx(file.path(export_dir, "events_length_triple.xlsx"))



#~ Position ----

# Are the frame-shifting events near extremities of the genes?



gcoords <- wb_load_gene_coords(289) |>
  select(gene_id,
         gstart = start,
         gend = end)




pos_SE <- coords_all$coords[[ which(coords_all$event_type == "SE") ]] |>
  mutate(gene_id = str_match(event_id, "^(WBGene[0-9]+);")[,2],
         .before = 1) |>
  mutate(length_is_triple = (exon_length %% 3) == 0) |>
  left_join(gcoords,
            by = "gene_id") |>
  mutate(distance_from_start = if_else(strand == "+",
                                    exon_start - gstart,
                                    gend - exon_end),
         distance_from_end = if_else(strand == "+",
                                     gend - exon_end,
                                     exon_start - gstart)) |>
  select(event_id, distance_from_start, distance_from_end, length_is_triple) |>
  full_join(
    features_long |>
      filter(event_type == "SE") |>
      select(event_type, event_id, has_ds, is_detectable) |>
      distinct(),
    by = "event_id"
  ) |>
  filter(event_id %in% se_exons_in_cds)
stopifnot( !any(is.na(pos_SE)) )




pos_A5 <- coords_all$coords[[ which(coords_all$event_type == "A5") ]] |>
  mutate(gene_id = str_match(event_id, "^(WBGene[0-9]+);")[,2],
         .before = 1) |>
  mutate(length_is_triple = (overhang_length %% 3) == 0) |>
  left_join(gcoords,
            by = "gene_id") |>
  mutate(distance_from_start = if_else(strand == "+",
                                       overhang_start - gstart,
                                       gend - overhang_end),
         distance_from_end = if_else(strand == "+",
                                     gend - overhang_end,
                                     overhang_start - gstart)) |>
  select(event_id, distance_from_start, distance_from_end, length_is_triple) |>
  full_join(
    features_long |>
      filter(event_type == "A5") |>
      select(event_type, event_id, has_ds, is_detectable) |>
      distinct(),
    by = "event_id"
  ) |>
  filter(event_id %in% a5_overhang_in_cds)
stopifnot( !any(is.na(pos_A5)) )



pos_A3 <- coords_all$coords[[ which(coords_all$event_type == "A3") ]] |>
  mutate(gene_id = str_match(event_id, "^(WBGene[0-9]+);")[,2],
         .before = 1) |>
  mutate(length_is_triple = (overhang_length %% 3) == 0) |>
  left_join(gcoords,
            by = "gene_id") |>
  mutate(distance_from_start = if_else(strand == "+",
                                       overhang_start - gstart,
                                       gend - overhang_end),
         distance_from_end = if_else(strand == "+",
                                     gend - overhang_end,
                                     overhang_start - gstart)) |>
  select(event_id, distance_from_start, distance_from_end, length_is_triple) |>
  full_join(
    features_long |>
      filter(event_type == "A3") |>
      select(event_type, event_id, has_ds, is_detectable) |>
      distinct(),
    by = "event_id"
  ) |>
  filter(event_id %in% a3_overhang_in_cds)
stopifnot( !any(is.na(pos_A3)) )





bind_rows(pos_SE, pos_A5, pos_A3) |>
  filter(is_detectable) |>
  mutate(distance_from_extremity = pmin(distance_from_start, distance_from_end),
         event_type = case_match(event_type,
                                 "A3" ~ "Alt. 3' ss",
                                 "A5" ~ "Alt. 5' ss",
                                 "SE" ~ "Cassette exon")) |>
  ggplot() +
  theme_classic() +
  scale_y_continuous(labels = scales::label_number(scale = 1e4, accuracy = 1)) +
  scale_fill_brewer(type = "qual",
                    labels = c("Frame-shifting","Frame-preserving")) +
  guides(fill = guide_legend(title = NULL )) +
  theme(legend.position = "top") +
  ylab("Density (×10,000)") +
  # scale_x_log10() +
  scale_x_continuous(limits = c(0, 5000)) +
  xlab("Distance from gene extremities (bp, truncated)") +
  facet_grid(rows = vars(event_type), scales = "free_y") +
  geom_density(aes(x = distance_from_extremity, fill = length_is_triple),
               alpha = .5)

# ggsave("events_length_triple_position.pdf", path = export_dir,
#        width = 16, height = 10, units = "cm")

# > Fig 6F

# save source data

# bind_rows(pos_SE, pos_A5, pos_A3) |>
#   mutate(distance_from_extremity = pmin(distance_from_start, distance_from_end)) |>
#   select(event_type, event_id, length_is_triple,
#          distance_from_gene_start = distance_from_start,
#          distance_from_gene_end = distance_from_end,
#          distance_from_extremity,
#          is_detectable,
#          has_ds) |>
#   writexl::write_xlsx(file.path(export_dir, "events_length_triple_position_sourceData.xlsx"))
  


pos_SE |>
  # filter(is_detectable) |>
  mutate(distance_from_extremity = pmin(distance_from_start, distance_from_end)) |>
  wilcox.test(distance_from_extremity ~ length_is_triple, data = _)

pos_A3 |>
  # filter(is_detectable) |>
  mutate(distance_from_extremity = pmin(distance_from_start, distance_from_end)) |>
  wilcox.test(distance_from_extremity ~ length_is_triple, data = _)

pos_A5 |>
  # filter(is_detectable) |>
  mutate(distance_from_extremity = pmin(distance_from_start, distance_from_end)) |>
  wilcox.test(distance_from_extremity ~ length_is_triple, data = _)

list(pos_SE, pos_A3, pos_A5) |>
  setNames(c("SE","A3","A5")) |>
  map_dbl(~{
    .x |>
      filter(is_detectable) |>
      mutate(distance_from_extremity = pmin(distance_from_start, distance_from_end)) |>
      wilcox.test(distance_from_extremity ~ length_is_triple, data = _) |>
      pluck("p.value")
  }) |>
  p.adjust()




# Microexons ----



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
            is_detectable = any(detectable),
            .by = c(exon_length, event_id))

lim_high <- 27
skipped_exons |>
  summarize(has_ds = any(is_ds),
            is_detectable = any(detectable),
            .by = c(exon_length, event_id)) |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = exon_length,
                     fill = has_ds),
                 bins = 30,
                 color = 'white') +
  scale_x_log10(labels = scales::label_comma()) +
  geom_vline(aes(xintercept = lim_high),
             color = 'grey50', linetype = "dashed", linewidth = .8) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  geom_text(data = tibble(
    x = c(11, 3000),
    y = c(90,90),
    label = c(paste0("Microexons\n(≤ 27 bp)\n",
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
#        width = 12, height = 14, units = "cm")


# export all cassette exons

# skipped_exons |>
#   summarize(has_ds = any(is_ds),
#             is_detectable = any(detectable),
#             .by = c(exon_length, event_id)) |>
#   separate_wider_regex(event_id,
#                        patterns = c(gene_id = "^WBGene[0-9]{8}", ";",
#                                     "SE", ":",
#                                     chr = "[IVX]+", ":",
#                                     start = "[0-9]+", "-[0-9]+:[0-9]+-",
#                                     end = "[0-9]+", ":", "[+-]$"),
#                        cols_remove = FALSE) |>
#   mutate(gene_name = i2s(gene_id, gids),
#          coordinates = paste0(chr,":",start,"-",end)) |>
#   mutate(is_microexon = exon_length <= 27) |>
#   select(gene_name, event_id, is_detectable, has_ds, exon_length, is_microexon, coordinates) |>
#   writexl::write_xlsx(file.path(export_dir, "microexons_sourceData.xlsx"))






# test

skp_ex_test <- skipped_exons |>
  summarize(has_ds = any(is_ds),
            is_detectable = any(detectable),
            .by = c(exon_length, event_id)) |>
  filter(is_detectable) |>
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
  mutate(type = if_else(is_microexon, "microexon", "longer exon") |>
           factor(levels = c("microexon", "longer exon"))) |>
  ggplot(aes(x = type, y = prop_ds, label = gene_name)) +
  theme_classic() +
  geom_point(position = pos, alpha = .2) +
  # ggrepel::geom_text_repel(box.padding = 1.2, position = pos, color = 'grey40') +
  scale_y_continuous(labels = scales::label_percent()) +
  xlab(NULL) + ylab("Proportion of neuron pairs DAS")

# save plot

# ggsave("microexons_ds_pairs.pdf", path = export_dir,
#        width = 7, height = 9, units = "cm")

# source data

# skipped_ex_by_pairs |>
#   mutate(is_microexon = exon_length <= 27,
#          prop_ds = nb_ds/nb_tot) |>
#   select(event_id, exon_length, is_microexon, nb_ds, nb_non_ds, nb_tot, prop_ds) |>
#   writexl::write_xlsx(file.path(export_dir, "microexon_pairs_ds_sourceData.xlsx"))



## Test

skipped_ex_by_pairs_filt <- skipped_ex_by_pairs |>
  filter(nb_tot > 10) |>
  mutate(is_microexon = exon_length <= 27,
         prop_ds = nb_ds/nb_tot)

skipped_ex_by_pairs_filt |>
  summarize(mean_prop_ds = mean(prop_ds),
            .by = is_microexon)

wilcox.test(prop_ds ~ is_microexon, data = skipped_ex_by_pairs_filt)





# Gini ----

dpsi_gini <- dpsi |>
  filter(event_type != "RI") |>
  filter(!is.na(dPSI)) |>
  summarize(max_dpsi = max(abs(dPSI)),
            gini = DescTools::Gini(abs(dPSI)),
            nb_pairs = n(),
            .by = c("event_id", "gene_id","event_type"))





#~ Plot Gini ----


dpsi_gini_plot <- dpsi_gini |>
  filter(nb_pairs > 5) |>
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

dpsi_gini_plot |>
  ggplot() +
  theme_classic() +
  ylab(expression(Gini~index~group("(",group("|",Delta*PSI,"|"),")"))) +
  xlab(expression(max~group("(",group("|",Delta*PSI,"|"),")"))) +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  theme(legend.position = "none") +
  annotate("segment", x = .5, xend = 1, y = .5, linetype = 'dashed', color = 'grey') +
  geom_vline(aes(xintercept = .5), linetype = 'dashed', color = 'grey') +
  facet_wrap(~event_type) +
  geom_point(aes(x = max_dpsi, y = gini, color = event_type),
             alpha = .3) +
  geom_text(data = dpsi_gini_plot |>
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

# ggsave("events_gini.pdf", path = export_dir,
#        width = 10, height = 8, units = "cm", scale = 1.8)

# souce data

# dpsi_gini |>
#   mutate(gene_name = i2s(gene_id, gids)) |>
#   select(event_type, gene_name, event_id, max_dpsi, gini_dpsi = gini, nb_measured_pairs = nb_pairs) |>
#   writexl::write_xlsx(file.path(export_dir, "events_gini_sourceData.xlsx"))
  

#~ Gini for microexons ----

# https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette#8197703
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


uex_gini <- dpsi_gini |>
  filter(event_type == "SE") |>
  left_join(skipped_exons_lengths,
            by = "event_id") |>
  mutate(type = if_else(exon_length <= lim_high,
                        "Microexon", "Longer exon"))

uex_gini_plot <- uex_gini |>
  filter(nb_pairs > 5) |>
  mutate(gene_name = i2s(gene_id, gids))

uex_gini |>
  ggplot() +
  theme_classic() +
  ylab(expression(Gini~index~group("(",group("|",Delta*PSI,"|"),")"))) +
  xlab(expression(max~group("(",group("|",Delta*PSI,"|"),")"))) +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  theme(legend.position = "none") +
  annotate("segment", x = .5, xend = 1, y = .5, linetype = 'dashed', color = 'grey') +
  geom_vline(aes(xintercept = .5), linetype = 'dashed', color = 'grey') +
  facet_wrap(~type) +
  geom_point(aes(x = max_dpsi, y = gini),
             color = gg_color_hue(6)[[5]],
             alpha = .3) +
  geom_text(data = uex_gini_plot |>
              summarize(low_dpsi = sum( max_dpsi <= .5),
                        low_gini = sum(max_dpsi > .5 & gini <= .5),
                        high_gini = sum(max_dpsi > .5 & gini > .5),
                        .by = type) |>
              pivot_longer(-type,
                           names_to = "class",
                           values_to = "number") |>
              mutate(prop = round(100*number / sum(number), 1),
                     .by = type) |>
              mutate(label = paste0(number, "\n(",prop,"%)")) |>
              left_join(tibble(x = c(0, 1, 1),
                               y = c(.1, .1, .9),
                               hjust = c(0,1,1),
                               class = c("low_dpsi","low_gini","high_gini"))),
            aes(x=x,y=y,label=label,hjust = hjust))


# ggsave("gini_microexons.pdf", path = export_dir,
#        width = 6, height = 4, units = "cm", scale = 2)

# ---> source data

# uex_gini |>
#   filter(nb_pairs > 5) |>
#   mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE),
#          .after = gene_id) |>
#   writexl::write_xlsx(file.path(export_dir, "sourceData_gini_microexons.xlsx"))















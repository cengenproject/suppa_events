


# Inits ----

library(tidyverse)


library(wbData)

tx2g <- wb_load_tx2gene(289)
gids <- wb_load_gene_ids(289)

export_dir <- "data/outs/240425_fig"






# Load ----


dpsi_dta <- read.table("data/240304_dpsi/240304.dpsi") |>
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
         neurB != "Ref")



# filter ----



# filter from thresholded
# gene_expression_table <-  read.delim("../majiq/data/2024-03-05_alec_integration/bsn12_subtracted_integrated_binarized_expression_withVDDD_FDR0.05_030424.tsv") |>
#   as.data.frame() |>
#   rownames_to_column("gene_id") |>
#   pivot_longer(-gene_id,
#                names_to = "neuron_id",
#                values_to = "expressed") |>
#   mutate(is_expressed = expressed == 1L) |> select(-expressed)

gene_expression_table <- cengenDataSC::cengen_sc_2_bulk |>
  as.data.frame() |>
  rownames_to_column("gene_id") |>
  pivot_longer(-gene_id,
               names_to = "neuron_id",
               values_to = "expressed") |>
  mutate(is_expressed = expressed == 1L) |> select(-expressed)


neurons_here <- unique(gene_expression_table$neuron_id)
genes_with_known_expr <- unique(gene_expression_table$gene_id)



# annotate expression

dpsi <- dpsi |>
  filter(neurA %in% neurons_here,
         neurB %in% neurons_here,
         gene_id %in% genes_with_known_expr) |>
  left_join(gene_expression_table |> rename(expr_in_neurA = is_expressed),
            by = c("gene_id", neurA = "neuron_id")) |>
  left_join(gene_expression_table |> rename(expr_in_neurB = is_expressed),
            by = c("gene_id", neurB = "neuron_id"))

dpsi <- dpsi |>
  mutate(detectable = expr_in_neurA & expr_in_neurB) |>
  select( -expr_in_neurA, -expr_in_neurB) |>
  mutate(is_ds = p.val < .1 & detectable & abs(dPSI) > .3)


# qs::qsave(dpsi, "intermediates/240425_dpsi/filt_dpsi.qs")




# psi <- read.delim("data/240301b_psiPerEvent.psi") |>
#   rownames_to_column("event_id") |>
#   as_tibble() |>
#   separate_wider_regex(event_id,
#                        patterns = c(gene_id = "^WBGene[0-9]{8}", ";",
#                                     event_type = "[SEA53MXRIFL]{2}", "\\:",
#                                     event_coordinates = "[IXV]+\\:[0-9:\\-]+:[+-]$"),
#                        cols_remove = FALSE)
# 
# 
# 
# 
# 
# 
# psi_lg <- pivot_longer(psi,
#                        -c(gene_id, event_type, event_coordinates, event_id),
#                        names_to = "sample_id",
#                        values_to = "PSI") |>
#   mutate(neuron_id = str_match(sample_id, "^([A-Zef0-9]{2,4})r[0-9]{1,4}")[,2])
# 
# 
# 
# tpm <- read.delim("data/231208_str_q_tx_TPM.tsv") |>
#   rownames_to_column("transcript_id") |>
#   as_tibble() |>
#   mutate(gene_id = wb_tx2g(transcript_id, tx2g, warn_missing = TRUE),
#          .after = "transcript_id")
# 
# 
# 
# stopifnot(all.equal(sort(unique(psi_lg$event_id)),
#                     sort(unique(dpsi$event_id))))
# 
# stopifnot(all.equal(
#   sort(unique(psi_lg$neuron_id)) |>
#     setdiff(c("ADF", "M4", "Ref")),
#   sort(unique(
#     union(dpsi$neurA,
#           dpsi$neurB))
#   ) |> setdiff("Ref")
# ))
# 
# 
# 
# #~ check a few random examples ----
# my_ev <- sample(psi_lg$event_id, 1)
# my_neurA <- sample(psi_lg$neuron_id |> setdiff(c("ADF", "M4")), 1)
# my_neurB <- sample(psi_lg$neuron_id |> setdiff(c("ADF", "M4")), 1)
# 
# 
# psi_lg |>
#   filter(event_id == my_ev,
#          neuron_id %in% c(my_neurA, my_neurB)) |>
#   summarize(mean_PSI = mean(PSI, na.rm = TRUE),
#             .by = neuron_id) |>
#   (\(x) print(x))() |>
#   pull(mean_PSI) |> diff()


my_ev <- sample(filt_dpsi$event_id[filt_dpsi$detectable & abs(filt_dpsi$dPSI) > .3], 1)
my_neurA <- sample(dpsi$neurA, 1)
my_neurB <- sample(dpsi$neurB, 1)

dpsi |>
  filter(event_id == my_ev,
         neurA %in% c(my_neurA, my_neurB),
         neurB %in% c(my_neurA, my_neurB)) |>
  mutate(detectable = expr_in_neurA & expr_in_neurB) |>
  select( -expr_in_neurA, -expr_in_neurB) |>
  mutate(is_ds = p.val < .1 & detectable & abs(dPSI) > .3) |>
  select(gene_name, event_type, neuron_pair, dPSI,p.val, detectable, is_ds)








#/ ======= Load ====== / ----

# Cleaned up version
# see below ("check a few events in genome browser"): low dPSI is often meaningless and in the noise

dpsi <- qs::qread("intermediates/240425_dpsi/filt_dpsi.qs")






# #~ PSI ----
# 
# psi_lg <- read.delim("data/240301b_psiPerEvent.psi") |>
#   rownames_to_column("event_id") |>
#   as_tibble() |>
#   separate_wider_regex(event_id,
#                        patterns = c(gene_id = "^WBGene[0-9]{8}", ";",
#                                     event_type = "[SEA53MXRIFL]{2}", "\\:",
#                                     event_coordinates = "[IXV]+\\:[0-9:\\-]+:[+-]$"),
#                        cols_remove = FALSE) |>
#   pivot_longer(-c(gene_id, event_type, event_coordinates, event_id),
#                        names_to = "sample_id",
#                        values_to = "PSI") |>
#   mutate(neuron_id = str_match(sample_id, "^([A-Zef0-9]{2,4})r[0-9]{1,4}")[,2]) |>
#   filter(neuron_id %in% neurons_here) |>
#   left_join(gene_expression_table |>
#               rename(expressed = is_expressed),
#             by = c("gene_id", "neuron_id"))





# Number DS ----
# note: already corrected for multiple testing (option -gc)

dpsi$p.val |> hist(breaks = 70)

table(dpsi$p.val < .05)

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

# ggsave("event_type_volcanoes.png", path = export_dir,
#        width = 30, height = 12, units = "cm")


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




# type vs DS
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

# ggsave("ds_per_type.pdf", path = export_dir,
#        width = 21, height = 8, units = "cm")



# type vs DS vs detectable


dpsi |>
  summarize(gene_expressed = any(detectable),
            has_ds = any(is_ds),
            .by = c(event_id, event_type)) |>
  mutate(
    category = case_when(
    !gene_expressed ~ "gene not expressed",
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
  geom_bar(aes(x = event_type, fill = category)) +
  scale_fill_manual(values = c("grey90", "grey30", "darkred")) +
  xlab(NULL) + ylab("Number of events") +
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  scale_y_continuous(labels = scales::label_comma())

# ggsave("ds_per_type.pdf", path = export_dir,
#        width = 16, height = 9, units = "cm")





# Check a few examples in browser

dpsi |>
  filter(p.val < .05 & detectable &
           dPSI > .3) |>
  slice_sample(n = 1) |>
  as.data.frame()

dpsi |>
  filter(p.val < .05 & detectable) |>
  pull(dPSI) |>
  (\(x) table(x > .3))()



dpsi |>
  filter(p.val < .05 & detectable) |> pull(dPSI) |> abs() |> hist()





# Extract coordinates of event ----
source("R/extract_event_coordinates.R")

coords_all <- dpsi |>
  select(event_id, event_type, gene_id, gene_name, event_coordinates) |>
  distinct() |>
  nest(.by = event_type) |>
  mutate(coords = map2(event_type, data,
                       ~ extract_coords(.x, .y[["event_coordinates"]]))) |>
  mutate(coords = map2(data, coords, ~bind_cols(.x["event_id"],
                                                .y)))
# qs::qsave(coords_all, "intermediates/240425_dpsi/240425_all_coords.qs")




# Check that equivalent to previous

# all.equal(d_a3,
#           d_a3b)
# 
# 
# all.equal(d_mx,
#           dpsi |>
#             filter(event_type == "MX") |>
#             left_join(coords_all$coords[[which(coords_all$event_type == "MX")]],
#                       by = "event_id") |>
#             select(all_of(names(d_mx))))






# GC content ----

# Import after running "sequence_properties.R"


seq_gc <- qs::qread("intermediates/240425_dpsi/seq_properties.qs") |>
  imap(~ {
    .x |>
      add_column(type = .y) |>
      pivot_longer(-c(event_id, type),
                   names_pattern = "^([a-z_]+)_(width|gc)$",
                   names_to = c("measurement", ".value")) |>
      mutate(percent_GC = 100 * gc/width)
  })



# Conservation ----

# Import after running "sequence_conservation.R"

seq_dir <- "intermediates/240305_event_coords/"
cons_files <- list.files(seq_dir, pattern = "*_cons.qs",
                        full.names = FALSE)
cons_types <- str_match(cons_files,
                       "^240305_([as35eflmxri]{2})_cons\\.qs$")[,2] |>
  toupper()


seq_cons <- cons_files |>
  set_names(cons_types) |>
  imap(\(.f, .t){
    qs::qread(file.path(seq_dir, .f)) |>
      add_column(type = .t)
  })






#>> prep plot ----

#~  A3 ----
d_a3 <- dpsi |>
  filter(event_type == "A3") |>
  left_join(coords_all$coords[[which(coords_all$event_type == "A3")]],
            by = "event_id") |>
  summarize(has_ds = any(is_ds),
            is_detectable = any(detectable),
            .by = c(intron_length, overhang_length,
                    event_id)) |>
  pivot_longer(-c(event_id, has_ds, is_detectable),
               names_to = "feature",
               names_pattern = "^([a-z]+)_length$",
               values_to = "length") |>
  left_join(seq_cons[["A3"]] |>
              pivot_longer(-c(event_id, type),
                           names_pattern = "^([a-z_]+)_conservation$",
                           names_to = "feature",
                           values_to = "conservation"),
            by = c("event_id", "feature")
  ) |> 
  left_join(seq_gc[["A3"]] |> rename(feature = measurement),
            by = c("event_id", "feature")
  )

#~ A5 ----

d_a5 <- dpsi |>
  filter(event_type == "A5") |>
  left_join(coords_all$coords[[which(coords_all$event_type == "A5")]],
            by = "event_id") |>
  summarize(has_ds = any(is_ds),
            is_detectable = any(detectable),
            .by = c(intron_length, overhang_length,
                    event_id)) |>
  pivot_longer(-c(event_id, has_ds, is_detectable),
               names_to = "feature",
               names_pattern = "^([a-z]+)_length$",
               values_to = "length") |>
  left_join(seq_cons[["A5"]] |>
              pivot_longer(-c(event_id, type),
                           names_pattern = "^([a-z_]+)_conservation$",
                           names_to = "feature",
                           values_to = "conservation"),
            by = c("event_id", "feature")
  ) |>
  left_join(seq_gc[["A5"]] |> rename(feature = measurement),
            by = c("event_id", "feature"))


#~ AF ----



d_af <- dpsi |>
  filter(event_type == "AF") |>
  left_join(coords_all$coords[[which(coords_all$event_type == "AF")]],
            by = "event_id") |>
  summarize(has_ds = any(is_ds),
            is_detectable = any(detectable),
            .by = c(event_id,
                    distal_exon_length, distal_intron_length,
                    proximal_exon_length, proximal_intron_length)) |>
  pivot_longer(-c(event_id, has_ds, is_detectable),
               names_to = "feature",
               names_pattern = "^([a-z_]+)_length$",
               values_to = "length") |>
  left_join(seq_cons[["AF"]] |>
              pivot_longer(-c(event_id, type),
                           names_pattern = "^([a-z_]+)_conservation$",
                           names_to = "feature",
                           values_to = "conservation"),
            by = c("event_id", "feature")
  ) |> 
  left_join(seq_gc[["AF"]] |> rename(feature = measurement),
            by = c("event_id", "feature"))



#~ AL ----


d_al <- dpsi |>
  filter(event_type == "AL") |>
  left_join(coords_all$coords[[which(coords_all$event_type == "AL")]],
            by = "event_id") |>
  summarize(has_ds = any(is_ds),
            is_detectable = any(detectable),
            .by = c(distal_exon_length, distal_intron_length,
                    proximal_exon_length, proximal_intron_length,
                    event_id)) |>
  pivot_longer(-c(event_id, has_ds, is_detectable),
               names_to = "feature",
               names_pattern = "^([a-z_]+)_length$",
               values_to = "length")  |>
  left_join(seq_cons[["AL"]] |>
              pivot_longer(-c(event_id, type),
                           names_pattern = "^([a-z_]+)_conservation$",
                           names_to = "feature",
                           values_to = "conservation"),
            by = c("event_id", "feature")
  ) |> 
  left_join(seq_gc[["AL"]] |> rename(feature = measurement),
            by = c("event_id", "feature"))

#~ MX ----


d_mx <- dpsi |>
  filter(event_type == "MX") |>
  left_join(coords_all$coords[[which(coords_all$event_type == "MX")]],
            by = "event_id") |>
  summarize(has_ds = any(is_ds),
            is_detectable = any(detectable),
            .by = c(first_exon_length, first_up_intron_length,
                    first_dn_intron_length, second_exon_length,
                    second_up_intron_length, second_dn_intron_length,
                    event_id)) |>
  pivot_longer(-c(event_id, has_ds, is_detectable),
               names_to = "feature",
               names_pattern = "^([a-z_]+)_length$",
               values_to = "length") |>
  mutate(feature = recode(feature,
                          first_up_intron = "first_upstream_intron",
                          first_dn_intron = "first_downstream_intron",
                          second_up_intron = "second_upstream_intron",
                          second_dn_intron = "second_downstream_intron")) |>
  left_join(seq_cons[["MX"]] |>
              pivot_longer(-c(event_id, type),
                           names_pattern = "^([a-z_]+)_conservation$",
                           names_to = "feature",
                           values_to = "conservation"),
            by = c("event_id", "feature")
  ) |>
  left_join(seq_gc[["MX"]] |> rename(feature = measurement),
            by = c("event_id", "feature"))



#~ RI ----


d_ri <- dpsi |>
  filter(event_type == "RI") |>
  left_join(coords_all$coords[[which(coords_all$event_type == "RI")]],
            by = "event_id") |>
  summarize(has_ds = any(is_ds),
            is_detectable = any(detectable),
            .by = c(intron_length, upstream_exon_length,
                    downstream_exon_length,
                    event_id)) |>
  pivot_longer(-c(event_id, has_ds, is_detectable),
               names_to = "feature",
               names_pattern = "^([a-z_]+)_length$",
               values_to = "length") |>
  left_join(seq_cons[["RI"]] |>
              pivot_longer(-c(event_id, type),
                           names_pattern = "^([a-z_]+)_conservation$",
                           names_to = "feature",
                           values_to = "conservation"),
            by = c("event_id", "feature")
  ) |> 
  left_join(seq_gc[["RI"]] |> rename(feature = measurement),
            by = c("event_id", "feature"))


#~ SE ----


d_se <- dpsi |>
  filter(event_type == "SE") |>
  left_join(coords_all$coords[[which(coords_all$event_type == "SE")]],
            by = "event_id") |>
  summarize(has_ds = any(is_ds),
            is_detectable = any(detectable),
            .by = c(upstream_intron_length, downstream_intron_length,
                    exon_length,
                    event_id)) |>
  pivot_longer(-c(event_id, has_ds, is_detectable),
               names_to = "feature",
               names_pattern = "^([a-z_]+)_length$",
               values_to = "length") |>
  left_join(seq_cons[["SE"]] |>
              pivot_longer(-c(event_id, type),
                           names_pattern = "^([a-z_]+)_conservation$",
                           names_to = "feature",
                           values_to = "conservation"),
            by = c("event_id", "feature")
  ) |> 
  left_join(seq_gc[["SE"]] |> rename(feature = measurement),
            by = c("event_id", "feature"))



# qs::qsave(list(d_a3, d_a5, d_af, d_al, d_mx, d_ri, d_se),
#           "intermediates/240307_events_full.qs")




# >>> Plots ----

#~~  A3 ----
gg_l_a3 <- d_a3 |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x = has_ds, y = length, color = has_ds)) +
  facet_wrap(~feature, scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Length (bp)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))


gg_gc_a3 <- d_a3 |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x = has_ds, y = percent_GC, color = has_ds)) +
  facet_wrap(~feature, scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Percent GC") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))

gg_cs_a3 <- d_a3 |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x = has_ds, y = conservation, color = has_ds)) +
  facet_wrap(~feature, scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Conservation score") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))




#~~ A5 ----
            
gg_l_a5 <- d_a5 |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = length, color = has_ds)) +
  facet_wrap(~feature, scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Length (bp)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))


gg_gc_a5 <- d_a5 |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = percent_GC, color = has_ds)) +
  facet_wrap(~feature, scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Percent GC") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))

gg_cs_a5 <- d_a5 |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = conservation, color = has_ds)) +
  facet_wrap(~feature, scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Conservation score") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))




#~~ AF ----
gg_l_af <- d_af |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = length, color = has_ds)) +
  facet_grid(cols = vars(feature), scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Length (bp)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))


gg_gc_af <- d_af |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = percent_GC, color = has_ds)) +
  facet_grid(cols = vars(feature), scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Percent GC") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))

gg_cs_af <- d_af |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = conservation, color = has_ds)) +
  facet_grid(cols = vars(feature), scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Conservation score") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))



#~~ AL ----

gg_l_al <- d_al |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = length, color = has_ds)) +
  facet_grid(cols = vars(feature), scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Length (bp)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))


gg_gc_al <- d_al |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = percent_GC, color = has_ds)) +
  facet_grid(cols = vars(feature), scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Percent GC") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))


gg_cs_al <- d_al |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = conservation, color = has_ds)) +
  facet_grid(cols = vars(feature), scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Conservation score") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))



#~~ MX ----

gg_l_mx <- d_mx |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = length, color = has_ds)) +
  facet_grid(cols = vars(feature), scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Length (bp)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))


gg_gc_mx <- d_mx |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = percent_GC, color = has_ds)) +
  facet_grid(cols = vars(feature), scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Percent GC") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))

gg_cs_mx <- d_mx |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = conservation, color = has_ds)) +
  facet_grid(cols = vars(feature), scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Conservation score") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))



#~~ RI ----

gg_l_ri <- d_ri |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = length, color = has_ds)) +
  facet_wrap(~feature, scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Length (bp)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))



gg_gc_ri <- d_ri |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = percent_GC, color = has_ds)) +
  facet_wrap(~feature, scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Percent GC") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))



gg_cs_ri <- d_ri |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = conservation, color = has_ds)) +
  facet_wrap(~feature, scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Conservation score") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))




#~~ SE ----

gg_l_se <- d_se |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = length, color = has_ds)) +
  facet_wrap(~feature, scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Length (bp)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))



gg_gc_se <- d_se |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = percent_GC, color = has_ds)) +
  facet_wrap(~feature, scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Percent GC") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))




gg_cs_se <- d_se |>
  filter(is_detectable) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = conservation, color = has_ds)) +
  facet_wrap(~feature, scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Conservation score") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))





#~ save ----


gr_a3 <- gridExtra::gtable_rbind(
  egg::gtable_frame(ggplotGrob(gg_l_a3)),
  egg::gtable_frame(ggplotGrob(gg_gc_a3)),
  egg::gtable_frame(ggplotGrob(gg_cs_a3))
)
gr_a5 <- gridExtra::gtable_rbind(
  egg::gtable_frame(ggplotGrob(gg_l_a5)),
  egg::gtable_frame(ggplotGrob(gg_gc_a5)),
  egg::gtable_frame(ggplotGrob(gg_cs_a5))
)
gr_af <- gridExtra::gtable_rbind(
  egg::gtable_frame(ggplotGrob(gg_l_af)),
  egg::gtable_frame(ggplotGrob(gg_gc_af)),
  egg::gtable_frame(ggplotGrob(gg_cs_af))
)
gr_al <- gridExtra::gtable_rbind(
  egg::gtable_frame(ggplotGrob(gg_l_al)),
  egg::gtable_frame(ggplotGrob(gg_gc_al)),
  egg::gtable_frame(ggplotGrob(gg_cs_al))
)
gr_mx <- gridExtra::gtable_rbind(
  egg::gtable_frame(ggplotGrob(gg_l_mx)),
  egg::gtable_frame(ggplotGrob(gg_gc_mx)),
  egg::gtable_frame(ggplotGrob(gg_cs_mx))
)
gr_ri <- gridExtra::gtable_rbind(
  egg::gtable_frame(ggplotGrob(gg_l_ri)),
  egg::gtable_frame(ggplotGrob(gg_gc_ri)),
  egg::gtable_frame(ggplotGrob(gg_cs_ri))
)
gr_se <- gridExtra::gtable_rbind(
  egg::gtable_frame(ggplotGrob(gg_l_se)),
  egg::gtable_frame(ggplotGrob(gg_gc_se)),
  egg::gtable_frame(ggplotGrob(gg_cs_se))
)

gr_row_up <- gridExtra::gtable_cbind(gr_a3, gr_a5, gr_af, gr_al)
gr_row_dn <- gridExtra::gtable_cbind(gr_mx, gr_ri, gr_se)



gr <-gridExtra::arrangeGrob(gr_row_up, gr_row_dn)

# ggsave("upper_row.png", plot = gr_row_up, path = export_dir,
#        width = 35, height = 14, units = "cm")
# 
# ggsave("lower_row.png", plot = gr_row_dn, path = export_dir,
#        width = 35, height = 14, units = "cm")


# ggsave("all_together.png", plot = gr, path = export_dir,
#        width = 40, height = 60, units = "cm")




#~~~ Individual plots

# ggsave("lengths_a3.pdf", plot = gg_l_a3, path = export_dir,
#        width = 9, height = 6, units = "cm")
# 
# ggsave("percGC_a3.pdf", plot = gg_gc_a3, path = export_dir,
#        width = 9, height = 6, units = "cm")
# 
# ggsave("cons_a3.pdf", plot = gg_cs_a3, path = export_dir,
#        width = 9, height = 6, units = "cm")
# 
# 
# 
# ggsave("lengths_a5.pdf", plot = gg_l_a5, path = export_dir,
#        width = 9, height = 6, units = "cm")
# 
# ggsave("percGC_a5.pdf", plot = gg_gc_a5, path = export_dir,
#        width = 9, height = 6, units = "cm")
# ggsave("cons_a5.pdf", plot = gg_cs_a5, path = export_dir,
#        width = 9, height = 6, units = "cm")
# 
# 
# ggsave("lengths_af.pdf", plot = gg_l_af, path = export_dir,
#        width = 18, height = 6, units = "cm")
# 
# ggsave("percGC_af.pdf", plot = gg_gc_af, path = export_dir,
#        width = 18, height = 6, units = "cm")
 # ggsave("cons_af.pdf", plot = gg_cs_af, path = export_dir,
 #       width = 18, height = 6, units = "cm")
 # 
# 
# 
# 
# ggsave("lengths_al.pdf", plot = gg_l_al, path = export_dir,
#        width = 18, height = 6, units = "cm")
# 
# ggsave("percGC_al.pdf", plot = gg_gc_al, path = export_dir,
#        width = 18, height = 6, units = "cm")
# ggsave("cons_al.pdf", plot = gg_cs_al, path = export_dir,
#        width = 18, height = 6, units = "cm")
# 
# 
# 
# ggsave("lengths_mx.pdf", plot = gg_l_mx, path = export_dir,
#        width = 18, height = 6, units = "cm")
# 
# ggsave("percGC_mx.pdf", plot = gg_gc_mx, path = export_dir,
#        width = 18, height = 6, units = "cm")
# ggsave("cons_mx.pdf", plot = gg_cs_mx, path = export_dir,
#        width = 18, height = 6, units = "cm")
# 
# 
# 
# ggsave("lengths_ri.pdf", plot = gg_l_ri, path = export_dir,
#        width = 12, height = 6, units = "cm")
# 
# ggsave("percGC_ri.pdf", plot = gg_gc_ri, path = export_dir,
#        width = 12, height = 6, units = "cm")
# ggsave("cons_ri.pdf", plot = gg_cs_ri, path = export_dir,
#        width = 12, height = 6, units = "cm")
# 
# 
# 
# 
# ggsave("lengths_se.pdf", plot = gg_l_se, path = export_dir,
#        width = 12, height = 6, units = "cm")
# 
# ggsave("percGC_se.pdf", plot = gg_gc_se, path = export_dir,
#        width = 12, height = 6, units = "cm")
# ggsave("cons_se.pdf", plot = gg_cs_se, path = export_dir,
#               width = 12, height = 6, units = "cm")
       
       
# Stat tests ----

tests <- map_dfr(list(d_a3, d_a5, d_af, d_al, d_mx, d_ri, d_se),
        ~ .x |>
          rename(event_type = type.x) |>
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
                       values_to = "p_val")) |>
  mutate(padj = p.adjust(p_val))



table(tests$padj < .1)

hist(tests$p_val, breaks = 40)
hist(tests$padj, breaks = 40)

tests |>
  filter(padj < 0.05) |>
  mutate(sig = cut(padj,
                   breaks =   c(.1,  .05,  .01, .001,  0),
                   labels = rev(c("#", "*", "**", "***"))))


medians <- map_dfr(
  list(d_a3, d_a5, d_af, d_al, d_mx, d_ri, d_se),
  
  ~ .x |>
    rename(event_type = type.x) |>
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
)


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




#~ Focus on some of the significant changes ----


# se exon length

d_se |>
  filter(is_detectable) |>
  filter(feature == "exon") |>
  ggplot() +
  theme_classic() +
  geom_density(aes(x = length, fill = has_ds),
               alpha = .5) +
  scale_x_log10() +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_y_continuous(limits = c(0,1.2)) +
  theme(legend.position = "none") +
  xlab("Skipped exon length (bp)")

# ggsave("se_exon_length_density.pdf", path = export_dir,
#        width = 15, height = 8, units = "cm")



# AF length

d_af |>
  filter(is_detectable) |>
  filter(feature == "distal_exon") |>
  ggplot() +
  theme_classic() +
  geom_density(aes(x = length, fill = has_ds),
               alpha = .5) +
  scale_x_log10() +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_y_continuous(limits = c(0,1.2)) +
  theme(legend.position = "none") +
  xlab("Alternative first exon, distal exon length (bp)")

# ggsave("af_distal_exon_length_density.pdf", path = export_dir,
#        width = 15, height = 8, units = "cm")


d_af |>
  filter(is_detectable) |>
  filter(feature == "distal_intron") |>
  ggplot() +
  theme_classic() +
  geom_density(aes(x = length, fill = has_ds),
               alpha = .5) +
  scale_x_log10(limits = c(100, 40000)) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_y_continuous(limits = c(0,1.2)) +
  theme(legend.position = "none") +
  xlab("Alternative first exon, distal intron length (bp)")

# ggsave("af_distal_intron_length_density.pdf", path = export_dir,
#        width = 15, height = 8, units = "cm")


d_af |>
  filter(is_detectable) |>
  filter(feature == "proximal_exon") |>
  ggplot() +
  theme_classic() +
  geom_density(aes(x = length, fill = has_ds),
               alpha = .5) +
  scale_x_log10() +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_y_continuous(limits = c(0,1.2)) +
  theme(legend.position = "none") +
  xlab("Alternative first exon, proximal exon length (bp)")

# ggsave("af_proximal_exon_length_density.pdf", path = export_dir,
#        width = 15, height = 8, units = "cm")




# AF GC content

d_af |>
  filter(is_detectable) |>
  filter(feature == "distal_intron") |>
  ggplot() +
  theme_classic() +
  geom_density(aes(x = percent_GC, fill = has_ds),
               alpha = .5) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  theme(legend.position = "none") +
  xlab("Alternative first exon, distal intron GC content (%)")

# ggsave("af_distal_intron_gc_density.pdf", path = export_dir,
#        width = 15, height = 8, units = "cm")


d_af |>
  filter(is_detectable) |>
  filter(feature == "proximal_intron") |>
  ggplot() +
  theme_classic() +
  geom_density(aes(x = percent_GC, fill = has_ds),
               alpha = .5) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  theme(legend.position = "none") +
  xlab("Alternative first exon, proximal intron GC content (%)")

# ggsave("af_proximal_intron_gc_density.pdf", path = export_dir,
#        width = 15, height = 8, units = "cm")




# AF conservation


d_af |>
  filter(is_detectable) |>
  filter(feature == "distal_intron") |>
  ggplot() +
  theme_classic() +
  geom_density(aes(x = conservation, fill = has_ds),
               alpha = .5, bounds = c(0,1)) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_x_continuous(limits = c(0,1)) +
  theme(legend.position = "none") +
  xlab("Alternative first exon, distal intron conservation score")

# ggsave("af_distal_intron_conservation_density.pdf", path = export_dir,
#        width = 15, height = 8, units = "cm")



d_af |>
  filter(is_detectable) |>
  filter(feature == "proximal_intron") |>
  ggplot() +
  theme_classic() +
  geom_density(aes(x = conservation, fill = has_ds),
               alpha = .5, bounds = c(0,1)) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_x_continuous(limits = c(0,1)) +
  theme(legend.position = "none") +
  xlab("Alternative first exon, proximal intron conservation score")

# ggsave("af_proximal_intron_conservation_density.pdf", path = export_dir,
#        width = 15, height = 8, units = "cm")













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
            .by = c(exon_length, event_id))

skipped_exons2 |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = exon_length,
                     fill = has_ds),
                 bins = 40,
                 color = 'white') +
  scale_x_log10() +
  geom_vline(aes(xintercept = 27),
             color = 'grey10', linetype = "dashed") +
  scale_fill_manual(values = c("grey30", "darkred")) +
  geom_text(data = tibble(
    x = c(11, 1000),
    y = c(60,60),
    label = c(paste0("Microexons (<= 27 bp)\n",
                     sum(skipped_exons2$has_ds[skipped_exons2$exon_length <= 27]),
                     "/",
                     sum(skipped_exons2$exon_length <= 27),
                     " (", round(100*mean(skipped_exons2$has_ds[skipped_exons2$exon_length <= 27]))," %)"),
              paste0("Other exons (> 27 bp)\n",
                     sum(skipped_exons2$has_ds[skipped_exons2$exon_length > 27]),
                     "/",
                     sum(skipped_exons2$exon_length > 27),
                     " (", round(100*mean(skipped_exons2$has_ds[skipped_exons2$exon_length > 27]))," %)"))
  ),
  aes(x=x,y=y,label=label)) +
  theme(legend.position = "none") +
  xlab("Exon length (bp)") +
  ylab("Number of exons")

# ggsave("microexons.pdf", path = export_dir,
#        width = 20, height = 12, units = "cm")

skipped_exons2 |>
  mutate(microexon = exon_length <= 30) |>
  # summarize(prop_ds = mean(has_ds),
  #           .by = microexon) |>
  ggplot() +
  theme_classic() +
  geom_bar(aes(x = microexon,
               fill = has_ds)) +
  scale_fill_manual(values = c("grey30", "darkred"))










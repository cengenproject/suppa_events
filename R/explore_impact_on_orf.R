

# Idea: for each event (focus on SE first), does it affect the ORF?



# [[ version 3 ]] ----
# Inits ----



library(GenomicFeatures)
library(plyranges)

library(tidyverse)

source("R/extract_event_coordinates.R")


#~ load ----
# | load events ----
dpsi <- qs::qread("intermediates/240918/240920_dpsi_neurons.qs")


coords_all <- dpsi |>
  select(event_id, event_type, gene_id, gene_name, event_coordinates) |>
  distinct() |>
  nest(.by = event_type) |>
  mutate(coords = map2(event_type, data,
                       ~ extract_coords(.x, .y[["event_coordinates"]]))) |>
  mutate(coords = map2(data, coords,
                       ~bind_cols(.x["event_id"],
                                  .y)))


# | load annotation ----

gtf_path <- wbData::wb_get_gtf_path(289)


txdb <- local(wbData::wb_load_TxDb(289))

cds <- GenomicFeatures::cds(txdb)




# match events to cds ----


# | Alt 3' ss ----


# match the overhang

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

length(a3_overhang_in_cds)




# | Alt 5' ss ----


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

length(a5_overhang_in_cds)



# | SE ----

# we get the coordinates of the exon (inclusion)

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

exons_in_cds <- join_overlap_inner_directed(se_inclusion, cds) |>
  as.data.frame() |>
  pull(event_id) |>
  unique()

length(exons_in_cds)




#~ Check impact ----


proportions_ndf <-  qs::qread("intermediates/240918/241011_features.qs") |>
  filter((event_type == "A5" & feature == "overhang" & (event_id %in% a5_overhang_in_cds)) |
           (event_type == "A3" & feature == "overhang" & (event_id %in% a3_overhang_in_cds)) |
           (event_type == "SE" & feature == "exon" & (event_id %in% exons_in_cds)) ) |>
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


xx <- proportions_l[[1]] |>
  filter(is_detectable) |>
  select(n_frame, n_not_frame)
prop.test(as.matrix(xx))
fisher.test(as.matrix(xx))

prop_pfs_tab <- proportions_l |>
  imap(~ .x |> add_column(event_type = .y, .before = 1)) |>
  map(~{
    test <- .x |>
      filter(is_detectable) |>
      select(n_frame, n_not_frame) |>
      chisq.test()
    
    .x |>
      add_column(pvalue = c(NA,NA, test[["p.value"]] ))
    
  }) |>
  bind_rows()

prop_pfs_tab$p_adj <- NA

prop_pfs_tab$p_adj[!is.na(prop_pfs_tab$pvalue)] = p.adjust(
  prop_pfs_tab$pvalue[!is.na(prop_pfs_tab$pvalue)]
)
prop_pfs_tab



#~ Position ----

# Are the frame-shifting events near extremities of the genes?



gcoords <- wbData::wb_load_gene_coords(289) |>
  select(gene_id,
         gstart = start,
         gend = end)




pos_SE <- coords_all$coords[[ which(coords_all$event_type == "SE") ]] |>
  mutate(gene_id = str_match(event_id, "^(WBGene[0-9]+);")[,2],
         .before = 1) |>
  filter(event_id %in% exons_in_cds) |>
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
      filter(event_type == "SE", event_id %in% exons_in_cds) |>
      select(event_type, event_id, has_ds, is_detectable) |>
      distinct(),
    by = "event_id"
  )
stopifnot( !any(is.na(pos_SE)) )




pos_A5 <- coords_all$coords[[ which(coords_all$event_type == "A5") ]] |>
  mutate(gene_id = str_match(event_id, "^(WBGene[0-9]+);")[,2],
         .before = 1) |>
  filter(event_id %in% a5_overhang_in_cds) |>
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
      filter(event_type == "A5", event_id %in% a5_overhang_in_cds) |>
      select(event_type, event_id, has_ds, is_detectable) |>
      distinct(),
    by = "event_id"
  )
stopifnot( !any(is.na(pos_A5)) )



pos_A3 <- coords_all$coords[[ which(coords_all$event_type == "A3") ]] |>
  mutate(gene_id = str_match(event_id, "^(WBGene[0-9]+);")[,2],
         .before = 1) |>
  filter(event_id %in% a3_overhang_in_cds) |>
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
      filter(event_type == "A3", event_id %in% a3_overhang_in_cds) |>
      select(event_type, event_id, has_ds, is_detectable) |>
      distinct(),
    by = "event_id"
  )
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



# Interaction isoform switch ----

# we ask if the alternative event is associated with an isoform switch



#~ annot impact ----



introns <- GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE) |>
  unlist()

introns$transcript_id <- names(introns)

introns$prot_id <- str_match(introns$transcript_id,
                             paste0(seq_id = "^[A-Z0-9cel_]+\\.t?[0-9]{1,4}",
                                    prot_id = "((?:[a-z])?)",
                                    tx_id = "(?:\\.[0-9]{1,2})?$"))[,2]
names(introns) <- NULL


exons <- GenomicFeatures::exonsBy(txdb, "tx", use.names = TRUE) |>
  unlist()

exons$transcript_id <- names(exons)

exons$prot_id <- str_match(exons$transcript_id,
                             paste0(seq_id = "^[A-Z0-9cel_]+\\.t?[0-9]{1,4}",
                                    prot_id = "((?:[a-z])?)",
                                    tx_id = "(?:\\.[0-9]{1,2})?$"))[,2]
names(exons) <- NULL






#~ Alt 3' ss ----


# we match the introns: we already have the coordinates of the short intron
# we get the coordinates of the long intron by adding the overhang


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


a3_short_intron <- a3_gr |>
  filter(feature == "intron")

a3_long_intron <- a3_gr |>
  group_by(event_id) |>
  reduce_ranges_directed()



a3_short_intron_prots <- dplyr::left_join(as.data.frame(a3_short_intron),
                                          as.data.frame(introns),
                                          by = c("seqnames", "start", "end", "width", "strand"),
                                          relationship = "many-to-many") |>
  summarize(prot_ids = list(unique(prot_id)),
            .by = event_id) |>
  as_tibble()


a3_long_intron_prots <- dplyr::left_join(as.data.frame(a3_long_intron),
                                         as.data.frame(introns),
                                         by = c("seqnames", "start", "end", "width", "strand"),
                                         relationship = "many-to-many") |>
  summarize(prot_ids = list(unique(prot_id)),
            .by = event_id) |>
  as_tibble()

stopifnot(nrow(a3_short_intron_prots) == nrow(a3_long_intron_prots) &&
            all.equal(a3_short_intron_prots$event_id, a3_long_intron_prots$event_id) )

a3_isoforms <- left_join(a3_short_intron_prots,
                         a3_long_intron_prots,
                         by = "event_id") |>
  summarize(iso_switch = length( intersect(unlist(prot_ids.x), unlist(prot_ids.y)) ) == 0,
            no_switch = identical(unlist(prot_ids.x), unlist(prot_ids.y)),
            .by = event_id) |>
  transmute(event_id,
            isoform = case_when(iso_switch ~ "switch",
                                no_switch ~ "no switch",
                                .default = "overlap"))


# # without transmute, check that no TRUE-TRUE
# table(a3_isoforms$iso_switch, a3_isoforms$no_switch)


table(a3_isoforms$isoform)

a3_isoforms |>
  count(isoform) |>
  mutate(`%` = round( 100 * n / sum(n) , digits = 1))









#~ Alt 5' ss ----

# same as A3

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


a5_short_intron <- a5_gr |>
  filter(feature == "intron")

a5_long_intron <- a5_gr |>
  group_by(event_id) |>
  reduce_ranges_directed()



a5_short_intron_prots <- dplyr::left_join(as.data.frame(a5_short_intron),
                                          as.data.frame(introns),
                                          by = c("seqnames", "start", "end", "width", "strand"),
                                          relationship = "many-to-many") |>
  summarize(prot_ids = list(unique(prot_id)),
            .by = event_id) |>
  as_tibble()


a5_long_intron_prots <- dplyr::left_join(as.data.frame(a5_long_intron),
                                         as.data.frame(introns),
                                         by = c("seqnames", "start", "end", "width", "strand"),
                                         relationship = "many-to-many") |>
  summarize(prot_ids = list(unique(prot_id)),
            .by = event_id) |>
  as_tibble()

stopifnot(nrow(a5_short_intron_prots) == nrow(a5_long_intron_prots) &&
            all.equal(a5_short_intron_prots$event_id, a5_long_intron_prots$event_id) )

a5_isoforms <- left_join(a5_short_intron_prots,
                         a5_long_intron_prots,
                         by = "event_id") |>
  summarize(iso_switch = length( intersect(unlist(prot_ids.x), unlist(prot_ids.y)) ) == 0,
            no_switch = identical(unlist(prot_ids.x), unlist(prot_ids.y)),
            .by = event_id) |>
  transmute(event_id,
            isoform = case_when(iso_switch ~ "switch",
                                no_switch ~ "no switch",
                                .default = "overlap"))

# table(a5_isoforms$iso_switch, a5_isoforms$no_switch)


table(a5_isoforms$isoform)

a5_isoforms |>
  count(isoform) |>
  mutate(`%` = round( 100 * n / sum(n) , digits = 1))













#~ SE ----

# we get the coordinates of the exon (inclusion)

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

se_skipping <- se_gr |>
  filter(feature == "intron")





se_inclusion_prots <- dplyr::left_join(as.data.frame(se_inclusion),
                                          as.data.frame(exons),
                                          by = c("seqnames", "start", "end", "width", "strand"),
                                          relationship = "many-to-many") |>
  summarize(prot_ids = list(unique(prot_id)),
            .by = event_id) |>
  as_tibble()


se_skipping_prots <- dplyr::left_join(as.data.frame(se_skipping),
                                         as.data.frame(introns),
                                         by = c("seqnames", "start", "end", "width", "strand"),
                                         relationship = "many-to-many") |>
  summarize(prot_ids = list(unique(prot_id)),
            .by = event_id) |>
  as_tibble()

stopifnot(nrow(se_inclusion_prots) == nrow(se_skipping_prots) &&
            all.equal(se_inclusion_prots$event_id, se_skipping_prots$event_id) )


se_isoforms <- left_join(se_inclusion_prots,
                         se_skipping_prots,
                         by = "event_id") |>
  summarize(iso_switch = length( intersect(unlist(prot_ids.x), unlist(prot_ids.y)) ) == 0,
            no_switch = identical(unlist(prot_ids.x), unlist(prot_ids.y)),
            .by = event_id) |>
  transmute(event_id,
            isoform = case_when(iso_switch ~ "switch",
                                no_switch ~ "no switch",
                                .default = "overlap"))

# table(se_isoforms$iso_switch, se_isoforms$no_switch)


table(se_isoforms$isoform)

se_isoforms |>
  count(isoform) |>
  mutate(`%` = round( 100 * n / sum(n) , digits = 1))


se_isoforms |>
  filter(isoform == "overlap") |>
  slice_sample(n = 4)






#~ Alt first ----

# using introns (note, could also use exons, but there are small differences, and the introns
# better reflect the measurement of junction-spanning reads; the differences are mostly due
# to exons that share their splice site (part of the event) but
# have a different end, e.g. ajm-1, atf-7, clh-3)


af_gr <- coords_all |>
  filter(event_type == "AF") |>
  chuck("coords", 1L) |>
  select(event_id, chr, strand,
         distal_intron_start, distal_intron_end,
         proximal_intron_start, proximal_intron_end) |>
  pivot_longer(cols = contains("_intron_"),
               names_to = c("feature", ".value"),
               names_pattern = "^(distal|proximal)_intron_(start|end)$") |>
  rename(seqnames = chr) |>
  as_granges()


af_distal <- af_gr |>
  filter(feature == "distal")

af_proximal <- af_gr |>
  filter(feature == "proximal")






af_distal_prots <- dplyr::left_join(as.data.frame(af_distal),
                                    as.data.frame(introns),
                                    by = c("seqnames", "start", "end", "width", "strand"),
                                    relationship = "many-to-many") |>
  summarize(prot_ids = list(unique(prot_id)),
            .by = event_id) |>
  as_tibble()


af_proximal_prots <- dplyr::left_join(as.data.frame(af_proximal),
                                      as.data.frame(introns),
                                      by = c("seqnames", "start", "end", "width", "strand"),
                                      relationship = "many-to-many") |>
  summarize(prot_ids = list(unique(prot_id)),
            .by = event_id) |>
  as_tibble()

stopifnot(nrow(af_distal_prots) == nrow(af_proximal_prots) &&
            all.equal(af_distal_prots$event_id, af_proximal_prots$event_id) )


af_isoforms <- left_join(af_distal_prots,
                         af_proximal_prots,
                         by = "event_id") |>
  summarize(iso_switch = length( intersect(unlist(prot_ids.x), unlist(prot_ids.y)) ) == 0,
            no_switch = identical(unlist(prot_ids.x), unlist(prot_ids.y)),
            .by = event_id) |>
  transmute(event_id,
            isoform = case_when(iso_switch ~ "switch",
                                no_switch ~ "no switch",
                                .default = "overlap"))

# table(af_isoforms$iso_switch, af_isoforms$no_switch)


table(af_isoforms$isoform)












#~ Alt last ----



al_gr <- coords_all |>
  filter(event_type == "AL") |>
  chuck("coords", 1L) |>
  select(event_id, chr, strand,
         distal_intron_start, distal_intron_end,
         proximal_intron_start, proximal_intron_end) |>
  pivot_longer(cols = contains("_intron_"),
               names_to = c("feature", ".value"),
               names_pattern = "^(distal|proximal)_intron_(start|end)$") |>
  rename(seqnames = chr) |>
  as_granges()


al_distal <- al_gr |>
  filter(feature == "distal")

al_proximal <- al_gr |>
  filter(feature == "proximal")






al_distal_prots <- dplyr::left_join(as.data.frame(al_distal),
                                    as.data.frame(introns),
                                    by = c("seqnames", "start", "end", "width", "strand"),
                                    relationship = "many-to-many") |>
  summarize(prot_ids = list(unique(prot_id)),
            .by = event_id) |>
  as_tibble()


al_proximal_prots <- dplyr::left_join(as.data.frame(al_proximal),
                                      as.data.frame(introns),
                                      by = c("seqnames", "start", "end", "width", "strand"),
                                      relationship = "many-to-many") |>
  summarize(prot_ids = list(unique(prot_id)),
            .by = event_id) |>
  as_tibble()

stopifnot(nrow(al_distal_prots) == nrow(al_proximal_prots) &&
            all.equal(al_distal_prots$event_id, al_proximal_prots$event_id) )


al_isoforms <- left_join(al_distal_prots,
                                 al_proximal_prots,
                                 by = "event_id") |>
  summarize(iso_switch = length( intersect(unlist(prot_ids.x), unlist(prot_ids.y)) ) == 0,
            no_switch = identical(unlist(prot_ids.x), unlist(prot_ids.y)),
            .by = event_id) |>
  transmute(event_id,
            isoform = case_when(iso_switch ~ "switch",
                                no_switch ~ "no switch",
                                .default = "overlap"))

# table(al_isoforms$iso_switch, al_isoforms$no_switch)


table(al_isoforms$isoform)




#~ Check ----


event_impact <- bind_rows(
  a3_isoforms |> add_column(event_type = "A3", .before = 1),
  a5_isoforms |> add_column(event_type = "A5", .before = 1),
  se_isoforms |> add_column(event_type = "SE", .before = 1),
  af_isoforms |> add_column(event_type = "AF", .before = 1),
  al_isoforms |> add_column(event_type = "AL", .before = 1)
)

table(event_impact$event_type, event_impact$isoform)

events_DAS <- qs::qread("intermediates/240918/241011_features.qs") |>
  select(event_type, event_id, is_detectable, has_ds) |>
  filter(event_type != "MX") |>
  distinct()

#~~ by DAS ----
# display all numbers
left_join(events_DAS,
          event_impact,
          by = c("event_type", "event_id")) |>
  mutate(
    category = case_when(
      ! is_detectable ~ "not measured",
      has_ds ~ "DAS",
      .default = "not DAS") |>
      factor(levels = c("not measured", "not DAS", "DAS") |> rev()),
    event_type = case_match(event_type,
                            "A3" ~ "Alt. 3' ss",
                            "A5" ~ "Alt. 5' ss",
                            "SE" ~ "Cassette exon",
                            "AF" ~ "Alt. first exon",
                            "AL" ~ "Alt. last exon")
  ) |>
  count(event_type, category, isoform) |>
  ggplot() +
  theme_classic() +
  facet_wrap(~event_type) +
  geom_col(aes(x = category, y = n, fill = isoform),
           position = "dodge") +
  geom_text(aes( x = category, y = n+20, label = n, group = isoform),
            position = position_dodge(width = .9))



left_join(events_DAS,
          event_impact,
          by = c("event_type", "event_id")) |>
  mutate(
    category = case_when(
      ! is_detectable ~ "not measured",
      has_ds ~ "DAS",
      .default = "not DAS") |>
      factor(levels = c("not measured", "not DAS", "DAS") |> rev()),
    event_type = case_match(event_type,
                            "A3" ~ "Alt. 3' ss",
                            "A5" ~ "Alt. 5' ss",
                            "SE" ~ "Cassette exon",
                            "AF" ~ "Alt. first exon",
                            "AL" ~ "Alt. last exon")
  ) |>
  filter(category != "not measured",
         isoform != "overlap") |>
  count(event_type, category, isoform) |>
  ggplot() +
  theme_classic() +
  facet_wrap(~event_type) +
  geom_col(aes(x = category, y = n, fill = isoform),
           position = "fill")




#~~ by switch ----
# display all numbers
left_join(events_DAS,
          event_impact,
          by = c("event_type", "event_id")) |>
  mutate(
    category = case_when(
      ! is_detectable ~ "not measured",
      has_ds ~ "DAS",
      .default = "not DAS") |>
      factor(levels = c("not measured", "not DAS", "DAS") |> rev()),
    event_type = case_match(event_type,
                            "A3" ~ "Alt. 3' ss",
                            "A5" ~ "Alt. 5' ss",
                            "SE" ~ "Cassette exon",
                            "AF" ~ "Alt. first exon",
                            "AL" ~ "Alt. last exon")
  ) |>
  count(event_type, category, isoform) |>
  ggplot() +
  theme_classic() +
  facet_wrap(~event_type) +
  geom_col(aes( x = isoform, y = n, fill = category),
           position = "dodge") +
  geom_text(aes( x = isoform, y = n+20, label = n, group = category),
            position = position_dodge(width = .9))

left_join(events_DAS,
          event_impact,
          by = c("event_type", "event_id")) |>
  mutate(
    category = case_when(
      ! is_detectable ~ "not measured",
      has_ds ~ "DAS",
      .default = "not DAS") |>
      factor(levels = c("not measured", "not DAS", "DAS") |> rev()),
    event_type = case_match(event_type,
                            "A3" ~ "Alt. 3' ss",
                            "A5" ~ "Alt. 5' ss",
                            "SE" ~ "Cassette exon",
                            "AF" ~ "Alt. first exon",
                            "AL" ~ "Alt. last exon")
  ) |>
  filter(category != "not measured",
         isoform != "overlap") |>
  count(event_type, category, isoform) |>
  ggplot() +
  theme_classic() +
  coord_flip() +
  facet_grid(rows = vars(event_type)) +
  geom_col(aes( x = isoform, y = n, fill = category),
           position = "fill") +
  geom_text(aes( x = isoform, y = n, label = n, group = category),
            position = position_fill(vjust = .5))


left_join(events_DAS,
          event_impact,
          by = c("event_type", "event_id")) |>
  mutate(
    category = case_when(
      ! is_detectable ~ "not measured",
      has_ds ~ "DAS",
      .default = "not DAS") |>
      factor(levels = c("not measured", "not DAS", "DAS") |> rev()),
    event_type = case_match(event_type,
                            "A3" ~ "Alt. 3' ss",
                            "A5" ~ "Alt. 5' ss",
                            "SE" ~ "Cassette exon",
                            "AF" ~ "Alt. first exon",
                            "AL" ~ "Alt. last exon")
  ) |>
  filter(category != "not measured",
         isoform != "overlap") |>
  count(event_type, category, isoform) |>
  group_by(event_type) |>
  group_split() |>
  map_dfr(~{
    tibble(
      event_type = .x[["event_type"]][[1]],
      prop_das_noswitch = .x[["n"]][[2]] / ( .x[["n"]][[4]] + .x[["n"]][[2]] ),
      prop_das_switch = .x[["n"]][[1]] / ( .x[["n"]][[3]] + .x[["n"]][[1]] ),
      pval = (.x |>
                pull(n) |>
                matrix(nrow = 2) |>
                fisher.test() |>
                chuck("p.value") )
    )
  }) |>
  mutate(p_adj = p.adjust(pval))



#~~ only switch ----

# display grouping switch and overlap
assembled <- left_join(events_DAS,
          event_impact,
          by = c("event_type", "event_id")) |>
  mutate(
    category = case_when(
      ! is_detectable ~ "not measured",
      has_ds ~ "DAS",
      .default = "not DAS") |>
      factor(levels = c("not measured", "not DAS", "DAS") |> rev()),
    
    isoform = case_match(isoform,
                         "switch" ~ "switch",
                         c("no switch", "overlap") ~ "other"),
    
    event_type = case_match(event_type,
                            "A3" ~ "Alt. 3' ss",
                            "A5" ~ "Alt. 5' ss",
                            "SE" ~ "Cassette exon",
                            "AF" ~ "Alt. first exon",
                            "AL" ~ "Alt. last exon")
  ) |>
  filter(category != "not measured") |>
  count(event_type, category, isoform)

bind_rows(
  assembled |> filter(isoform == "switch"),
  assembled |>
    summarize(n = sum(n),
              .by = -c(isoform, n)) |>
    add_column(isoform = "total")
) |>
  ggplot() +
  theme_classic() +
  coord_flip() +
  facet_grid(rows = vars(event_type)) +
  geom_col(aes( x = isoform, y = n, fill = category),
           position = "fill") +
  geom_text(aes( x = isoform, y = n, label = n, group = category),
            position = position_fill(vjust = .5))



bind_rows(
  assembled |> filter(isoform == "switch"),
  assembled |>
    summarize(n = sum(n),
              .by = -c(isoform, n)) |>
    add_column(isoform = "total")
) |>
  group_by(event_type) |>
  group_split() |>
  map_dfr(~{
    tibble(
      event_type = .x[["event_type"]][[1]],
      prop_das_switch = .x[["n"]][[1]] / ( .x[["n"]][[2]] + .x[["n"]][[1]] ),
      prop_das_total = .x[["n"]][[3]] / ( .x[["n"]][[3]] + .x[["n"]][[4]] ),
      pval = (.x |>
                pull(n) |>
                matrix(nrow = 2) |>
                t() |>
                prop.test() |>
                chuck("p.value") )
    )
    
  }) |>
  mutate(p_adj = p.adjust(pval))



# ______ ----

proportions_ndf <-  qs::qread("intermediates/240918/241011_features.qs") |>
  filter((event_type == "A5" & feature == "overhang" & (event_id %in% a5_overhang_in_cds)) |
           (event_type == "A3" & feature == "overhang" & (event_id %in% a3_overhang_in_cds)) |
           (event_type == "SE" & feature == "exon" & (event_id %in% exons_in_cds)) ) |>
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


xx <- proportions_l[[1]] |>
  filter(is_detectable) |>
  select(n_frame, n_not_frame)
prop.test(as.matrix(xx))
fisher.test(as.matrix(xx))

prop_pfs_tab <- proportions_l |>
  imap(~ .x |> add_column(event_type = .y, .before = 1)) |>
  map(~{
    test <- .x |>
      filter(is_detectable) |>
      select(n_frame, n_not_frame) |>
      chisq.test()
    
    .x |>
      add_column(pvalue = c(NA,NA, test[["p.value"]] ))
    
  }) |>
  bind_rows()

prop_pfs_tab$p_adj <- NA

prop_pfs_tab$p_adj[!is.na(prop_pfs_tab$pvalue)] = p.adjust(
  prop_pfs_tab$pvalue[!is.na(prop_pfs_tab$pvalue)]
)
prop_pfs_tab






#~~ Check impact ---




proportions_ndf <- features_long |>
  filter(!prot_same) |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  summarize(n_frame = sum(length_is_triple),
            n_not_frame = sum(!length_is_triple),
            .by = c(has_ds, is_detectable)) |>
  mutate(prop_in_frame = round( 100 * n_frame / (n_frame + n_not_frame) )) |>
  arrange(is_detectable, has_ds)


proportions_ndf


proportions_ndf |>
  select(-prop_in_frame) |>
  pivot_longer(-c(has_ds, is_detectable),
               names_to = "frame",
               values_to = "count") |>
  mutate(
    category = case_when(
      ! is_detectable ~ "not measured",
      has_ds ~ "DAS",
      .default = "not DAS") |>
      factor(levels = c("not measured", "not DAS", "DAS") |> rev()),
    frame = factor(frame, levels = c("n_not_frame", "n_frame"))
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
  geom_col(aes(x = category, y = count, fill = frame))



#~ test

proportions_ndf |>
  filter(is_detectable) |>
  select(n_frame, n_not_frame) |>
  chisq.test()











# _____________ ----

#~~ separate ----
features_long <- qs::qread("intermediates/240918/241011_features.qs") |>
  filter( event_type == "A3" & feature == "overhang" )


proportions_ndf <- features_long |>
  filter(event_id %in% a3_overhang_in_cds) |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  summarize(n_frame = sum(length_is_triple),
            n_not_frame = sum(!length_is_triple),
            .by = c(has_ds, is_detectable)) |>
  mutate(prop_in_frame = round( 100 * n_frame / (n_frame + n_not_frame) )) |>
  arrange(is_detectable, has_ds)

proportions_ndf


proportions_ndf |>
  select(-prop_in_frame) |>
  pivot_longer(-c(has_ds, is_detectable),
               names_to = "frame",
               values_to = "count") |>
  mutate(
    category = case_when(
      ! is_detectable ~ "not measured",
      has_ds ~ "DAS",
      .default = "not DAS") |>
      factor(levels = c("not measured", "not DAS", "DAS") |> rev()),
    frame = factor(frame, levels = c("n_not_frame", "n_frame"))
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
  geom_col(aes(x = category, y = count, fill = frame))



#~~ test ----

proportions_ndf |>
  filter(is_detectable) |>
  select(n_frame, n_not_frame) |>
  chisq.test()







# [[ first version ]] ----

# Inits ----


library(tidyverse)

dpsi <- qs::qread("intermediates/240918/240920_dpsi_neurons.qs")


source("R/extract_event_coordinates.R")


gtf_path <- wbData::wb_get_gtf_path(289)


#~ load ----
coords_all <- dpsi |>
  select(event_id, event_type, gene_id, gene_name, event_coordinates) |>
  distinct() |>
  nest(.by = event_type) |>
  mutate(coords = map2(event_type, data,
                       ~ extract_coords(.x, .y[["event_coordinates"]]))) |>
  mutate(coords = map2(data, coords,
                       ~bind_cols(.x["event_id"],
                                  .y)))




#~ Preprocess sets of coordinates ----

coords_exon <- rtracklayer::readGFF(gtf_path,
                                    filter = list(type = "exon")) |>
  rename(exon_start = start,
         exon_end = end) |>
  as_tibble() |>
  transmute(gene_id,
            transcript_id,
            exon_3p_end = if_else(strand == "+", exon_end, exon_start),
            exon_5p_start = if_else(strand == "+", exon_start, exon_end))




# # for alternative approach
# coords_cds <- rtracklayer::readGFF(gtf_path,
#                                    filter = list(type = "CDS")) |>
#   rename(exon_start = start,
#          exon_end = end) |>
#   as_tibble()
# 
# coords_cds_match <- coords_cds |>
#   select(gene_id, exon_start, exon_end) |>
#   add_column(type = "CDS") |>
#   distinct()
# 
# 
# 
# coords_start_codon <- rtracklayer::readGFF(gtf_path,
#                                    filter = list(type = "start_codon")) |>
#   as_tibble()
# 
# coords_start_codon_match <- coords_start_codon |>
#   select(gene_id, codon_start = start, codon_end = end) |>
#   add_column(type = "start_codon") |>
#   distinct()
# 
# 
# 
# coords_stop_codon <- rtracklayer::readGFF(gtf_path,
#                                            filter = list(type = "stop_codon")) |>
#   as_tibble()
# 
# coords_stop_codon_match <- coords_stop_codon |>
#   select(gene_id, codon_start = start, codon_end = end) |>
#   add_column(type = "stop_codon") |>
#   distinct()
# 
# 
# 
# 
# 
# coords_utr <- rtracklayer::readGFF(gtf_path,
#                                    filter = list(type = c("three_prime_utr", "five_prime_utr"))) |>
#   as_tibble()
# 
# coords_utr_match <- coords_utr |>
#   select(gene_id, strand, utr_start = start, utr_end = end) |>
#   add_column(type = "UTR") |>
#   distinct()
# 
# coords_utr_intron <- coords_utr |>
#   arrange(transcript_id, start, end) |>
#   mutate(next_start = if_else( strand == "+", lead(start), lag(start) ),
#          next_end = if_else( strand == "+", lead(end), lag(end) ),
#          .by = c(transcript_id, type)) |>
#   filter(!is.na(next_start)) |>
#   transmute(
#     gene_id,
#     intron_start = if_else( strand == "+", end + 1L, next_end + 1L ),
#     intron_end = if_else( strand == "+", next_start - 1L, start - 1L ),
#     intron_3p_end = if_else( strand == "+", intron_end, intron_start ),
#     intron_5p_start = if_else( strand == "+", intron_start, intron_end)
#   )










# SE ----




event_coords <- coords_all |>
  filter(event_type == "SE") |>
  chuck("coords", 1L) |>
  select(-matches("^c[1-4]$")) |>
  separate_wider_delim(event_id,
                       delim = ";",
                       names = c("gene_id", NA),
                       cols_remove = FALSE)

# some events share the same exon with different flanking introns
exon_coords <- event_coords |>
  select(gene_id, strand, exon_start, exon_end) |>
  distinct()





# Four main cases:
# 1. exon matches a CDS exon exactly
# 2. exon matches a UTR exon exactly
# 3/4. one side matches CDS, the other UTR (i.e. contains a stop or start codon)
# We are loosing 40 events:
# * many have a non-coding transcript (e.g. cyl-1)
# * some exons end within a stop codon (e.g. inx-10)

exons_in_cds <- exon_coords |>
  inner_join(coords_cds_match,
            by = c("gene_id",
                   "exon_start", "exon_end"))

exons_in_utr <- exon_coords |>
  inner_join(coords_utr_match,
             by = c("gene_id",
                    "exon_start" = "utr_start",
                    "exon_end" = "utr_end",
                    "strand"))



# just looking for start/stop codon is not enough: we also need a match of UTR
# so we need a start codon + upstream utr intron
# or a stop codon + downstream utr intron
# note that here up/downstream = 5p/3p in tx coordinates, not genomic coordinates

exons_on_start_flanked_by_utr_intron <- coords_start_codon_match |>
  inner_join(exon_coords,
             by = join_by(gene_id, within(codon_start, codon_end, exon_start, exon_end))) |>
  mutate(exon_5p_start_min1 = if_else(strand == "+", exon_start - 1L, exon_end + 1L)) |>
  semi_join(coords_utr_intron,
            by = c("exon_5p_start_min1" = "intron_3p_end"))


exons_on_stop_flanked_by_utr_intron <- coords_stop_codon_match |>
  inner_join(exon_coords,
             by = join_by(gene_id, within(codon_start, codon_end, exon_start, exon_end))) |>
  mutate(exon_3p_end_plus1 = if_else(strand == "+", exon_end + 1L, exon_start - 1L)) |>
  semi_join(coords_utr_intron,
            by = c("exon_3p_end_plus1" = "intron_5p_start"))





all_matched_exons <- bind_rows(
  exons_in_cds,
  exons_in_utr,
  exons_on_start_flanked_by_utr_intron,
  exons_on_stop_flanked_by_utr_intron
) |>
  select(gene_id, exon_start, exon_end, type)

stopifnot(nrow(anti_join(all_matched_exons, exon_coords)) == 0)


# exons we are not considering (e.g. non-coding transcripts):
exons_ignored <- anti_join(exon_coords, all_matched_exons,
                           by = join_by(gene_id, exon_start, exon_end))



# in some cases several exons with same coordinates but different coding potential
# from RNA-Seq we can't tell which transcript is affected
exons_with_unclear_impact <- all_matched_exons |>
  count(gene_id, exon_start, exon_end) |> 
  filter(n > 1) |> arrange(desc(n)) |>
  select(gene_id, exon_start, exon_end)


exons_with_impact <- all_matched_exons |>
  anti_join(exons_with_unclear_impact,
            by = join_by(gene_id, exon_start, exon_end))

stopifnot( nrow(exon_coords) ==
             nrow(exons_with_impact) + nrow(exons_ignored) + nrow(exons_with_unclear_impact) )

table(exons_with_impact$type)


event_types <- event_coords |>
  left_join(exons_with_impact,
            by = c("gene_id", "exon_start", "exon_end")) |>
  select(event_id, exon_length, type)

table(event_types$type)


#~ Check impact ----
features_long <- qs::qread("intermediates/240918/241011_features.qs") |>
  filter( event_type == "SE" & feature == "exon" ) |>
  left_join(event_types,
          by = c("event_id", "length" = "exon_length"))


proportions_ndf <- features_long |>
  # filter(type == "CDS") |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  summarize(n_frame = sum(length_is_triple),
            n_not_frame = sum(!length_is_triple),
            .by = c(has_ds, is_detectable)) |>
  mutate(prop_in_frame = round( 100 * n_frame / (n_frame + n_not_frame) )) |>
  arrange(is_detectable, has_ds)



proportions_ndf |>
  select(-prop_in_frame) |>
  pivot_longer(-c(has_ds, is_detectable),
               names_to = "frame",
               values_to = "count") |>
  mutate(
    category = case_when(
      ! is_detectable ~ "not measured",
      has_ds ~ "DAS",
      .default = "not DAS") |>
      factor(levels = c("not measured", "not DAS", "DAS") |> rev()),
    frame = factor(frame, levels = c("n_not_frame", "n_frame"))
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
  geom_col(aes(x = category, y = count, fill = frame))



#~ test ----

proportions_ndf |>
  filter(is_detectable) |>
  select(n_frame, n_not_frame) |>
  chisq.test()












# A3ss ----




#~ annotate impact on protein ----


a3_coords <- coords_all |>
  filter(event_type == "A3") |>
  chuck("coords", 1L) |>
  select(-matches("^c[1-4]$")) |>
  separate_wider_delim(event_id,
                       delim = ";",
                       names = c("gene_id", NA),
                       cols_remove = FALSE)

a3_coords_match <- a3_coords |>
  mutate(shorter_exon_5p_start = if_else(strand == "+", overhang_end + 1L, overhang_start - 1L),
         longer_exon_5p_start = if_else(strand == "+", overhang_start, overhang_end)) |>
  select(event_id, gene_id,
         ends_with("exon_5p_start")) |>
  pivot_longer(cols = ends_with("exon_5p_start"),
               names_to = c("exon_selected", NA,NA,NA),
               names_sep = "_",
               values_to = "exon_5p_start")






a3_corresponding_prot <- left_join(a3_coords_match,
                                   coords_exon,
                                   by = join_by(gene_id, exon_5p_start),
                                   relationship = "many-to-many") |>
  separate_wider_regex(transcript_id,
                       patterns = c(seq_id = "^[A-Z0-9cel_]+\\.t?[0-9]{1,4}",
                                    prot_id = "(?:[a-z])?",
                                    tx_id = "(?:\\.[0-9]{1,2})?$"),
                       cols_remove = FALSE) |>
  select(event_id, exon_selected, prot_id) |>
  distinct() |>
  summarize(prot_ids = list(sort(prot_id)),
            .by = c(event_id, exon_selected)) |>
  summarize(iso_switch = length( intersect(prot_ids[[1]], prot_ids[[2]]) ) == 0,
            no_switch = identical(prot_ids[[1]], prot_ids[[2]]),
            n = length(prot_ids),
            .by = event_id)

stopifnot(all(a3_corresponding_prot$n == 2))


table(a3_corresponding_prot$prot_same)


a3_isoform <- a3_corresponding_prot |>
  transmute(event_id,
            isoform = case_when(iso_switch ~ "switch",
                                no_switch ~ "no switch",
                                .default = "overlap"))

table(a3_isoform$isoform)

a3_isoform |>
  count(isoform) |>
  mutate(`%` = round( 100 * n / sum(n) , digits = 1))


#~~ alternative approach ----

# # While the first approach above uses exon names, here we match the coordinates to CDS or UTR.
# # This approach appears more "correct", but turns out less conservative:
# # for example see X:15,072,953..15,073,153 where we have an overhang that's in a UTR,
# # but affects which protein isoform is expressed
# 
# 
# a3_coords_events <- coords_all |>
#   filter(event_type == "A3") |>
#   chuck("coords",1L) |>
#   select(-matches("^c[1-4]$")) |>
#   separate_wider_delim(event_id,
#                        delim = ";",
#                        names = c("gene_id", NA),
#                        cols_remove = FALSE)
# 
# # some events share the same overhang with different flanking exons
# a3_coords_overhangs <- a3_coords_events |>
#   select(gene_id, overhang_start, overhang_end) |>
#   distinct()
# 
# 
# 
# 
# 
# # Two cases:
# # 1. overhang inside a CDS exon
# # 2. overhang inside UTR exon
# 
# 
# overhangs_in_cds <- a3_coords_overhangs |>
#   inner_join(coords_cds_match,
#              by = join_by(gene_id,
#                           between(overhang_start, exon_start, exon_end),
#                           between(overhang_end, exon_start, exon_end)))
# 
# overhangs_in_utr <- a3_coords_overhangs |>
#   inner_join(coords_utr_match,
#              by = join_by(gene_id,
#                           between(overhang_start, utr_start, utr_end),
#                           between(overhang_end, utr_start, utr_end)))
# 
# overhang_with_stop_codon <- coords_stop_codon_match |>
#   inner_join(a3_coords_overhangs,
#              by = join_by(gene_id,
#                           within(codon_start, codon_end, overhang_start, overhang_end)))
# 
# overhang_with_start_codon <- coords_start_codon_match |>
#   inner_join(a3_coords_overhangs,
#              by = join_by(gene_id,
#                           within(codon_start, codon_end, overhang_start, overhang_end)))
# 
# 
# 
# 
# all_matched_overhangs <- bind_rows(
#   overhangs_in_cds,
#   overhangs_in_utr,
#   overhang_with_stop_codon,
#   overhang_with_start_codon
#   ) |>
#   select(gene_id, overhang_start, overhang_end, type)
# 
# stopifnot(nrow(anti_join(all_matched_overhangs, a3_coords_overhangs)) == 0)
# 
# 
# # a3 we are not considering (e.g. non-coding transcripts):
# a3_ignored <- anti_join(a3_coords_overhangs, all_matched_overhangs,
#                            by = join_by(gene_id, overhang_start, overhang_end))
# 
# 
# 
# # in some cases several overhangs with same coordinates but different coding potential
# # from RNA-Seq we can't tell which transcript is affected
# overhangs_with_unclear_impact <- all_matched_overhangs |>
#   count(gene_id, overhang_start, overhang_end) |> 
#   filter(n > 1) |> arrange(desc(n)) |>
#   select(gene_id, overhang_start, overhang_end)
# 
# 
# overhangs_with_impact <- all_matched_overhangs |>
#   anti_join(overhangs_with_unclear_impact,
#             by = join_by(gene_id, overhang_start, overhang_end))
# 
# stopifnot( nrow(a3_coords_overhangs) ==
#              nrow(overhangs_with_impact) + nrow(a3_ignored) + nrow(overhangs_with_unclear_impact) )
# 
# table(overhangs_with_impact$type)
# 
# 
# a3_event_types <- a3_coords_events |>
#   left_join(overhangs_with_impact,
#             by = c("gene_id", "overhang_start", "overhang_end")) |>
#   select(event_id, overhang_length, type)
# 
# table(a3_event_types$type)
# 
# 
# # Compare apporaches
# 
# comp <- left_join(
#   a3_corresponding_prot,
#   a3_event_types,
#   by = "event_id"
# )
# 
# table(comp$prot_same, comp$type, useNA = 'ifany')
# 
# left_join(
#   a3_corresponding_prot,
#   a3_event_types,
#   by = "event_id"
# ) |>
#   filter(!prot_same & type == "UTR")
# 
# 
# 
# left_join(a3_coords_match,
#           coords_exon,
#           by = join_by(gene_id, exon_5p_start),
#           relationship = "many-to-many") |>
#   separate_wider_regex(transcript_id,
#                        patterns = c(seq_id = "^[A-Z0-9cel_]+\\.t?[0-9]{1,4}",
#                                     prot_id = "(?:[a-z])?",
#                                     tx_id = "(?:\\.[0-9]{1,2})?$"),
#                        cols_remove = FALSE) |>
#   select(event_id, exon_selected, prot_id) |>
#   distinct() |>
#   summarize(prot_ids = paste(sort(prot_id), collapse = ","),
#             .by = c(event_id, exon_selected)) |>
#   filter(event_id == "WBGene00004781;A3:III:6451774-6451821:6451766-6451821:-")
# 
# 
# 
# left_join(
#   a3_corresponding_prot,
#   a3_event_types,
#   by = "event_id"
# ) |>
#   filter(startsWith(event_id, "WBGene00004781"))
# 






#~ Check impact ----
features_long <- qs::qread("intermediates/240918/241011_features.qs") |>
  filter( event_type == "A3" & feature == "overhang" ) |>
  left_join(a3_corresponding_prot |> select(-n),
            by = c("event_id"))


proportions_ndf <- features_long |>
  filter(iso_switch) |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  summarize(n_frame = sum(length_is_triple),
            n_not_frame = sum(!length_is_triple),
            .by = c(has_ds, is_detectable)) |>
  mutate(prop_in_frame = round( 100 * n_frame / (n_frame + n_not_frame) )) |>
  arrange(is_detectable, has_ds)

proportions_ndf


proportions_ndf |>
  select(-prop_in_frame) |>
  pivot_longer(-c(has_ds, is_detectable),
               names_to = "frame",
               values_to = "count") |>
  mutate(
    category = case_when(
      ! is_detectable ~ "not measured",
      has_ds ~ "DAS",
      .default = "not DAS") |>
      factor(levels = c("not measured", "not DAS", "DAS") |> rev()),
    frame = factor(frame, levels = c("n_not_frame", "n_frame"))
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
  geom_col(aes(x = category, y = count, fill = frame))



#~~ test ----

proportions_ndf |>
  filter(is_detectable) |>
  select(n_frame, n_not_frame) |>
  chisq.test()





#~ check impact vs DAS ----

# previous comparison was on length of overhang (multiple of 3) within the protein
# instead ask if "changing the protein" and "being DAS" are linked



a3_das_prot <- qs::qread("intermediates/240918/241011_features.qs") |>
  filter( event_type == "A3" ) |>
  select(event_id, has_ds, is_detectable) |>
  distinct() |>
  left_join(a3_isoform,
            by = c("event_id"))



proportions_ndf <- a3_das_prot |>
  summarize(n_no_switch = sum(isoform == "no switch"),
            n_switch = sum(isoform == "switch"),
            .by = c(has_ds, is_detectable)) |>
  mutate(prop_no_switch = round( 100 * n_no_switch / (n_switch + n_no_switch) )) |>
  arrange(is_detectable, has_ds)

proportions_ndf


proportions_ndf |>
  select(-prop_no_switch) |>
  pivot_longer(-c(has_ds, is_detectable),
               names_to = "isoform",
               values_to = "count") |>
  mutate(
    category = case_when(
      ! is_detectable ~ "not measured",
      has_ds ~ "DAS",
      .default = "not DAS") |>
      factor(levels = c("not measured", "not DAS", "DAS") |> rev()),
    isoform = factor(isoform, levels = c("n_switch", "n_no_switch"))
  ) |>
  ggplot() +
  theme_classic() +
  guides(fill = guide_legend(title = NULL )) +
  theme(legend.position = "top") +
  scale_fill_brewer(type = "qual",
                    labels = c("isoform-switching", "No isoform change")) +
  xlab(NULL) +
  ylab("Number of events") +
  coord_flip() +
  geom_col(aes(x = category, y = count, fill = isoform))



#~ test ----

proportions_ndf |>
  filter(is_detectable) |>
  select(n_no_switch, n_switch) |>
  chisq.test()



#~ transposed ----


a3_das_prot |>
  summarize(not_measured = sum(! is_detectable),
            not_DAS = sum(is_detectable & !has_ds),
            DAS = sum(has_ds),
            .by = c(isoform)) |>
  mutate(`% DAS` = round( 100 * DAS / (DAS + not_DAS) )) |>
  arrange(isoform)



a3_das_prot |>
  summarize(not_measured = sum(! is_detectable),
            not_DAS = sum(is_detectable & !has_ds),
            DAS = sum(has_ds),
            .by = c(isoform)) |>
  pivot_longer(-isoform,
               names_to = "DAS",
               values_to = "n") |>
  filter(DAS != "not_measured") |>
  ggplot() +
  theme_classic() +
  guides(fill = guide_legend(title = NULL )) +
  theme(legend.position = "top") +
  scale_fill_brewer(type = "qual") +
  xlab(NULL) +
  ylab("Number of events") +
  coord_flip() +
  geom_col(aes(x = isoform, y = n, fill = DAS))





#~ test ----

a3_das_prot |>
  summarize(not_measured = sum(! is_detectable),
            not_DAS = sum(is_detectable & !has_ds),
            DAS = sum(has_ds),
            .by = c(isoform)) |>
  filter(isoform != "overlap") |>
  select(not_DAS, DAS) |>
  chisq.test()








# A5ss ----

#~ annotate impact on protein ----


a5_coords <- coords_all |>
  filter(event_type == "A5") |>
  chuck("coords", 1L) |>
  select(-matches("^c[1-4]$")) |>
  separate_wider_delim(event_id,
                       delim = ";",
                       names = c("gene_id", NA),
                       cols_remove = FALSE)

a5_coords_match <- a5_coords |>
  mutate(shorter_exon_3p_end = if_else(strand == "+", overhang_start - 1L, overhang_end + 1L),
         longer_exon_3p_end = if_else(strand == "+", overhang_end, overhang_start)) |>
  select(event_id, gene_id,
         ends_with("exon_3p_end")) |>
  pivot_longer(cols = ends_with("exon_3p_end"),
               names_to = c("exon_selected", NA,NA,NA),
               names_sep = "_",
               values_to = "exon_3p_end")






a5_corresponding_prot <- left_join(a5_coords_match,
                                   coords_exon,
                                   by = join_by(gene_id, exon_3p_end),
                                   relationship = "many-to-many") |>
  separate_wider_regex(transcript_id,
                       patterns = c(seq_id = "^[A-Z0-9cel_]+\\.t?[0-9]{1,4}",
                                    prot_id = "(?:[a-z])?",
                                    tx_id = "(?:\\.[0-9]{1,2})?$"),
                       cols_remove = FALSE) |>
  select(event_id, exon_selected, prot_id) |>
  distinct() |>
  summarize(prot_ids = list(sort(prot_id)),
            .by = c(event_id, exon_selected)) |>
  summarize(prot_same = length( intersect(prot_ids[[1]], prot_ids[[2]]) ) > 0,
            n = length(prot_ids),
            .by = event_id)

stopifnot(all(a5_corresponding_prot$n == 2))


table(a5_corresponding_prot$prot_same)




# ~~ alternative approach ----
# a5_event_coords <- coords_all |>
#   filter(event_type == "A5") |>
#   chuck("coords",1L) |>
#   select(-matches("^c[1-4]$")) |>
#   separate_wider_delim(event_id,
#                        delim = ";",
#                        names = c("gene_id", NA),
#                        cols_remove = FALSE)
# 
# # some events share the same overhang with different flanking exons
# a5_coords_overhang <- a5_event_coords |>
#   select(gene_id, overhang_start, overhang_end) |>
#   distinct()
# 
# 
# 
# 
# 
# # Two cases:
# # 1. overhang inside a CDS exon
# # 2. overhang inside UTR exon
# 
# 
# overhangs_in_cds <- a5_coords_overhang |>
#   inner_join(coords_cds_match,
#              by = join_by(gene_id,
#                           between(overhang_start, exon_start, exon_end),
#                           between(overhang_end, exon_start, exon_end)))
# 
# overhangs_in_utr <- a5_coords_overhang |>
#   inner_join(coords_utr_match,
#              by = join_by(gene_id,
#                           between(overhang_start, utr_start, utr_end),
#                           between(overhang_end, utr_start, utr_end)))
# 
# overhang_with_stop_codon <- coords_stop_codon_match |>
#   inner_join(a5_coords_overhang,
#              by = join_by(gene_id,
#                           within(codon_start, codon_end, overhang_start, overhang_end)))
# 
# 
# overhang_with_start_codon <- coords_start_codon_match |>
#   inner_join(a5_coords_overhang,
#              by = join_by(gene_id,
#                           within(codon_start, codon_end, overhang_start, overhang_end)))
# 
# 
# 
# 
# 
# all_matched_overhangs <- bind_rows(
#   overhangs_in_cds,
#   overhangs_in_utr,
#   overhang_with_stop_codon,
#   overhang_with_start_codon
# ) |>
#   select(gene_id, overhang_start, overhang_end, type)
# 
# stopifnot(nrow(anti_join(all_matched_overhangs, a5_coords_overhang)) == 0)
# 
# 
# # a5 we are not considering (e.g. non-coding transcripts):
# a5_ignored <- anti_join(a5_coords_overhang, all_matched_overhangs,
#                         by = join_by(gene_id, overhang_start, overhang_end))
# 
# 
# 
# # in some cases several overhangs with same coordinates but different coding potential
# # from RNA-Seq we can't tell which transcript is affected
# overhangs_with_unclear_impact <- all_matched_overhangs |>
#   count(gene_id, overhang_start, overhang_end) |> 
#   filter(n > 1) |> arrange(desc(n)) |>
#   select(gene_id, overhang_start, overhang_end)
# 
# 
# overhangs_with_impact <- all_matched_overhangs |>
#   anti_join(overhangs_with_unclear_impact,
#             by = join_by(gene_id, overhang_start, overhang_end))
# 
# stopifnot( nrow(a5_coords_overhang) ==
#              nrow(overhangs_with_impact) + nrow(a5_ignored) + nrow(overhangs_with_unclear_impact) )
# 
# table(overhangs_with_impact$type)
# 
# 
# a5_event_types <- a5_event_coords |>
#   left_join(overhangs_with_impact,
#             by = c("gene_id", "overhang_start", "overhang_end")) |>
#   select(event_id, overhang_length, type)
# 
# table(a5_event_types$type)
# 
# 
# 
# # Compare apporaches
# 
# comp <- left_join(
#   a5_corresponding_prot,
#   a5_event_types,
#   by = "event_id"
# )
# 
# table(comp$prot_same, comp$type, useNA = 'ifany')
# 
# left_join(
#   a5_corresponding_prot,
#   a5_event_types,
#   by = "event_id"
# ) |>
#   filter(!prot_same & type == "UTR")
# 
# 
# 
# left_join(a5_coords_match,
#           coords_exon,
#           by = join_by(gene_id, exon_3p_end),
#           relationship = "many-to-many") |>
#   separate_wider_regex(transcript_id,
#                        patterns = c(seq_id = "^[A-Z0-9cel_]+\\.t?[0-9]{1,4}",
#                                     prot_id = "(?:[a-z])?",
#                                     tx_id = "(?:\\.[0-9]{1,2})?$"),
#                        cols_remove = FALSE) |>
#   select(event_id, exon_selected, prot_id) |>
#   distinct() |>
#   summarize(prot_ids = paste(sort(prot_id), collapse = ","),
#             .by = c(event_id, exon_selected)) |>
#   filter(event_id == "WBGene00001518;A5:III:6319606-6320319:6319456-6320319:+")
# 
# 
# 
# left_join(
#   a3_corresponding_prot,
#   a3_event_types,
#   by = "event_id"
# ) |>
#   filter(startsWith(event_id, "WBGene00004781"))







#~ Check impact ----
features_long <- qs::qread("intermediates/240918/241011_features.qs") |>
  filter( event_type == "A5" & feature == "overhang" ) |>
  left_join(a5_corresponding_prot |> select(-n),
            by = c("event_id"))


proportions_ndf <- features_long |>
  filter(!prot_same) |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  summarize(n_frame = sum(length_is_triple),
            n_not_frame = sum(!length_is_triple),
            .by = c(has_ds, is_detectable)) |>
  mutate(prop_in_frame = round( 100 * n_frame / (n_frame + n_not_frame) )) |>
  arrange(is_detectable, has_ds)


proportions_ndf


proportions_ndf |>
  select(-prop_in_frame) |>
  pivot_longer(-c(has_ds, is_detectable),
               names_to = "frame",
               values_to = "count") |>
  mutate(
    category = case_when(
      ! is_detectable ~ "not measured",
      has_ds ~ "DAS",
      .default = "not DAS") |>
      factor(levels = c("not measured", "not DAS", "DAS") |> rev()),
    frame = factor(frame, levels = c("n_not_frame", "n_frame"))
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
  geom_col(aes(x = category, y = count, fill = frame))



#~ test ----

proportions_ndf |>
  filter(is_detectable) |>
  select(n_frame, n_not_frame) |>
  chisq.test()













#~~~~~~~~ ----



# Idea: for each event (focus on SE first), does it affect the ORF?


# Inits ----

library(GenomicFeatures)
library(plyranges)

library(tidyverse)



dpsi <- qs::qread("intermediates/240918/240920_dpsi_neurons.qs")


source("R/extract_event_coordinates.R")


gtf_path <- wbData::wb_get_gtf_path(289)


#~ load ----
coords_all <- dpsi |>
  select(event_id, event_type, gene_id, gene_name, event_coordinates) |>
  distinct() |>
  nest(.by = event_type) |>
  mutate(coords = map2(event_type, data,
                       ~ extract_coords(.x, .y[["event_coordinates"]]))) |>
  mutate(coords = map2(data, coords,
                       ~bind_cols(.x["event_id"],
                                  .y)))




#~ Preprocess sets of coordinates ----

# coords_exon <- rtracklayer::readGFF(gtf_path,
#                                     filter = list(type = "exon")) |>
#   rename(exon_start = start,
#          exon_end = end) |>
#   as_tibble() |>
#   transmute(gene_id,
#             transcript_id,
#             exon_3p_end = if_else(strand == "+", exon_end, exon_start),
#             exon_5p_start = if_else(strand == "+", exon_start, exon_end))


txdb <- local(wbData::wb_load_TxDb(289))

cds <- GenomicFeatures::cds(txdb)




# ~~~~~~~~~~~~ ----

# [[ Older approaches ]] ----


# SE ----

#~ annotate impact on protein ----
event_coords <- coords_all |>
  filter(event_type == "SE") |>
  chuck("coords", 1L) |>
  select(-matches("^c[1-4]$")) |>
  separate_wider_delim(event_id,
                       delim = ";",
                       names = c("gene_id", NA),
                       cols_remove = FALSE)

# some events share the same exon with different flanking introns
exon_coords <- event_coords |>
  select(gene_id, strand, exon_start, exon_end) |>
  distinct()





se_coords <- coords_all |>
  filter(event_type == "SE") |>
  chuck("coords", 1L) |>
  select(-matches("^c[1-4]$")) |>
  separate_wider_delim(event_id,
                       delim = ";",
                       names = c("gene_id", NA),
                       cols_remove = FALSE)

se_coords |>
  filter(gene_id == "WBGene00000018") |> as.data.frame()


se_coords_match <- se_coords |>
  mutate(upstream_exon_3p_end = if_else(strand == "+", upstream_intron_start - 1L, downstream_intron_end + 1L),
         cassette_exon_5p_start = if_else(strand == "+", exon_start, exon_end),
         cassette_exon_3p_end = if_else(strand == "+", exon_end, exon_start),
         downstream_exon_5p_start = if_else(strand == "+", downstream_intron_end + 1L, upstream_intron_start - 1L)) |>
  select(event_id, gene_id,
         matches("exon_[35]p_")) |>
  pivot_longer(cols = -c(event_id, gene_id),
               names_to = c("exon_selected", NA,NA,NA),
               names_sep = "_",
               values_to = "exon_5p_start")










#~ alt approach ----


coords_cds <- rtracklayer::readGFF(gtf_path,
                                   filter = list(type = "CDS")) |>
  rename(exon_start = start,
         exon_end = end) |>
  as_tibble()

coords_cds_match <- coords_cds |>
  select(gene_id, exon_start, exon_end) |>
  add_column(type = "CDS") |>
  distinct()

coords_utr <- rtracklayer::readGFF(gtf_path,
                                   filter = list(type = c("three_prime_utr", "five_prime_utr"))) |>
  as_tibble()
coords_utr_match <- coords_utr |>
  select(gene_id, strand, utr_start = start, utr_end = end) |>
  add_column(type = "UTR") |>
  distinct()

coords_start_codon <- rtracklayer::readGFF(gtf_path,
                                           filter = list(type = "start_codon")) |>
  as_tibble()
coords_start_codon_match <- coords_start_codon |>
  select(gene_id, codon_start = start, codon_end = end) |>
  add_column(type = "start_codon") |>
  distinct()
coords_stop_codon <- rtracklayer::readGFF(gtf_path,
                                          filter = list(type = "stop_codon")) |>
  as_tibble()
coords_stop_codon_match <- coords_stop_codon |>
  select(gene_id, codon_start = start, codon_end = end) |>
  add_column(type = "stop_codon") |>
  distinct()

coords_utr_intron <- coords_utr |>
  arrange(transcript_id, start, end) |>
  mutate(next_start = if_else( strand == "+", lead(start), lag(start) ),
         next_end = if_else( strand == "+", lead(end), lag(end) ),
         .by = c(transcript_id, type)) |>
  filter(!is.na(next_start)) |>
  transmute(
    gene_id,
    intron_start = if_else( strand == "+", end + 1L, next_end + 1L ),
    intron_end = if_else( strand == "+", next_start - 1L, start - 1L ),
    intron_3p_end = if_else( strand == "+", intron_end, intron_start ),
    intron_5p_start = if_else( strand == "+", intron_start, intron_end)
  )


# Four main cases:
# 1. exon matches a CDS exon exactly
# 2. exon matches a UTR exon exactly
# 3/4. one side matches CDS, the other UTR (i.e. contains a stop or start codon)
# We are loosing 40 events:
# * many have a non-coding transcript (e.g. cyl-1)
# * some exons end within a stop codon (e.g. inx-10)

exons_in_cds <- exon_coords |>
  inner_join(coords_cds_match,
             by = c("gene_id",
                    "exon_start", "exon_end"))

exons_in_utr <- exon_coords |>
  inner_join(coords_utr_match,
             by = c("gene_id",
                    "exon_start" = "utr_start",
                    "exon_end" = "utr_end",
                    "strand"))

# just looking for start/stop codon is not enough: we also need a match of UTR
# so we need a start codon + upstream utr intron
# or a stop codon + downstream utr intron
# note that here up/downstream = 5p/3p in tx coordinates, not genomic coordinates

exons_on_start_flanked_by_utr_intron <- coords_start_codon_match |>
  inner_join(exon_coords,
             by = join_by(gene_id, within(codon_start, codon_end, exon_start, exon_end))) |>
  mutate(exon_5p_start_min1 = if_else(strand == "+", exon_start - 1L, exon_end + 1L)) |>
  semi_join(coords_utr_intron,
            by = c("exon_5p_start_min1" = "intron_3p_end"))


exons_on_stop_flanked_by_utr_intron <- coords_stop_codon_match |>
  inner_join(exon_coords,
             by = join_by(gene_id, within(codon_start, codon_end, exon_start, exon_end))) |>
  mutate(exon_3p_end_plus1 = if_else(strand == "+", exon_end + 1L, exon_start - 1L)) |>
  semi_join(coords_utr_intron,
            by = c("exon_3p_end_plus1" = "intron_5p_start"))





all_matched_exons <- bind_rows(
  exons_in_cds,
  exons_in_utr,
  exons_on_start_flanked_by_utr_intron,
  exons_on_stop_flanked_by_utr_intron
) |>
  select(gene_id, exon_start, exon_end, type)

stopifnot(nrow(anti_join(all_matched_exons, exon_coords)) == 0)


# exons we are not considering (e.g. non-coding transcripts):
exons_ignored <- anti_join(exon_coords, all_matched_exons,
                           by = join_by(gene_id, exon_start, exon_end))



# in some cases several exons with same coordinates but different coding potential
# from RNA-Seq we can't tell which transcript is affected
exons_with_unclear_impact <- all_matched_exons |>
  count(gene_id, exon_start, exon_end) |> 
  filter(n > 1) |> arrange(desc(n)) |>
  select(gene_id, exon_start, exon_end)


exons_with_impact <- all_matched_exons |>
  anti_join(exons_with_unclear_impact,
            by = join_by(gene_id, exon_start, exon_end))

stopifnot( nrow(exon_coords) ==
             nrow(exons_with_impact) + nrow(exons_ignored) + nrow(exons_with_unclear_impact) )

table(exons_with_impact$type)


event_types <- event_coords |>
  left_join(exons_with_impact,
            by = c("gene_id", "exon_start", "exon_end")) |>
  select(event_id, exon_length, type)

table(event_types$type)


event_types |>
  filter(type != "CDS")


#~ based on CDS within ----









#~ Check impact ----
features_long <- qs::qread("intermediates/240918/241011_features.qs") |>
  filter( event_type == "SE" & feature == "exon" ) |>
  left_join(event_types,
            by = c("event_id", "length" = "exon_length"))


proportions_ndf <- features_long |>
  # filter(type == "CDS") |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  summarize(n_frame = sum(length_is_triple),
            n_not_frame = sum(!length_is_triple),
            .by = c(has_ds, is_detectable)) |>
  mutate(prop_in_frame = round( 100 * n_frame / (n_frame + n_not_frame) )) |>
  arrange(is_detectable, has_ds)



proportions_ndf |>
  select(-prop_in_frame) |>
  pivot_longer(-c(has_ds, is_detectable),
               names_to = "frame",
               values_to = "count") |>
  mutate(
    category = case_when(
      ! is_detectable ~ "not measured",
      has_ds ~ "DAS",
      .default = "not DAS") |>
      factor(levels = c("not measured", "not DAS", "DAS") |> rev()),
    frame = factor(frame, levels = c("n_not_frame", "n_frame"))
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
  geom_col(aes(x = category, y = count, fill = frame))



#~ test ----

proportions_ndf |>
  filter(is_detectable) |>
  select(n_frame, n_not_frame) |>
  chisq.test()












# A3ss ----




#~ annotate impact on protein ----


a3_coords <- coords_all |>
  filter(event_type == "A3") |>
  chuck("coords", 1L) |>
  select(-matches("^c[1-4]$")) |>
  separate_wider_delim(event_id,
                       delim = ";",
                       names = c("gene_id", NA),
                       cols_remove = FALSE)

a3_coords_match <- a3_coords |>
  mutate(shorter_exon_5p_start = if_else(strand == "+", overhang_end + 1L, overhang_start - 1L),
         longer_exon_5p_start = if_else(strand == "+", overhang_start, overhang_end)) |>
  select(event_id, gene_id,
         ends_with("exon_5p_start")) |>
  pivot_longer(cols = ends_with("exon_5p_start"),
               names_to = c("exon_selected", NA,NA,NA),
               names_sep = "_",
               values_to = "exon_5p_start")






a3_corresponding_prot <- left_join(a3_coords_match,
                                   coords_exon,
                                   by = join_by(gene_id, exon_5p_start),
                                   relationship = "many-to-many") |>
  separate_wider_regex(transcript_id,
                       patterns = c(seq_id = "^[A-Z0-9cel_]+\\.t?[0-9]{1,4}",
                                    prot_id = "(?:[a-z])?",
                                    tx_id = "(?:\\.[0-9]{1,2})?$"),
                       cols_remove = FALSE) |>
  select(event_id, exon_selected, prot_id) |>
  distinct() |>
  summarize(prot_ids = list(sort(prot_id)),
            .by = c(event_id, exon_selected)) |>
  summarize(iso_switch = length( intersect(prot_ids[[1]], prot_ids[[2]]) ) == 0,
            no_switch = identical(prot_ids[[1]], prot_ids[[2]]),
            n = length(prot_ids),
            .by = event_id)

stopifnot(all(a3_corresponding_prot$n == 2))


table(a3_corresponding_prot$prot_same)


a3_isoform <- a3_corresponding_prot |>
  transmute(event_id,
            isoform = case_when(iso_switch ~ "switch",
                                no_switch ~ "no switch",
                                .default = "overlap"))

table(a3_isoform$isoform)

a3_isoform |>
  count(isoform) |>
  mutate(`%` = round( 100 * n / sum(n) , digits = 1))


#~~ alternative approach ----

# # While the first approach above uses exon names, here we match the coordinates to CDS or UTR.
# # This approach appears more "correct", but turns out less conservative:
# # for example see X:15,072,953..15,073,153 where we have an overhang that's in a UTR,
# # but affects which protein isoform is expressed
# 
# 
# a3_coords_events <- coords_all |>
#   filter(event_type == "A3") |>
#   chuck("coords",1L) |>
#   select(-matches("^c[1-4]$")) |>
#   separate_wider_delim(event_id,
#                        delim = ";",
#                        names = c("gene_id", NA),
#                        cols_remove = FALSE)
# 
# # some events share the same overhang with different flanking exons
# a3_coords_overhangs <- a3_coords_events |>
#   select(gene_id, overhang_start, overhang_end) |>
#   distinct()
# 
# 
# 
# 
# 
# # Two cases:
# # 1. overhang inside a CDS exon
# # 2. overhang inside UTR exon
# 
# 
# overhangs_in_cds <- a3_coords_overhangs |>
#   inner_join(coords_cds_match,
#              by = join_by(gene_id,
#                           between(overhang_start, exon_start, exon_end),
#                           between(overhang_end, exon_start, exon_end)))
# 
# overhangs_in_utr <- a3_coords_overhangs |>
#   inner_join(coords_utr_match,
#              by = join_by(gene_id,
#                           between(overhang_start, utr_start, utr_end),
#                           between(overhang_end, utr_start, utr_end)))
# 
# overhang_with_stop_codon <- coords_stop_codon_match |>
#   inner_join(a3_coords_overhangs,
#              by = join_by(gene_id,
#                           within(codon_start, codon_end, overhang_start, overhang_end)))
# 
# overhang_with_start_codon <- coords_start_codon_match |>
#   inner_join(a3_coords_overhangs,
#              by = join_by(gene_id,
#                           within(codon_start, codon_end, overhang_start, overhang_end)))
# 
# 
# 
# 
# all_matched_overhangs <- bind_rows(
#   overhangs_in_cds,
#   overhangs_in_utr,
#   overhang_with_stop_codon,
#   overhang_with_start_codon
#   ) |>
#   select(gene_id, overhang_start, overhang_end, type)
# 
# stopifnot(nrow(anti_join(all_matched_overhangs, a3_coords_overhangs)) == 0)
# 
# 
# # a3 we are not considering (e.g. non-coding transcripts):
# a3_ignored <- anti_join(a3_coords_overhangs, all_matched_overhangs,
#                            by = join_by(gene_id, overhang_start, overhang_end))
# 
# 
# 
# # in some cases several overhangs with same coordinates but different coding potential
# # from RNA-Seq we can't tell which transcript is affected
# overhangs_with_unclear_impact <- all_matched_overhangs |>
#   count(gene_id, overhang_start, overhang_end) |> 
#   filter(n > 1) |> arrange(desc(n)) |>
#   select(gene_id, overhang_start, overhang_end)
# 
# 
# overhangs_with_impact <- all_matched_overhangs |>
#   anti_join(overhangs_with_unclear_impact,
#             by = join_by(gene_id, overhang_start, overhang_end))
# 
# stopifnot( nrow(a3_coords_overhangs) ==
#              nrow(overhangs_with_impact) + nrow(a3_ignored) + nrow(overhangs_with_unclear_impact) )
# 
# table(overhangs_with_impact$type)
# 
# 
# a3_event_types <- a3_coords_events |>
#   left_join(overhangs_with_impact,
#             by = c("gene_id", "overhang_start", "overhang_end")) |>
#   select(event_id, overhang_length, type)
# 
# table(a3_event_types$type)
# 
# 
# # Compare apporaches
# 
# comp <- left_join(
#   a3_corresponding_prot,
#   a3_event_types,
#   by = "event_id"
# )
# 
# table(comp$prot_same, comp$type, useNA = 'ifany')
# 
# left_join(
#   a3_corresponding_prot,
#   a3_event_types,
#   by = "event_id"
# ) |>
#   filter(!prot_same & type == "UTR")
# 
# 
# 
# left_join(a3_coords_match,
#           coords_exon,
#           by = join_by(gene_id, exon_5p_start),
#           relationship = "many-to-many") |>
#   separate_wider_regex(transcript_id,
#                        patterns = c(seq_id = "^[A-Z0-9cel_]+\\.t?[0-9]{1,4}",
#                                     prot_id = "(?:[a-z])?",
#                                     tx_id = "(?:\\.[0-9]{1,2})?$"),
#                        cols_remove = FALSE) |>
#   select(event_id, exon_selected, prot_id) |>
#   distinct() |>
#   summarize(prot_ids = paste(sort(prot_id), collapse = ","),
#             .by = c(event_id, exon_selected)) |>
#   filter(event_id == "WBGene00004781;A3:III:6451774-6451821:6451766-6451821:-")
# 
# 
# 
# left_join(
#   a3_corresponding_prot,
#   a3_event_types,
#   by = "event_id"
# ) |>
#   filter(startsWith(event_id, "WBGene00004781"))
# 






#~ Check impact ----
features_long <- qs::qread("intermediates/240918/241011_features.qs") |>
  filter( event_type == "A3" & feature == "overhang" ) |>
  left_join(a3_corresponding_prot |> select(-n),
            by = c("event_id"))


proportions_ndf <- features_long |>
  filter(iso_switch) |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  summarize(n_frame = sum(length_is_triple),
            n_not_frame = sum(!length_is_triple),
            .by = c(has_ds, is_detectable)) |>
  mutate(prop_in_frame = round( 100 * n_frame / (n_frame + n_not_frame) )) |>
  arrange(is_detectable, has_ds)

proportions_ndf


proportions_ndf |>
  select(-prop_in_frame) |>
  pivot_longer(-c(has_ds, is_detectable),
               names_to = "frame",
               values_to = "count") |>
  mutate(
    category = case_when(
      ! is_detectable ~ "not measured",
      has_ds ~ "DAS",
      .default = "not DAS") |>
      factor(levels = c("not measured", "not DAS", "DAS") |> rev()),
    frame = factor(frame, levels = c("n_not_frame", "n_frame"))
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
  geom_col(aes(x = category, y = count, fill = frame))



#~~ test ----

proportions_ndf |>
  filter(is_detectable) |>
  select(n_frame, n_not_frame) |>
  chisq.test()





#~ check impact vs DAS ----

# previous comparison was on length of overhang (multiple of 3) within the protein
# instead ask if "changing the protein" and "being DAS" are linked



a3_das_prot <- qs::qread("intermediates/240918/241011_features.qs") |>
  filter( event_type == "A3" ) |>
  select(event_id, has_ds, is_detectable) |>
  distinct() |>
  left_join(a3_isoform,
            by = c("event_id"))



proportions_ndf <- a3_das_prot |>
  summarize(n_no_switch = sum(isoform == "no switch"),
            n_switch = sum(isoform == "switch"),
            .by = c(has_ds, is_detectable)) |>
  mutate(prop_no_switch = round( 100 * n_no_switch / (n_switch + n_no_switch) )) |>
  arrange(is_detectable, has_ds)

proportions_ndf


proportions_ndf |>
  select(-prop_no_switch) |>
  pivot_longer(-c(has_ds, is_detectable),
               names_to = "isoform",
               values_to = "count") |>
  mutate(
    category = case_when(
      ! is_detectable ~ "not measured",
      has_ds ~ "DAS",
      .default = "not DAS") |>
      factor(levels = c("not measured", "not DAS", "DAS") |> rev()),
    isoform = factor(isoform, levels = c("n_switch", "n_no_switch"))
  ) |>
  ggplot() +
  theme_classic() +
  guides(fill = guide_legend(title = NULL )) +
  theme(legend.position = "top") +
  scale_fill_brewer(type = "qual",
                    labels = c("isoform-switching", "No isoform change")) +
  xlab(NULL) +
  ylab("Number of events") +
  coord_flip() +
  geom_col(aes(x = category, y = count, fill = isoform))



#~ test ----

proportions_ndf |>
  filter(is_detectable) |>
  select(n_no_switch, n_switch) |>
  chisq.test()



#~ transposed ----


a3_das_prot |>
  summarize(not_measured = sum(! is_detectable),
            not_DAS = sum(is_detectable & !has_ds),
            DAS = sum(has_ds),
            .by = c(isoform)) |>
  mutate(`% DAS` = round( 100 * DAS / (DAS + not_DAS) )) |>
  arrange(isoform)



a3_das_prot |>
  summarize(not_measured = sum(! is_detectable),
            not_DAS = sum(is_detectable & !has_ds),
            DAS = sum(has_ds),
            .by = c(isoform)) |>
  pivot_longer(-isoform,
               names_to = "DAS",
               values_to = "n") |>
  filter(DAS != "not_measured") |>
  ggplot() +
  theme_classic() +
  guides(fill = guide_legend(title = NULL )) +
  theme(legend.position = "top") +
  scale_fill_brewer(type = "qual") +
  xlab(NULL) +
  ylab("Number of events") +
  coord_flip() +
  geom_col(aes(x = isoform, y = n, fill = DAS))





#~ test ----

a3_das_prot |>
  summarize(not_measured = sum(! is_detectable),
            not_DAS = sum(is_detectable & !has_ds),
            DAS = sum(has_ds),
            .by = c(isoform)) |>
  filter(isoform != "overlap") |>
  select(not_DAS, DAS) |>
  chisq.test()








# A5ss ----

#~ annotate impact on protein ----


a5_coords <- coords_all |>
  filter(event_type == "A5") |>
  chuck("coords", 1L) |>
  select(-matches("^c[1-4]$")) |>
  separate_wider_delim(event_id,
                       delim = ";",
                       names = c("gene_id", NA),
                       cols_remove = FALSE)

a5_coords_match <- a5_coords |>
  mutate(shorter_exon_3p_end = if_else(strand == "+", overhang_start - 1L, overhang_end + 1L),
         longer_exon_3p_end = if_else(strand == "+", overhang_end, overhang_start)) |>
  select(event_id, gene_id,
         ends_with("exon_3p_end")) |>
  pivot_longer(cols = ends_with("exon_3p_end"),
               names_to = c("exon_selected", NA,NA,NA),
               names_sep = "_",
               values_to = "exon_3p_end")






a5_corresponding_prot <- left_join(a5_coords_match,
                                   coords_exon,
                                   by = join_by(gene_id, exon_3p_end),
                                   relationship = "many-to-many") |>
  separate_wider_regex(transcript_id,
                       patterns = c(seq_id = "^[A-Z0-9cel_]+\\.t?[0-9]{1,4}",
                                    prot_id = "(?:[a-z])?",
                                    tx_id = "(?:\\.[0-9]{1,2})?$"),
                       cols_remove = FALSE) |>
  select(event_id, exon_selected, prot_id) |>
  distinct() |>
  summarize(prot_ids = list(sort(prot_id)),
            .by = c(event_id, exon_selected)) |>
  summarize(prot_same = length( intersect(prot_ids[[1]], prot_ids[[2]]) ) > 0,
            n = length(prot_ids),
            .by = event_id)

stopifnot(all(a5_corresponding_prot$n == 2))


table(a5_corresponding_prot$prot_same)




# ~~ alternative approach ----
# a5_event_coords <- coords_all |>
#   filter(event_type == "A5") |>
#   chuck("coords",1L) |>
#   select(-matches("^c[1-4]$")) |>
#   separate_wider_delim(event_id,
#                        delim = ";",
#                        names = c("gene_id", NA),
#                        cols_remove = FALSE)
# 
# # some events share the same overhang with different flanking exons
# a5_coords_overhang <- a5_event_coords |>
#   select(gene_id, overhang_start, overhang_end) |>
#   distinct()
# 
# 
# 
# 
# 
# # Two cases:
# # 1. overhang inside a CDS exon
# # 2. overhang inside UTR exon
# 
# 
# overhangs_in_cds <- a5_coords_overhang |>
#   inner_join(coords_cds_match,
#              by = join_by(gene_id,
#                           between(overhang_start, exon_start, exon_end),
#                           between(overhang_end, exon_start, exon_end)))
# 
# overhangs_in_utr <- a5_coords_overhang |>
#   inner_join(coords_utr_match,
#              by = join_by(gene_id,
#                           between(overhang_start, utr_start, utr_end),
#                           between(overhang_end, utr_start, utr_end)))
# 
# overhang_with_stop_codon <- coords_stop_codon_match |>
#   inner_join(a5_coords_overhang,
#              by = join_by(gene_id,
#                           within(codon_start, codon_end, overhang_start, overhang_end)))
# 
# 
# overhang_with_start_codon <- coords_start_codon_match |>
#   inner_join(a5_coords_overhang,
#              by = join_by(gene_id,
#                           within(codon_start, codon_end, overhang_start, overhang_end)))
# 
# 
# 
# 
# 
# all_matched_overhangs <- bind_rows(
#   overhangs_in_cds,
#   overhangs_in_utr,
#   overhang_with_stop_codon,
#   overhang_with_start_codon
# ) |>
#   select(gene_id, overhang_start, overhang_end, type)
# 
# stopifnot(nrow(anti_join(all_matched_overhangs, a5_coords_overhang)) == 0)
# 
# 
# # a5 we are not considering (e.g. non-coding transcripts):
# a5_ignored <- anti_join(a5_coords_overhang, all_matched_overhangs,
#                         by = join_by(gene_id, overhang_start, overhang_end))
# 
# 
# 
# # in some cases several overhangs with same coordinates but different coding potential
# # from RNA-Seq we can't tell which transcript is affected
# overhangs_with_unclear_impact <- all_matched_overhangs |>
#   count(gene_id, overhang_start, overhang_end) |> 
#   filter(n > 1) |> arrange(desc(n)) |>
#   select(gene_id, overhang_start, overhang_end)
# 
# 
# overhangs_with_impact <- all_matched_overhangs |>
#   anti_join(overhangs_with_unclear_impact,
#             by = join_by(gene_id, overhang_start, overhang_end))
# 
# stopifnot( nrow(a5_coords_overhang) ==
#              nrow(overhangs_with_impact) + nrow(a5_ignored) + nrow(overhangs_with_unclear_impact) )
# 
# table(overhangs_with_impact$type)
# 
# 
# a5_event_types <- a5_event_coords |>
#   left_join(overhangs_with_impact,
#             by = c("gene_id", "overhang_start", "overhang_end")) |>
#   select(event_id, overhang_length, type)
# 
# table(a5_event_types$type)
# 
# 
# 
# # Compare apporaches
# 
# comp <- left_join(
#   a5_corresponding_prot,
#   a5_event_types,
#   by = "event_id"
# )
# 
# table(comp$prot_same, comp$type, useNA = 'ifany')
# 
# left_join(
#   a5_corresponding_prot,
#   a5_event_types,
#   by = "event_id"
# ) |>
#   filter(!prot_same & type == "UTR")
# 
# 
# 
# left_join(a5_coords_match,
#           coords_exon,
#           by = join_by(gene_id, exon_3p_end),
#           relationship = "many-to-many") |>
#   separate_wider_regex(transcript_id,
#                        patterns = c(seq_id = "^[A-Z0-9cel_]+\\.t?[0-9]{1,4}",
#                                     prot_id = "(?:[a-z])?",
#                                     tx_id = "(?:\\.[0-9]{1,2})?$"),
#                        cols_remove = FALSE) |>
#   select(event_id, exon_selected, prot_id) |>
#   distinct() |>
#   summarize(prot_ids = paste(sort(prot_id), collapse = ","),
#             .by = c(event_id, exon_selected)) |>
#   filter(event_id == "WBGene00001518;A5:III:6319606-6320319:6319456-6320319:+")
# 
# 
# 
# left_join(
#   a3_corresponding_prot,
#   a3_event_types,
#   by = "event_id"
# ) |>
#   filter(startsWith(event_id, "WBGene00004781"))







#~ Check impact ----
features_long <- qs::qread("intermediates/240918/241011_features.qs") |>
  filter( event_type == "A5" & feature == "overhang" ) |>
  left_join(a5_corresponding_prot |> select(-n),
            by = c("event_id"))


proportions_ndf <- features_long |>
  filter(!prot_same) |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  summarize(n_frame = sum(length_is_triple),
            n_not_frame = sum(!length_is_triple),
            .by = c(has_ds, is_detectable)) |>
  mutate(prop_in_frame = round( 100 * n_frame / (n_frame + n_not_frame) )) |>
  arrange(is_detectable, has_ds)


proportions_ndf


proportions_ndf |>
  select(-prop_in_frame) |>
  pivot_longer(-c(has_ds, is_detectable),
               names_to = "frame",
               values_to = "count") |>
  mutate(
    category = case_when(
      ! is_detectable ~ "not measured",
      has_ds ~ "DAS",
      .default = "not DAS") |>
      factor(levels = c("not measured", "not DAS", "DAS") |> rev()),
    frame = factor(frame, levels = c("n_not_frame", "n_frame"))
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
  geom_col(aes(x = category, y = count, fill = frame))



#~ test ----

proportions_ndf |>
  filter(is_detectable) |>
  select(n_frame, n_not_frame) |>
  chisq.test()








# Alt first ----

#~ annot impact ----

af_coords <- coords_all |>
  filter(event_type == "AF") |>
  chuck("coords", 1L) |>
  select(-matches("^c[1-6]$")) |>
  separate_wider_delim(event_id,
                       delim = ";",
                       names = c("gene_id", NA),
                       cols_remove = FALSE)

af_coords_match <- af_coords |>
  mutate(distal_3p_end = if_else(strand == "+", distal_exon_end, distal_exon_start),
         proximal_3p_end = if_else(strand == "+", proximal_exon_end, proximal_exon_start)) |>
  select(event_id, gene_id,
         distal_3p_end, proximal_3p_end) |>
  pivot_longer(cols = ends_with("3p_end"),
               names_to = c("exon_position", NA,NA),
               names_sep = "_",
               values_to = "exon_3p_end")






af_corresponding_prot <- left_join(af_coords_match,
                                   coords_exon,
                                   by = join_by(gene_id, exon_3p_end),
                                   relationship = "many-to-many") |>
  separate_wider_regex(transcript_id,
                       patterns = c(seq_id = "^[A-Z0-9cel_]+\\.t?[0-9]{1,4}",
                                    prot_id = "(?:[a-z])?",
                                    tx_id = "(?:\\.[0-9]{1,2})?$"),
                       cols_remove = FALSE) |>
  select(event_id, exon_position, prot_id) |>
  distinct() |>
  summarize(prot_ids = paste(sort(prot_id), collapse = ","),
            .by = c(event_id, exon_position)) |>
  summarize(prot_same = identical(prot_ids[[1]], prot_ids[[2]]),
            .by = event_id)

table(af_corresponding_prot$prot_same)



#~ Check impact ----
features_long <- qs::qread("intermediates/240918/241011_features.qs") |>
  filter( event_type == "AF" & feature == "overhang" ) |>
  left_join(af_corresponding_prot |> select(-n),
            by = c("event_id"))


proportions_ndf <- features_long |>
  filter(!prot_same) |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  summarize(n_frame = sum(length_is_triple),
            n_not_frame = sum(!length_is_triple),
            .by = c(has_ds, is_detectable)) |>
  mutate(prop_in_frame = round( 100 * n_frame / (n_frame + n_not_frame) )) |>
  arrange(is_detectable, has_ds)


proportions_ndf


proportions_ndf |>
  select(-prop_in_frame) |>
  pivot_longer(-c(has_ds, is_detectable),
               names_to = "frame",
               values_to = "count") |>
  mutate(
    category = case_when(
      ! is_detectable ~ "not measured",
      has_ds ~ "DAS",
      .default = "not DAS") |>
      factor(levels = c("not measured", "not DAS", "DAS") |> rev()),
    frame = factor(frame, levels = c("n_not_frame", "n_frame"))
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
  geom_col(aes(x = category, y = count, fill = frame))



#~ test ----

proportions_ndf |>
  filter(is_detectable) |>
  select(n_frame, n_not_frame) |>
  chisq.test()








# Alt last ----



al_coords <- coords_all |>
  filter(event_type == "AL") |>
  chuck("coords", 1L) |>
  select(-matches("^c[1-6]$")) |>
  separate_wider_delim(event_id,
                       delim = ";",
                       names = c("gene_id", NA),
                       cols_remove = FALSE)

al_coords_match <- al_coords |>
  mutate(distal_5p_start = if_else(strand == "+", distal_exon_start, distal_exon_end),
         proximal_5p_start = if_else(strand == "+", proximal_exon_start, proximal_exon_end)) |>
  select(event_id, gene_id,
         distal_5p_start, proximal_5p_start) |>
  pivot_longer(cols = ends_with("5p_start"),
               names_to = c("exon_position", NA,NA),
               names_sep = "_",
               values_to = "exon_5p_start")






al_corresponding_prot <- left_join(al_coords_match,
                                   coords_exon,
                                   by = join_by(gene_id, exon_5p_start),
                                   relationship = "many-to-many") |>
  separate_wider_regex(transcript_id,
                       patterns = c(seq_id = "^[A-Z0-9cel_]+\\.t?[0-9]{1,4}",
                                    prot_id = "(?:[a-z])?",
                                    tx_id = "(?:\\.[0-9]{1,2})?$"),
                       cols_remove = FALSE) |>
  select(event_id, exon_position, prot_id) |>
  distinct() |>
  summarize(prot_ids = paste(sort(prot_id), collapse = ","),
            .by = c(event_id, exon_position)) |>
  summarize(prot_same = identical(prot_ids[[1]], prot_ids[[2]]),
            .by = event_id)

table(al_corresponding_prot$prot_same)




al_corresponding_prot |>
  slice_sample(n = 10)

al_corresponding_prot |>
  filter(startsWith(event_id, "WBGene00005648"))


al_coords |>
  filter(gene_id == "WBGene00005648")











































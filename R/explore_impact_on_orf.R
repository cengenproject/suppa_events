

# Idea: for each event (focus on SE first), does it affect the ORF?


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




























# are frame-preserving events intercepting stop codons?

coords_stop_codon <- rtracklayer::readGFF(wbData::wb_get_gtf_path(289),
                                          filter = list(type = "stop_codon")) |>
  rename(seqnames = seqid) |>
  as_granges()
coords_start_codon <- rtracklayer::readGFF(wbData::wb_get_gtf_path(289),
                                          filter = list(type = "start_codon")) |>
  rename(seqnames = seqid) |>
  as_granges()




a3_with_stop <- join_overlap_intersect(a3_overhang_gr, coords_stop_codon)$event_id
a3_with_start <- join_overlap_intersect(a3_overhang_gr, coords_start_codon)$event_id


features_long |>
  filter(event_type == "A3" & feature == "overhang" ) |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  select(event_id, has_ds, is_detectable, length_is_triple) |>
  mutate(with_end = (event_id %in% a3_with_stop) | (event_id %in% a3_with_start)) |>
  count(is_detectable, has_ds, length_is_triple, with_end) |>
  pivot_wider(id_cols = 1:3,
              names_from = "with_end",
              values_from = "n",
              values_fill = 0)





a5_with_stop <- join_overlap_intersect(a5_overhang_gr, coords_stop_codon)$event_id
a5_with_start <- join_overlap_intersect(a5_overhang_gr, coords_start_codon)$event_id

features_long |>
  filter(event_type == "A5" & feature == "overhang" ) |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  select(event_id, has_ds, is_detectable, length_is_triple) |>
  mutate(with_end = (event_id %in% a5_with_stop) | (event_id %in% a5_with_start)) |>
  count(is_detectable, has_ds, length_is_triple, with_end) |>
  pivot_wider(id_cols = 1:3,
              names_from = "with_end",
              values_from = "n",
              values_fill = 0)






se_with_stop <- join_overlap_intersect(se_inclusion, coords_stop_codon)$event_id
se_with_start <- join_overlap_intersect(se_inclusion, coords_start_codon)$event_id

features_long |>
  filter(event_type == "SE" & feature == "exon" ) |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  select(event_id, has_ds, is_detectable, length_is_triple) |>
  mutate(with_end = (event_id %in% se_with_stop) | (event_id %in% se_with_start)) |>
  count(is_detectable, has_ds, length_is_triple, with_end) |>
  pivot_wider(id_cols = 1:3,
              names_from = "with_end",
              values_from = "n",
              values_fill = 0)



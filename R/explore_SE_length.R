# Inits ----

library(tidyverse)


library(wbData)

gids <- wb_load_gene_ids(289)

export_dir <- "data/outs/2408_fig"


dpsi <- qs::qread("intermediates/240425_dpsi/filt_dpsi_thres_integrated.qs")





source("R/extract_event_coordinates.R")

coords_all <- dpsi |>
  select(event_id, event_type, gene_id, gene_name, event_coordinates) |>
  distinct() |>
  nest(.by = event_type) |>
  mutate(coords = map2(event_type, data,
                       ~ extract_coords(.x, .y[["event_coordinates"]]))) |>
  mutate(coords = map2(data, coords, ~bind_cols(.x["event_id"],
                                                .y)))
                                                
se_coords <- coords_all |> filter(event_type == "SE") |> pull(coords)
se_coords <- se_coords[[1]]

d_se <- dpsi |>
  filter(event_type == "SE") |>
  filter(detectable) |>
  left_join(coords_all$coords[[which(coords_all$event_type == "SE")]],
            by = "event_id") |>
  summarize(has_ds = any(is_ds),
            .by = c(upstream_intron_length, downstream_intron_length,
                    exon_length,
                    event_id, gene_id, gene_name)) |>
  pivot_longer(-c(event_id,gene_id, gene_name, has_ds),
               names_to = "feature",
               names_pattern = "^([a-z_]+)_length$",
               values_to = "length")

d_se_med <- d_se |>
  filter(feature == "exon") |>
  summarize(median = median(length),
            .by = has_ds)

d_se |>
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



d_se |>
  filter(feature == "exon") |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  summarize(n_frame = sum(length_is_triple),
         n_not_frame = sum(!length_is_triple),
         .by = c(has_ds)) |>
  mutate(prop_in_frame = round( 100 * n_frame / (n_frame + n_not_frame) ))


d_se |>
  filter(feature == "exon") |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  filter(has_ds, !length_is_triple) |>
  slice_sample(n = 1) |>
  as.data.frame()




dpsi |>
  filter(gene_name == "par-3",
         event_type == "SE",
         detectable) |>
  filter(p.val < .1) |>
  arrange(desc(abs(dPSI))) |>
  View()




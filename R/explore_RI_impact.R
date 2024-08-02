# Inits ----

library(tidyverse)


library(wbData)

tx2g <- wb_load_tx2gene(289)
gids <- wb_load_gene_ids(289)

export_dir <- "data/outs/2408_fig"


dpsi <- qs::qread("intermediates/240425_dpsi/filt_dpsi.qs")





source("R/extract_event_coordinates.R")

coords_all <- dpsi |>
  select(event_id, event_type, gene_id, gene_name, event_coordinates) |>
  distinct() |>
  nest(.by = event_type) |>
  mutate(coords = map2(event_type, data,
                       ~ extract_coords(.x, .y[["event_coordinates"]]))) |>
  mutate(coords = map2(data, coords, ~bind_cols(.x["event_id"],
                                                .y)))

ri_coords <- coords_all |> filter(event_type == "RI") |> pull(coords)
ri_coords <- ri_coords[[1]]

d_ri <- dpsi |>
  filter(event_type == "RI") |>
  left_join(coords_all$coords[[which(coords_all$event_type == "RI")]],
            by = "event_id") |>
  summarize(has_ds = any(is_ds),
            is_detectable = any(detectable),
            .by = c(intron_length, upstream_exon_length,
                    downstream_exon_length,
                    event_id, gene_id, gene_name)) |>
  pivot_longer(-c(event_id,gene_id, gene_name, has_ds, is_detectable),
               names_to = "feature",
               names_pattern = "^([a-z_]+)_length$",
               values_to = "length")

d_ri_med <- d_ri |>
  filter(is_detectable,
         feature == "intron") |>
  summarize(median = median(length),
            .by = has_ds)

d_ri |>
  filter(is_detectable) |>
  filter(feature == "intron") |>
  mutate(has_ds = if_else(has_ds, "dAS", "non-dAS") |>
           factor(levels = c("non-dAS", "dAS"))) |>
  ggplot() +
  theme_classic() +
  ggridges::geom_density_ridges(aes(x = length, y = has_ds, fill = has_ds),
                                scale = 3) +
  geom_vline(aes(xintercept = median, color = has_ds),
             data = d_ri_med,
             linetype = 'dashed', linewidth = 1) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_color_manual(values = c("grey50", "red3")) +
  scale_y_discrete(expand = c(0.01, 0), labels = NULL) +
  # scale_y_continuous(limits = c(0,1.2)) +
  theme(legend.position = "none") +
  xlab("Retained intron length (bp)") + ylab(NULL)



d_ri |>
  filter(feature == "intron") |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  summarize(n_frame = sum(length_is_triple),
            n_not_frame = sum(!length_is_triple),
            .by = c(is_detectable, has_ds)) |>
  mutate(prop_in_frame = round( 100 * n_frame / (n_frame + n_not_frame) ))


d_ri |>
  filter(feature == "intron") |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  filter(is_detectable, has_ds, length_is_triple) |>
  slice_sample(n = 1) |>
  as.data.frame()



dpsi |>
  filter(gene_name == "pmr-1",
         event_type == "RI") |>
  View()


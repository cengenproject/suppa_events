# Inits ----

library(tidyverse)


library(wbData)

gids <- wb_load_gene_ids(289)

export_dir <- "data/outs/2408_fig"


dpsi <- qs::qread("intermediates/240425_dpsi/filt_dpsi_thres_integrated.qs")






##

source("R/extract_event_coordinates.R")

coords_all <- dpsi |>
  select(event_id, event_type, gene_id, gene_name, event_coordinates) |>
  distinct() |>
  nest(.by = event_type) |>
  mutate(coords = map2(event_type, data,
                       ~ extract_coords(.x, .y[["event_coordinates"]]))) |>
  mutate(coords = map2(data, coords, ~bind_cols(.x["event_id"],
                                                .y)))

a3_coords <- ( coords_all |> filter(event_type == "A3") |> pull(coords) )[[1]]
a5_coords <- ( coords_all |> filter(event_type == "A5") |> pull(coords) )[[1]]


#~  A3 ----
d_a3 <- dpsi |>
  filter(detectable) |>
  filter(event_type == "A3") |>
  left_join(coords_all$coords[[which(coords_all$event_type == "A3")]],
            by = "event_id") |>
  summarize(has_ds = any(is_ds),
            .by = c(intron_length, overhang_length,
                    event_id, gene_name)) |>
  pivot_longer(-c(event_id, gene_name, has_ds),
               names_to = "feature",
               names_pattern = "^([a-z]+)_length$",
               values_to = "length")


#~ A5 ----

d_a5 <- dpsi |>
  filter(detectable) |>
  filter(event_type == "A5") |>
  left_join(coords_all$coords[[which(coords_all$event_type == "A5")]],
            by = "event_id") |>
  summarize(has_ds = any(is_ds),
            .by = c(intron_length, overhang_length,
                    event_id, gene_name)) |>
  pivot_longer(-c(event_id, gene_name, has_ds),
               names_to = "feature",
               names_pattern = "^([a-z]+)_length$",
               values_to = "length")


d_a3_med <- d_a3 |>
  filter(feature == "overhang") |>
  summarize(median = median(length),
            .by = has_ds)

d_a5_med <- d_a5 |>
  filter(feature == "overhang") |>
  summarize(median = median(length),
            .by = has_ds)


d_a3 |>
  filter(feature == "overhang") |>
  mutate(has_ds = if_else(has_ds, "dAS", "non-dAS") |>
           factor(levels = c("non-dAS", "dAS"))) |>
  ggplot() +
  theme_classic() +
  ggridges::geom_density_ridges(aes(x = length, y = has_ds, fill = has_ds),
                                scale = 3) +
  geom_vline(aes(xintercept = median, color = has_ds),
             data = d_a3_med,
             linetype = 'dashed', linewidth = 1) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_color_manual(values = c("grey50", "red3")) +
  scale_y_discrete(expand = c(0.01, 0), labels = NULL) +
  # scale_y_continuous(limits = c(0,1.2)) +
  theme(legend.position = "none") +
  xlab("alt 3' Overhang length (bp)") + ylab(NULL)

d_a5 |>
  filter(feature == "overhang") |>
  mutate(has_ds = if_else(has_ds, "dAS", "non-dAS") |>
           factor(levels = c("non-dAS", "dAS"))) |>
  ggplot() +
  theme_classic() +
  ggridges::geom_density_ridges(aes(x = length, y = has_ds, fill = has_ds),
                                scale = 3) +
  geom_vline(aes(xintercept = median, color = has_ds),
             data = d_a5_med,
             linetype = 'dashed', linewidth = 1) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  scale_color_manual(values = c("grey50", "red3")) +
  scale_y_discrete(expand = c(0.01, 0), labels = NULL) +
  # scale_y_continuous(limits = c(0,1.2)) +
  theme(legend.position = "none") +
  xlab("alt 5' Overhang length (bp)") + ylab(NULL)



d_a3 |>
  filter(feature == "overhang") |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  summarize(n_frame = sum(length_is_triple),
            n_not_frame = sum(!length_is_triple),
            .by = has_ds) |>
  mutate(prop_in_frame = round( 100 * n_frame / (n_frame + n_not_frame) ))

# test
conting_a3 <- d_a3 |>
  filter(feature == "overhang") |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  summarize(n_frame = sum(length_is_triple),
            n_not_frame = sum(!length_is_triple),
            .by = has_ds) |>
  column_to_rownames("has_ds") |> as.matrix()

chisq.test(conting_a3)
prop.test(conting_a3)



d_a5 |>
  filter(feature == "overhang") |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  summarize(n_frame = sum(length_is_triple),
            n_not_frame = sum(!length_is_triple),
            .by = has_ds) |>
  mutate(prop_in_frame = round( 100 * n_frame / (n_frame + n_not_frame) ))

# test
conting_a5 <- d_a5 |>
  filter(feature == "overhang") |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  summarize(n_frame = sum(length_is_triple),
            n_not_frame = sum(!length_is_triple),
            .by = has_ds) |>
  column_to_rownames("has_ds") |> as.matrix()

chisq.test(conting_a5)


## A3
d_a3 |>
  filter(feature == "overhang") |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  filter(has_ds, !length_is_triple) |>
  slice_sample(n = 1) |>
  as.data.frame()



dpsi |>
  filter(gene_name == "M03C11.3",
         event_type == "A3",
         detectable) |>
  filter(p.val < .1) |>
  arrange(desc(abs(dPSI))) |>
  View()


## A5

d_a5 |>
  filter(feature == "overhang") |>
  mutate(length_is_triple = (length %% 3) == 0) |>
  filter(has_ds, !length_is_triple) |>
  slice_sample(n = 1) |>
  as.data.frame()




dpsi |>
  filter(gene_name == "ubr-5",
         event_type == "A5",
         detectable) |>
  filter(p.val < .1) |>
  arrange(desc(abs(dPSI))) |>
  View()



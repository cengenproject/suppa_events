


# Inits ----

library(tidyverse)


library(wbData)

tx2g <- wb_load_tx2gene(289)
gids <- wb_load_gene_ids(289)

export_dir <- "data/outs/240305_fig"


# Load ----


dpsi_dta <- read.table("data/240304_dpsi/240304.dpsi") |>
  rownames_to_column("event_id") |>
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
                       cols_remove = FALSE)

# qs::qsave(dpsi, "intermediates/240305_dpsi/dpsi.qs")


psi <- read.delim("data/240301b_psiPerEvent.psi") |>
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



tpm <- read.delim("data/231208_str_q_tx_TPM.tsv") |>
  rownames_to_column("transcript_id") |>
  as_tibble() |>
  mutate(gene_id = wb_tx2g(transcript_id, tx2g, warn_missing = TRUE),
         .after = "transcript_id")



stopifnot(all.equal(sort(unique(psi_lg$event_id)),
                    sort(unique(dpsi$event_id))))

stopifnot(all.equal(
  sort(unique(psi_lg$neuron_id)) |>
    setdiff(c("ADF", "M4", "Ref")),
  sort(unique(
    union(dpsi$neurA,
          dpsi$neurB))
  )
))



#~ check a few random examples ----
my_ev <- sample(psi_lg$event_id, 1)
my_neurA <- sample(psi_lg$neuron_id |> setdiff(c("ADF", "M4")), 1)
my_neurB <- sample(psi_lg$neuron_id |> setdiff(c("ADF", "M4")), 1)


psi_lg |>
  filter(event_id == my_ev,
         neuron_id %in% c(my_neurA, my_neurB)) |>
  summarize(mean_PSI = mean(PSI, na.rm = TRUE),
            .by = neuron_id) |>
  (\(x) print(x))() |>
  pull(mean_PSI) |> diff()

dpsi |>
  filter(event_id == my_ev,
         neurA %in% c(my_neurA, my_neurB),
         neurB %in% c(my_neurA, my_neurB))






#/ ======= Analysis ====== / ----



dpsi <- qs::qread("intermediates/240305_dpsi/dpsi.qs") |>
  filter(neurA != "Ref",
         neurB != "Ref")

neurons_here <- unique(union(dpsi$neurA,
                             dpsi$neurB))

# expression
gene_expression_table <- as.data.frame(cengenDataSC::cengen_sc_3_bulk > 0) |>
  rownames_to_column("gene_id") |>
  as_tibble() |>
  pivot_longer(-gene_id,
               names_to = "neuron",
               values_to = "gene_is_expressed")

dpsi <- dpsi |>
  left_join(gene_expression_table |> rename(expr_in_neurA = gene_is_expressed),
            by = c("gene_id", neurA = "neuron")) |>
  left_join(gene_expression_table |> rename(expr_in_neurB = gene_is_expressed),
            by = c("gene_id", neurB = "neuron")) |>
  mutate(detectable = expr_in_neurA & expr_in_neurB) |>
  select( -expr_in_neurA, -expr_in_neurB)


psi_lg <- psi_lg |>
  left_join(gene_expression_table |> rename(expressed = gene_is_expressed),
            by = c("gene_id", neuron_id = "neuron"))






# Number DS ----
# note: already corrected for multiple testing (option -gc)

dpsi$p.val |> hist(breaks = 70)

table(dpsi$p.val < .05)

dpsi |>
  slice_sample(n = 1e4) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = dPSI,
                 y = -log(p.val),
                 color = p.val < .05),
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



# type vs DS
dpsi |>
  summarize(has_ds = any(p.val < 0.05),
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

table(neurons_here %in% colnames(cengenDataSC::cengen_sc_3_bulk))





dpsi |>
  summarize(gene_expressed = any(detectable),
            has_ds = any(p.val < 0.05),
            .by = c(event_id, event_type)) |>
  mutate(category = case_when(
    !gene_expressed ~ "gene not expressed",
    !has_ds ~ "no differential splicing in neurons",
    .default = "differentially spliced",
  ) |> fct_inorder()) |>
  select(event_type, category) |>
  ggplot() +
  theme_classic() +
  geom_bar(aes(x = event_type, fill = category)) +
  scale_fill_manual(values = c("grey90", "grey30", "darkred")) +
  xlab(NULL) + ylab("Number of events") +
  theme(legend.position = "top",
        legend.title = element_blank())

# ggsave("ds_per_type.pdf", path = export_dir,
#        width = 21, height = 8, units = "cm")




#/ ===  Each type  === ----


# For all events, calling c1, c2, ... the first, second,... coordinate given
# See here for correspondence to exon borders:
# https://github.com/comprna/SUPPA#generation-of-transcript-events-and-local-alternative-splicing-events
#
# here, schematics of exon borders: `###` constitutive exon, `---` alt exon, `...` intron




#~ A3 ----



# 
#
#'  +  #########...........-----######
#             e1        s2    s3
#             c1/c3     c2    c4
#
#'  -  #########-----...........######
#             e1   e2          s3
#             c3   c1          c2/c4
#
#



d_a3 <- dpsi |>
  filter(event_type == "A3") |>
  separate_wider_regex(event_coordinates,
                       patterns = c(
                         chr = "^[IVX]{1,3}", ":",
                         c1 = "[0-9]+", "-",
                         c2 = "[0-9]+", ":",
                         c3 = "[0-9]+", "-",
                         c4 = "[0-9]+", ":",
                         strand = "[+-]$"
                       )) |>
  mutate(across(c1:c4, as.integer)) |>
  mutate(test = if_else(strand == "+",
                        all(c1 == c3),
                        all(c2 == c4)),
         .by = "strand") |>
  (\(x) {stopifnot(all(x[["test"]])); select(x, -test)})() |>
  mutate(
    intron_start = if_else(strand == "+",
                           c1 + 1,
                           c1 + 1),
    intron_end = if_else(strand == "+",
                         c2 - 1,
                         c2 - 1),
    intron_length = intron_end - intron_start + 1,
    
    overhang_start = if_else(strand == "+",
                             c2,
                             c3 + 1),
    overhang_end = if_else(strand == "+",
                           c4 - 1,
                           c1),
    overhang_length = overhang_end - overhang_start +1
    )

# d_a3 |>
#   select(event_id, chr, strand,
#          intron_start, intron_end, intron_length,
#          overhang_start, overhang_end, overhang_length) |>
#   distinct()|>
#   qs::qsave("intermediates/240305_event_coords/240305_a3_lengths.qs")





#~ A5 ----

#
#'  +  #########-----...........######
#             e1   e2          s3
#             c3   c1          c2/c4
#
#'  -  #########...........-----######
#             e1        s2    s3
#             c1/c3     c2    c4
#
#



d_a5 <- dpsi |>
  filter(event_type == "A5") |> #slice_sample(n = 5) |>
  separate_wider_regex(event_coordinates,
                       patterns = c(
                         chr = "^[IVX]{1,3}", ":",
                         c1 = "[0-9]+", "-",
                         c2 = "[0-9]+", ":",
                         c3 = "[0-9]+", "-",
                         c4 = "[0-9]+", ":",
                         strand = "[+-]$"
                       )) |>
  mutate(across(c1:c4, as.integer)) |>
  mutate(test = if_else(strand == "+",
                        all(c2 == c4),
                        all(c1 == c3)),
         .by = "strand") |>
  (\(x) {stopifnot(all(x[["test"]])); select(x, -test)})() |>
  mutate(
    intron_start = if_else(strand == "+",
                           c1 + 1,
                           c1 + 1),
    intron_end = if_else(strand == "+",
                         c2 - 1,
                         c2 - 1),
    intron_length = intron_end - intron_start + 1,
    
    overhang_start = if_else(strand == "+",
                             c3 + 1,
                             c2),
    overhang_end = if_else(strand == "+",
                           c1,
                           c4 - 1),
    overhang_length = overhang_end - overhang_start +1
  )

# d_a5 |>
#   select(event_id, chr, strand,
#          intron_start, intron_end, intron_length,
#          overhang_start, overhang_end, overhang_length) |>
#   distinct()|>
#   qs::qsave("intermediates/240305_event_coords/240305_a5_lengths.qs")







#~ AF ----

#
#'  +    -----.......-------.........######
#       s1   e1      s2   e2         s3
#       c1   c2      c4   c5        c3/c6
#
#'  -  #########...........-----..........-------
#             e1        s2    e2         s3     e3
#             c1/c4     c2    c3         c5     c6
#
#



d_af <- dpsi |>
  filter(event_type == "AF") |> #slice_sample(n = 5) |>
  separate_wider_regex(event_coordinates,
                       patterns = c(
                         chr = "^[IVX]{1,3}", ":",
                         c1 = "[0-9]+", "[:-]",
                         c2 = "[0-9]+", "[:-]",
                         c3 = "[0-9]+", ":",
                         c4 = "[0-9]+", "[:-]",
                         c5 = "[0-9]+", "[:-]",
                         c6 = "[0-9]+", ":",
                         strand = "[+-]$"
                       )) |>
  mutate(across(c1:c6, as.integer)) |>
  mutate(test = if_else(strand == "+",
                        all(c3 == c6),
                        all(c1 == c4)),
         .by = "strand") |>
  (\(x) {stopifnot(all(x[["test"]])); select(x, -test)})() |>
  mutate(
    # Distal exon
    #exon
    distal_exon_start = if_else(strand == "+",
                                   c1,
                                   c5),
    distal_exon_end = if_else(strand == "+",
                                 c2,
                                 c6),
    distal_exon_length = distal_exon_end - distal_exon_start + 1,
    
    #intron
    distal_intron_start = if_else(strand == "+",
                                     c2 + 1,
                                     c1 + 1),
    distal_intron_end = if_else(strand == "+",
                                   c3 - 1,
                                   c5 - 1),
    distal_intron_length = distal_intron_end - distal_intron_start + 1,
    
    # Proximal exon
    #exon
    proximal_exon_start = if_else(strand == "+",
                                    c4,
                                    c2),
    proximal_exon_end = if_else(strand == "+",
                                  c5,
                                  c3),
    proximal_exon_length = proximal_exon_end - proximal_exon_start + 1,
    
    #intron
    proximal_intron_start = if_else(strand == "+",
                                      c5 + 1,
                                      c1 + 1),
    proximal_intron_end = if_else(strand == "+",
                                    c3 - 1,
                                    c2 - 1),
    proximal_intron_length = proximal_intron_end - proximal_intron_start + 1
  )

# d_af |>
#   select(event_id, chr, strand,
#          proximal_exon_start, proximal_exon_end, proximal_exon_length,
#          proximal_intron_start, proximal_intron_end, proximal_intron_length,
#          distal_exon_start, distal_exon_end, distal_exon_length,
#          distal_intron_start, distal_intron_end, distal_intron_length) |>
#   distinct()|>
#   qs::qsave("intermediates/240305_event_coords/240305_af_lengths.qs")







#~ AL ----

#'  +  #########...........-----..........-------
#             e1        s2    e2         s3     e3
#             c1/c4     c2    c3         c5     c6
#
#'  -    -----.......-------.........######
#       s1   e1      s2   e2         s3
#       c1   c2      c4   c5        c3/c6
#
#
#



d_al <- dpsi |>
  filter(event_type == "AL") |> #slice_sample(n = 5) |>
  separate_wider_regex(event_coordinates,
                       patterns = c(
                         chr = "^[IVX]{1,3}", ":",
                         c1 = "[0-9]+", "[:-]",
                         c2 = "[0-9]+", "[:-]",
                         c3 = "[0-9]+", ":",
                         c4 = "[0-9]+", "[:-]",
                         c5 = "[0-9]+", "[:-]",
                         c6 = "[0-9]+", ":",
                         strand = "[+-]$"
                       )) |>
  mutate(across(c1:c6, as.integer)) |>
  mutate(test = if_else(strand == "+",
                        all(c1 == c4),
                        all(c3 == c6)),
         .by = "strand") |>
  (\(x) {stopifnot(all(x[["test"]])); select(x, -test)})() |>
  mutate(
    # Distal exon
    #exon
    distal_exon_start = if_else(strand == "+",
                                c5,
                                c1),
    distal_exon_end = if_else(strand == "+",
                              c6,
                              c2),
    distal_exon_length = distal_exon_end - distal_exon_start + 1,
    
    #intron
    distal_intron_start = if_else(strand == "+",
                                  c1 + 1,
                                  c2 + 1),
    distal_intron_end = if_else(strand == "+",
                                c5 - 1,
                                c3 - 1),
    distal_intron_length = distal_intron_end - distal_intron_start + 1,
    
    # Proximal exon
    #exon
    proximal_exon_start = if_else(strand == "+",
                                  c2,
                                  c4),
    proximal_exon_end = if_else(strand == "+",
                                c3,
                                c5),
    proximal_exon_length = proximal_exon_end - proximal_exon_start + 1,
    
    #intron
    proximal_intron_start = if_else(strand == "+",
                                    c1 + 1,
                                    c5 + 1),
    proximal_intron_end = if_else(strand == "+",
                                  c2 - 1,
                                  c3 - 1),
    proximal_intron_length = proximal_intron_end - proximal_intron_start + 1
  )

# d_al |>
#   select(event_id, chr, strand,
#          proximal_exon_start, proximal_exon_end, proximal_exon_length,
#          proximal_intron_start, proximal_intron_end, proximal_intron_length,
#          distal_exon_start, distal_exon_end, distal_exon_length,
#          distal_intron_start, distal_intron_end, distal_intron_length) |>
#   distinct() |>
#   qs::qsave("intermediates/240305_event_coords/240305_al_lengths.qs")










#~ MX ----

#'  +  #########..........-----..........-------........#####
#             e1        s2    e2         s3     e3      s4
#             c1/c5     c2    c3         c6     c7      c4/c8
#
# Minus strand identical
#
#
#


# note: one event weird coordinates, excluding


d_mx <- dpsi |>
  filter(event_type == "MX") |> 
  filter(event_id != "WBGene00010673;MX:IV:12574301-12574620:12576543-12576754:12573993-12576902:12577002-12577062:+") |>
  separate_wider_regex(event_coordinates,
                       patterns = c(
                         chr = "^[IVX]{1,3}", ":",
                         c1 = "[0-9]+", "-",
                         c2 = "[0-9]+", ":",
                         c3 = "[0-9]+", "-",
                         c4 = "[0-9]+", ":",
                         c5 = "[0-9]+", "-",
                         c6 = "[0-9]+", ":",
                         c7 = "[0-9]+", "-",
                         c8 = "[0-9]+", ":",
                         strand = "[+-]$"
                       )) |>
  mutate(across(c1:c8, as.integer)) |>
  mutate(test = all(c1 == c5 & c4 == c8)) |>
  (\(x) {stopifnot(all(x[["test"]])); select(x, -test)})() |>
  mutate(
    # First exon
    #exon
    first_exon_start = c2,
    first_exon_end = c3,
    first_exon_length = first_exon_end - first_exon_start + 1,
    
    # upstream intron
    first_up_intron_start = c1 + 1,
    first_up_intron_end = c2 - 1,
    first_up_intron_length = first_up_intron_end - first_up_intron_start + 1,
    
    # downstream intron
    first_dn_intron_start = c3 + 1,
    first_dn_intron_end = c4 - 1,
    first_dn_intron_length = first_dn_intron_end - first_dn_intron_start + 1,
    
    # Second exon
    #exon
    second_exon_start = c6,
    second_exon_end = c7,
    second_exon_length = second_exon_end - second_exon_start + 1,
    
    # second upstream intron
    second_up_intron_start = c1 + 1,
    second_up_intron_end = c6 - 1,
    second_up_intron_length = second_up_intron_end - second_up_intron_start + 1,
    
    # second downstream intron
    second_dn_intron_start = c7 + 1,
    second_dn_intron_end = c4 - 1,
    second_dn_intron_length = second_dn_intron_end - second_dn_intron_start + 1
  )

# d_mx |>
#   select(event_id, chr, strand,
#          first_exon_start, first_exon_end, first_exon_length,
#          first_up_intron_start, first_up_intron_end, first_up_intron_length,
#          first_dn_intron_start, first_dn_intron_end, first_dn_intron_length,
#          second_exon_start, second_exon_end, second_exon_length,
#          second_up_intron_start, second_up_intron_end, second_up_intron_length,
#          second_dn_intron_start, second_dn_intron_end, second_dn_intron_length) |>
#   distinct() |>
#   qs::qsave("intermediates/240305_event_coords/240305_mx_lengths.qs")





#~ RI ----

#'  +  #########------------------############
#     s1       e1                s2          e2
#     c1       c2                c3          c4
#
# Minus strand identical
#
#
#




d_ri <- dpsi |>
  filter(event_type == "RI") |> 
  separate_wider_regex(event_coordinates,
                       patterns = c(
                         chr = "^[IVX]{1,3}", ":",
                         c1 = "[0-9]+", ":",
                         c2 = "[0-9]+", "-",
                         c3 = "[0-9]+", ":",
                         c4 = "[0-9]+", ":",
                         strand = "[+-]$"
                       )) |>
  mutate(across(c1:c4, as.integer)) |>
  mutate(
    # Upstream exon
    upstream_exon_start = if_else(strand == "+", c1, c3),
    upstream_exon_end = if_else(strand == "+", c2, c4),
    upstream_exon_length = upstream_exon_end - upstream_exon_start + 1,
    
    # Downstream exon
    downstream_exon_start = if_else(strand == "+", c3, c1),
    downstream_exon_end = if_else(strand == "+", c4, c2),
    downstream_exon_length = downstream_exon_end - downstream_exon_start + 1,
    
    # Intron
    intron_start = c2 + 1,
    intron_end = c3 - 1,
    intron_length = intron_end - intron_start + 1
  )

# d_ri |>
#   select(event_id, chr, strand,
#          upstream_exon_start, upstream_exon_end, upstream_exon_length,
#          downstream_exon_start, downstream_exon_end, downstream_exon_length,
#          intron_start, intron_end, intron_length) |>
#   distinct() |>
#   qs::qsave("intermediates/240305_event_coords/240305_ri_lengths.qs")







#~ SE ----

#'  +  ########...........-----------.........############
#            e1          s2          e2       s3
#            c1          c2          c3       c4
#
# Minus strand identical
#
#
#




d_se <- dpsi |>
  filter(event_type == "SE") |> 
  separate_wider_regex(event_coordinates,
                       patterns = c(
                         chr = "^[IVX]{1,3}", ":",
                         c1 = "[0-9]+", "-",
                         c2 = "[0-9]+", ":",
                         c3 = "[0-9]+", "-",
                         c4 = "[0-9]+", ":",
                         strand = "[+-]$"
                       )) |>
  mutate(across(c1:c4, as.integer)) |>
  mutate(
    # Upstream intron
    upstream_intron_start = c1 + 1,
    upstream_intron_end = c2 - 1,
    upstream_intron_length = upstream_intron_end - upstream_intron_start + 1,
    
    # Downstream intron
    downstream_intron_start = c3 + 1,
    downstream_intron_end = c4 - 1,
    downstream_intron_length = downstream_intron_end - downstream_intron_start + 1,
    
    # Exon
    exon_start = c2,
    exon_end = c3,
    exon_length = exon_end - exon_start + 1
  )

# d_se |>
#   select(event_id, chr, strand,
#          upstream_intron_start, upstream_intron_end, upstream_intron_length,
#          downstream_intron_start, downstream_intron_end, downstream_intron_length,
#          exon_start, exon_end, exon_length) |>
#   distinct() |>
#   qs::qsave("intermediates/240305_event_coords/240305_se_lengths.qs")


















# Plot lengths (see more plots below) ----

d_a3 |>
  filter(detectable) |>
  summarize(has_ds = any(p.val < 0.05),
            .by = c(intron_length, overhang_length, event_id)) |>
  pivot_longer(-c(event_id, has_ds),
               names_to = "measure",
               values_to = "length") |>
  mutate(measure = case_when(
    measure == "intron_length" ~ "Intron",
    measure == "overhang_length" ~ "Overhang"
  )) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = length, color = has_ds)) +
  facet_wrap(~measure, scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Length (bp)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))


# ggsave("lengths_a3.pdf", path = export_dir,
#        width = 9, height = 6, units = "cm")




d_a5 |>
  filter(detectable) |>
  summarize(has_ds = any(p.val < 0.05),
            .by = c(intron_length, overhang_length, event_id)) |>
  pivot_longer(-c(event_id, has_ds),
               names_to = "measure",
               values_to = "length") |>
  mutate(measure = case_when(
    measure == "intron_length" ~ "Intron",
    measure == "overhang_length" ~ "Overhang"
  )) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = length, color = has_ds)) +
  facet_wrap(~measure, scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Length (bp)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))


# ggsave("lengths_a5.pdf", path = export_dir,
#        width = 9, height = 6, units = "cm")




d_af |>
  filter(detectable) |>
  summarize(has_ds = any(p.val < 0.05),
            .by = c(distal_exon_length, distal_intron_length,
                    proximal_exon_length, proximal_intron_length,
                    event_id)) |>
  pivot_longer(-c(event_id, has_ds),
               names_to = "measure",
               values_to = "length") |>
  mutate(measure = recode(measure,
    distal_exon_length = "Distal Exon",
    distal_intron_length = "Distal Intron",
    proximal_exon_length = "Proximal Exon",
    proximal_intron_length = "Proximal Intron"
  )) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = length, color = has_ds)) +
  facet_wrap(~measure, scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Length (bp)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))


# ggsave("lengths_af.pdf", path = export_dir,
#        width = 18, height = 12, units = "cm")





d_al |>
  filter(detectable) |>
  summarize(has_ds = any(p.val < 0.05),
            .by = c(distal_exon_length, distal_intron_length,
                    proximal_exon_length, proximal_intron_length,
                    event_id)) |>
  pivot_longer(-c(event_id, has_ds),
               names_to = "measure",
               values_to = "length") |>
  mutate(measure = recode(measure,
                          distal_exon_length = "Distal Exon",
                          distal_intron_length = "Distal Intron",
                          proximal_exon_length = "Proximal Exon",
                          proximal_intron_length = "Proximal Intron"
  )) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = length, color = has_ds)) +
  facet_wrap(~measure, scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Length (bp)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))


# ggsave("lengths_al.pdf", path = export_dir,
#        width = 18, height = 12, units = "cm")




d_mx |>
  filter(detectable) |>
  summarize(has_ds = any(p.val < 0.05),
            .by = c(first_exon_length, first_up_intron_length,
                    first_dn_intron_length, second_exon_length,
                    second_up_intron_length, second_dn_intron_length,
                    event_id)) |>
  pivot_longer(-c(event_id, has_ds),
               names_to = "measure",
               values_to = "length") |>
  mutate(measure = recode(measure,
                          first_exon_length = "First Exon",
                          first_up_intron_length = "First upstream Intron",
                          first_dn_intron_length = "First downstream Intron",
                          second_exon_length = "Second Exon",
                          second_up_intron_length = "Second Upstream Intron",
                          second_dn_intron_length = "Second Downstream Intron"
  )) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = length, color = has_ds)) +
  facet_wrap(~measure, scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Length (bp)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))


# ggsave("lengths_mx.pdf", path = export_dir,
#        width = 21, height = 12, units = "cm")






d_ri |>
  filter(detectable) |>
  summarize(has_ds = any(p.val < 0.05),
            .by = c(intron_length, upstream_exon_length,
                    downstream_exon_length,
                    event_id)) |>
  pivot_longer(-c(event_id, has_ds),
               names_to = "measure",
               values_to = "length") |>
  mutate(measure = recode(measure,
                          intron_length = "Intron",
                          upstream_exon_length = "Upstream Exon",
                          downstream_exon_length = "Downstream Exon"
  )) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = length, color = has_ds)) +
  facet_wrap(~measure, scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Length (bp)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))


# ggsave("lengths_ri.pdf", path = export_dir,
#        width = 21, height = 6, units = "cm")





d_se |>
  filter(detectable) |>
  summarize(has_ds = any(p.val < 0.05),
            .by = c(upstream_intron_length, downstream_intron_length,
                    exon_length,
                    event_id)) |>
  pivot_longer(-c(event_id, has_ds),
               names_to = "measure",
               values_to = "length") |>
  mutate(measure = recode(measure,
                          exon_length = "Exon",
                          upstream_intron_length = "Upstream Intron",
                          downstream_intron_length = "Downstream Intron"
  )) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x= has_ds, y = length, color = has_ds)) +
  facet_wrap(~measure, scales = "free_y") +
  # scale_y_log10() +
  xlab(NULL) + ylab("Length (bp)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey30", "darkred"))


# ggsave("lengths_se.pdf", path = export_dir,
#        width = 21, height = 6, units = "cm")



# GC content ----

# Import after running "sequence_properties.R"

seq_dir <- "intermediates/240305_event_coords/"
seq_files <- list.files(seq_dir, pattern = "*_seq.qs",
                        full.names = FALSE)
seq_types <- str_match(seq_files,
                       "^240305_([as35eflmxri]{2})_seq\\.qs$")[,2] |>
  toupper()


seq_gc <- seq_files |>
  set_names(seq_types) |>
  imap(\(.f, .t){
    qs::qread(file.path(seq_dir, .f)) |>
        add_column(type = .t) |>
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
d_a3 <- d_a3 |>
  summarize(has_ds = any(p.val < 0.05),
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

d_a5 <- d_a5 |>
  summarize(has_ds = any(p.val < 0.05),
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



d_af <- d_af |>
  summarize(has_ds = any(p.val < 0.05 & detectable),
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


d_al <- d_al |>
  summarize(has_ds = any(p.val < 0.05 & detectable),
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


d_mx <- d_mx |>
  summarize(has_ds = any(p.val < 0.05 & detectable),
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


d_ri <- d_ri |>
  summarize(has_ds = any(p.val < 0.05 & detectable),
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


d_se <- d_se |>
  summarize(has_ds = any(p.val < 0.05 & detectable),
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
#           "intermediates/240305_events_full.qs")




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





# save ----


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

ggsave("all_together.pdf", plot = gr, path = export_dir,
       width = 40, height = 60, units = "cm")


egg::gtable_frame(ggplotGrob(gg_l_a3)) |>
  grid::grid.draw()

egg::gtable_frame(ggplotGrob(gg_l_a3))
egg::gtable_frame(ggplotGrob(gg_gc_a3))
egg::gtable_frame(ggplotGrob(gg_cs_a3))

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

hist(tests$p_val)
hist(tests$padj, breaks = 70)

tests |>
  filter(padj < 0.1) |>
  mutate(sig = cut(padj,
                   breaks =c(.1,  .05,  .01, .001,  0),
                   labels =   c("#", "*", "**", "***")))


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




# Specificity ----



psi_by_neuron <- psi_lg |>
  filter(expressed) |>
  summarize(PSI_neuron = mean(PSI, na.rm = TRUE),
            .by = c("event_id", "gene_id","event_type",
                    "neuron_id"))

hist(psi_by_neuron$PSI_neuron, breaks = 100)

ggplot(psi_by_neuron) +
  theme_classic() +
  geom_density(aes(x = PSI_neuron, fill = event_type), alpha = .3)




centered_gini <- function(x) DescTools::Gini(abs(x - .5))
gcad <- function(x){
  x <- x[! is.na(x)]
  n <- length(x)
  abs_diffs <- DescTools::CombPairs(x, x) |>
    as.matrix() |>
    matrixStats::rowDiffs() |> 
    abs()
  
  if(max(abs_diffs) == 0) return(0)
  
  DescTools::Gini(abs_diffs)*n/(n^2)
}
gcsd <- function(x){
  x <- x[! is.na(x)]
  n <- length(x)
  diffs <- DescTools::CombPairs(x, x) |>
    as.matrix() |>
    matrixStats::rowDiffs()
  
  if(all(diffs == 0)) return(0)
  
  DescTools::Gini(diffs^2)*n/(n^2)
}


psi_var <- psi_by_neuron |>
  filter(!is.na(PSI_neuron)) |>
  summarize(var = var(PSI_neuron, na.rm = TRUE),
            # gmd = GiniDistance::gmd(PSI_neuron),
            # mad = mad(PSI_neuron),
            gcsd = gcsd(PSI_neuron),
            n = n(),
         .by = c("event_id", "gene_id","event_type")) |>
  filter(n > 20)



gg_var <- psi_var |>
  # filter(n > 30) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x = event_type, y = var),
                               alpha = .5) +
  ylab("Variance")


gg_gcsd <- psi_var |>
  # filter(n > 30) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x = event_type, y = gcsd),
                               alpha = .5) +
  ylab("Specificity")

patchwork::wrap_plots(gg_var, gg_gcsd, ncol = 1)

# ggsave("variance_specificity.pdf", path = export_dir,
#        width = 21, height = 12, units = "cm")






psi_var |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = var, y = gcsd))
psi_var |> filter(var>.1 & gmd < .32) |> pull(event_id) -> ev

psi_by_neuron |> filter(event_id %in% ev) |> 
  group_by(event_id) |>
  arrange(event_id, desc(PSI_neuron)) |>
  mutate(rank = rank(1-PSI_neuron)) |>
  ggplot() +
  theme_classic() +
  geom_line(aes(x = rank, y = PSI_neuron, color = event_id)) +
  theme(legend.position = "none")


psi_var |>
  ggplot() +
  theme_classic() +
  geom_boxplot(aes(x = event_type, y = gcad))


psi_var |>
  ggplot() +
  theme_classic() +
  geom_density(aes(x = gcad, fill = event_type), alpha = .3)



# look at least and most specific

extremal_var_events <- psi_var |>
  group_by(event_type) |>
  slice_max(order_by = gcsd, n = 5) |>
  pull(event_id)


mat_extr_var <- psi_by_neuron |>
  filter(event_id %in% extremal_var_events) |>
  pivot_wider(id_cols = event_id,
              names_from = neuron_id,
              values_from = PSI_neuron) |>
  column_to_rownames("event_id") |>
  as.matrix()

annot_df <- psi_by_neuron |>
  filter(event_id %in% extremal_var_events) |>
  select(-PSI_neuron, -neuron_id) |>
  distinct() |>
  left_join(psi_var,
            by = c("event_id", "gene_id", "event_type")) |>
  select(-gene_id) |>
  column_to_rownames("event_id") |>
  arrange(event_type, var)

pheatmap::pheatmap(mat_extr_var[rownames(annot_df),],
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   show_rownames = FALSE,
                   annotation_row = annot_df)


psi_by_neuron |>
  filter(event_id %in% extremal_tau_events) |>
  # count(event_id)
  arrange(event_id) |> View()


psi_by_neuron |>
  filter(event_id %in% extremal_var_events) |>
  group_by(event_id) |>
  arrange(event_id, desc(PSI_neuron)) |>
  mutate(rank = rank(1-PSI_neuron)) |>
  ggplot() +
  theme_classic() +
  geom_line(aes(x = rank, y = PSI_neuron, color = event_id)) +
  theme(legend.position = "none")
  



psi_var |>
  group_by(event_type) |>
  slice_max(order_by = gcsd, n = 5) |>
  pull(event_id)






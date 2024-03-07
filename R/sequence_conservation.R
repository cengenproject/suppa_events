
library(GenomicFeatures)
library(genomation)

bw <- rtracklayer::import("data/UCSC_fastcons/ce11.phastCons26way.bw")

seqlevels(bw) <- c("I", "II","III","IV", "M", "V","X")



get_score <- function(gr){
  problematic <- which(width(gr) == 1)
  width(gr)[problematic] <- 10
  
  score <- ScoreMatrixBin(target = bw,
                          windows = gr,
                          bin.num = 1, weight.col = "score")
  
  res <- as.numeric(score@.Data)
  res[problematic] <- NA
  res
}


#~ A3 ----

d_a3 <- qs::qread("intermediates/240305_event_coords/240305_a3.qs")

a3_intron_gr <- GRanges(seqnames = d_a3$chr,
                        ranges = IRanges(start = d_a3$intron_start,
                                         end = d_a3$intron_end),
                        strand = d_a3$strand)


a3_overhang_gr <- GRanges(seqnames = d_a3$chr,
                          ranges = IRanges(start = d_a3$overhang_start,
                                           end = d_a3$overhang_end),
                          strand = d_a3$strand)




a3_intron_score <- get_score(a3_intron_gr)

a3_overhang_score <- get_score(a3_intron_gr)



data.frame(event_id = d_a3$event_id,
           intron_conservation = as.numeric(a3_intron_score@.Data),
           overhang_conservation = as.numeric(a3_overhang_score@.Data)) |>
  qs::qsave("intermediates/240305_event_coords/240305_a3_cons.qs")





# A5 ----
d_a5 <- qs::qread("intermediates/240305_event_coords/240305_a5.qs")

a5_intron_gr <- GRanges(seqnames = d_a5$chr,
                        ranges = IRanges(start = d_a5$intron_start,
                                         end = d_a5$intron_end),
                        strand = d_a5$strand)


a5_overhang_gr <- GRanges(seqnames = d_a5$chr,
                          ranges = IRanges(start = d_a5$overhang_start,
                                           end = d_a5$overhang_end),
                          strand = d_a5$strand)



a5_intron_score <- get_score(a5_intron_gr)

a5_overhang_score <- get_score(a5_overhang_gr)




data.frame(event_id = d_a5$event_id,
           intron_conservation = as.numeric(a5_intron_score@.Data),
           overhang_conservation = as.numeric(a5_overhang_score@.Data)) |>
  qs::qsave("intermediates/240305_event_coords/240305_a5_cons.qs")



# AF ----
d_af <- qs::qread("intermediates/240305_event_coords/240305_af.qs")

af_prox_exon_gr <- GRanges(seqnames = d_af$chr,
                           ranges = IRanges(start = d_af$proximal_exon_start,
                                            end = d_af$proximal_exon_end),
                           strand = d_af$strand)

af_prox_intron_gr <- GRanges(seqnames = d_af$chr,
                             ranges = IRanges(start = d_af$proximal_intron_start,
                                              end = d_af$proximal_intron_end),
                             strand = d_af$strand)


af_dist_exon_gr <- GRanges(seqnames = d_af$chr,
                           ranges = IRanges(start = d_af$distal_exon_start,
                                            end = d_af$distal_exon_end),
                           strand = d_af$strand)

af_dist_intron_gr <- GRanges(seqnames = d_af$chr,
                             ranges = IRanges(start = d_af$distal_intron_start,
                                              end = d_af$distal_intron_end),
                             strand = d_af$strand)




af_prox_exon_score <- get_score(af_prox_exon_gr)

af_prox_intron_score <- get_score(af_prox_intron_gr)

af_dist_exon_score <- get_score(af_dist_exon_gr)

af_dist_intron_score <- get_score(af_dist_intron_gr)




data.frame(event_id = d_af$event_id,
           proximal_exon_conservation = af_prox_exon_score,
           proximal_intron_conservation = af_prox_intron_score,
           distal_exon_conservation = af_dist_exon_score,
           distal_intron_conservation = af_dist_intron_score
) |>
  qs::qsave("intermediates/240305_event_coords/240305_af_cons.qs")





# AL ----
d_al <- qs::qread("intermediates/240305_event_coords/240305_al.qs")

al_prox_ex_gr <- GRanges(seqnames = d_al$chr,
                         ranges = IRanges(start = d_al$proximal_exon_start,
                                          end = d_al$proximal_exon_end),
                         strand = d_al$strand)

al_prox_in_gr <- GRanges(seqnames = d_al$chr,
                         ranges = IRanges(start = d_al$proximal_intron_start,
                                          end = d_al$proximal_intron_end),
                         strand = d_al$strand)

al_dist_ex_gr <- GRanges(seqnames = d_al$chr,
                         ranges = IRanges(start = d_al$distal_exon_start,
                                          end = d_al$distal_exon_end),
                         strand = d_al$strand)

al_dist_in_gr <- GRanges(seqnames = d_al$chr,
                         ranges = IRanges(start = d_al$distal_intron_start,
                                          end = d_al$distal_intron_end),
                         strand = d_al$strand)



al_p_ex_seq <- get_score(al_prox_ex_gr)
al_p_i_seq <- get_score(al_prox_in_gr)
al_d_e_seq <- get_score(al_dist_ex_gr)
al_d_i_seq <- get_score(al_dist_in_gr)




data.frame(event_id = d_al$event_id,
           proximal_exon_conservation = al_p_ex_seq,

           proximal_intron_conservation = al_p_i_seq,

           distal_exon_conservation = al_d_e_seq,

           distal_intron_conservation = al_d_i_seq
           ) |>
  qs::qsave("intermediates/240305_event_coords/240305_al_cons.qs")




# MX ----
d_mx <- qs::qread("intermediates/240305_event_coords/240305_mx.qs")

mx_f_e <- GRanges(seqnames = d_mx$chr,
                  ranges = IRanges(start = d_mx$first_exon_start,
                                   end = d_mx$first_exon_end),
                  strand = d_mx$strand)

mx_f_ui <- GRanges(seqnames = d_mx$chr,
                   ranges = IRanges(start = d_mx$first_up_intron_start,
                                    end = d_mx$first_up_intron_end),
                   strand = d_mx$strand)
mx_f_di <- GRanges(seqnames = d_mx$chr,
                   ranges = IRanges(start = d_mx$first_dn_intron_start,
                                    end = d_mx$first_dn_intron_end),
                   strand = d_mx$strand)

mx_s_e <- GRanges(seqnames = d_mx$chr,
                  ranges = IRanges(start = d_mx$second_exon_start,
                                   end = d_mx$second_exon_end),
                  strand = d_mx$strand)
mx_s_ui <- GRanges(seqnames = d_mx$chr,
                   ranges = IRanges(start = d_mx$second_up_intron_start,
                                    end = d_mx$second_up_intron_end),
                   strand = d_mx$strand)
mx_s_di <- GRanges(seqnames = d_mx$chr,
                   ranges = IRanges(start = d_mx$second_dn_intron_start,
                                    end = d_mx$second_dn_intron_end),
                   strand = d_mx$strand)

mx_f_e_s <- get_score(mx_f_e)
mx_f_ui_s <- get_score(mx_f_ui)
mx_f_di_s <- get_score(mx_f_di)
mx_s_e_s <- get_score(mx_s_e)
mx_s_ui_s <- get_score(mx_s_ui)
mx_s_di_s <- get_score(mx_s_di)



data.frame(event_id = d_mx$event_id,
           first_exon_conservation = mx_f_e_s,
           first_upstream_intron_conservation = mx_f_ui_s,
           first_downstream_intron_conservation = mx_f_di_s,

           second_exon_conservation = mx_s_e_s,
           second_upstream_intron_conservation = mx_s_ui_s,
           second_downstream_intron_conservation = mx_s_di_s
           ) |>
  qs::qsave("intermediates/240305_event_coords/240305_mx_cons.qs")




# RI ----
d_ri <- qs::qread("intermediates/240305_event_coords/240305_ri.qs")


ri_ue_gr <- GRanges(seqnames = d_ri$chr,
                    ranges = IRanges(start = d_ri$upstream_exon_start,
                                     end = d_ri$upstream_exon_end),
                    strand = d_ri$strand)

ri_de_gr <- GRanges(seqnames = d_ri$chr,
                    ranges = IRanges(start = d_ri$downstream_exon_start,
                                     end = d_ri$downstream_exon_end),
                    strand = d_ri$strand)

ri_intron_gr <- GRanges(seqnames = d_ri$chr,
                        ranges = IRanges(start = d_ri$intron_start,
                                         end = d_ri$intron_end),
                        strand = d_ri$strand)


ri_de_s <- get_score(ri_ue_gr)
ri_ue_s <- get_score(ri_de_gr)
ri_i_s <- get_score(ri_intron_gr)



data.frame(event_id = d_ri$event_id,
           downstream_exon_conservation = ri_de_s,

           upstream_exon_conservation = ri_ue_s,

           intron_conservation = ri_i_s
             ) |>
  qs::qsave("intermediates/240305_event_coords/240305_ri_cons.qs")




# SE ----
d_se <- qs::qread("intermediates/240305_event_coords/240305_se.qs")

se_ui <- GRanges(seqnames = d_se$chr,
                 ranges = IRanges(start = d_se$upstream_intron_start,
                                  end = d_se$upstream_intron_end),
                 strand = d_se$strand)

se_di <- GRanges(seqnames = d_se$chr,
                 ranges = IRanges(start = d_se$downstream_intron_start,
                                  end = d_se$downstream_intron_end),
                 strand = d_se$strand)

se_e <- GRanges(seqnames = d_se$chr,
                ranges = IRanges(start = d_se$exon_start,
                                 end = d_se$exon_end),
                strand = d_se$strand)


se_ui_seq <- get_score(se_ui)
se_di_seq <- get_score(se_di)
se_e_seq <- get_score(se_e)


data.frame(event_id = d_se$event_id,
           upstream_intron_conservation = se_ui_seq,

           downstream_intron_conservation = se_di_seq,

           exon_conservation = se_e_seq
             ) |>
  qs::qsave("intermediates/240305_event_coords/240305_se_cons.qs")











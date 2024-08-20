
# Inits ----

library(wbData)
library(Biostrings)
library(GenomicFeatures)


gseq <- readDNAStringSet(wb_get_genome_path(289))

d_all <- qs::qread("intermediates/240814_dpsi/240814_all_coords.qs")


# functions ----

get_seq_a3 <- function(d_a3){
  
  
  a3_intron_gr <- GRanges(seqnames = d_a3$chr,
                          ranges = IRanges(start = d_a3$intron_start,
                                           end = d_a3$intron_end),
                          strand = d_a3$strand)
  
  a3_intron_seq <- gseq[a3_intron_gr]
  
  
  a3_overhang_gr <- GRanges(seqnames = d_a3$chr,
                            ranges = IRanges(start = d_a3$overhang_start,
                                             end = d_a3$overhang_end),
                            strand = d_a3$strand)
  
  a3_overhang_seq <- gseq[a3_overhang_gr]
  
  
  
  data.frame(event_id = d_a3$event_id,
             intron_width = width(a3_intron_seq),
             intron_gc = letterFrequency(a3_intron_seq,
                                         letters = "GC") |>
               as.numeric(),
             overhang_width = width(a3_overhang_seq),
             overhang_gc = letterFrequency(a3_overhang_seq,
                                           letters = "GC") |>
               as.numeric())
}

get_seq_a5 <- function(d_a5){
  
  
  a5_intron_gr <- GRanges(seqnames = d_a5$chr,
                          ranges = IRanges(start = d_a5$intron_start,
                                           end = d_a5$intron_end),
                          strand = d_a5$strand)
  
  a5_intron_seq <- gseq[a5_intron_gr]
  
  
  a5_overhang_gr <- GRanges(seqnames = d_a5$chr,
                            ranges = IRanges(start = d_a5$overhang_start,
                                             end = d_a5$overhang_end),
                            strand = d_a5$strand)
  
  a5_overhang_seq <- gseq[a5_overhang_gr]
  
  
  
  data.frame(event_id = d_a5$event_id,
             intron_width = width(a5_intron_seq),
             intron_gc = letterFrequency(a5_intron_seq,
                                         letters = "GC") |>
               as.numeric(),
             overhang_width = width(a5_overhang_seq),
             overhang_gc = letterFrequency(a5_overhang_seq,
                                           letters = "GC") |>
               as.numeric())
}



get_seq_af <- function(d_af){
  
  
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
  
  
  
  
  
  af_prox_exon_seq <- gseq[af_prox_exon_gr]
  af_prox_intron_seq <- gseq[af_prox_intron_gr]
  af_dist_exon_seq <- gseq[af_dist_exon_gr]
  af_dist_intron_seq <- gseq[af_dist_intron_gr]
  
  
  
  
  
  data.frame(event_id = d_af$event_id,
             proximal_exon_width = width(af_prox_exon_seq),
             proximal_exon_gc = letterFrequency(af_prox_exon_seq,
                                         letters = "GC") |>
               as.numeric(),

             proximal_intron_width = width(af_prox_intron_seq),
             proximal_intron_gc = letterFrequency(af_prox_intron_seq,
                                            letters = "GC") |>
               as.numeric(),

             distal_exon_width = width(af_dist_exon_seq),
             distal_exon_gc = letterFrequency(af_dist_exon_seq,
                                            letters = "GC") |>
               as.numeric(),

             distal_intron_width = width(af_dist_intron_seq),
             distal_intron_gc = letterFrequency(af_dist_intron_seq,
                                            letters = "GC") |>
               as.numeric()
             )
}
get_seq_al <- function(d_al){
  
  
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
  
  
  
  al_p_ex_seq <- gseq[al_prox_ex_gr]
  al_p_i_seq <- gseq[al_prox_in_gr]
  al_d_e_seq <- gseq[al_dist_ex_gr]
  al_d_i_seq <- gseq[al_dist_in_gr]
  
  
  
  
  data.frame(event_id = d_al$event_id,
             proximal_exon_width = width(al_p_ex_seq),
             proximal_exon_gc = letterFrequency(al_p_ex_seq,
                                         letters = "GC") |>
               as.numeric(),

             proximal_intron_width = width(al_p_i_seq),
             proximal_intron_gc = letterFrequency(al_p_i_seq,
                                            letters = "GC") |>
               as.numeric(),

             distal_exon_width = width(al_d_e_seq),
             distal_exon_gc = letterFrequency(al_d_e_seq,
                                            letters = "GC") |>
               as.numeric(),

             distal_intron_width = width(al_d_i_seq),
             distal_intron_gc = letterFrequency(al_d_i_seq,
                                            letters = "GC") |>
               as.numeric()
             )
}

get_seq_mx <- function(d_mx){
  
  
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
  
  mx_f_e_s <- gseq[mx_f_e]
  mx_f_ui_s <- gseq[mx_f_ui]
  mx_f_di_s <- gseq[mx_f_di]
  mx_s_e_s <- gseq[mx_s_e]
  mx_s_ui_s <- gseq[mx_s_ui]
  mx_s_di_s <- gseq[mx_s_di]
  
  
  
  data.frame(event_id = d_mx$event_id,
             first_exon_width = width(mx_f_e_s),
             first_exon_gc = letterFrequency(mx_f_e_s,
                                             letters = "GC") |>
               as.numeric(),
             first_upstream_intron_width = width(mx_f_ui_s),
             first_upstream_intron_gc = letterFrequency(mx_f_ui_s,
                                                        letters = "GC") |>
               as.numeric(),
             first_downstream_intron_width = width(mx_f_di_s),
             first_downstream_intron_gc = letterFrequency(mx_f_di_s,
                                                          letters = "GC") |>
               as.numeric(),
             
             second_exon_width = width(mx_s_e_s),
             second_exon_gc = letterFrequency(mx_s_e_s,
                                              letters = "GC") |>
               as.numeric(),
             second_upstream_intron_width = width(mx_s_ui_s),
             second_upstream_intron_gc = letterFrequency(mx_s_ui_s,
                                                         letters = "GC") |>
               as.numeric(),
             second_downstream_intron_width = width(mx_s_di_s),
             second_downstream_intron_gc = letterFrequency(mx_s_di_s,
                                                           letters = "GC") |>
               as.numeric()
  )
}


get_seq_ri <- function(d_ri){
  
  
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
  
  
  ri_de_s <- gseq[ri_ue_gr]
  ri_ue_s <- gseq[ri_de_gr]
  ri_i_s <- gseq[ri_intron_gr]
  
  
  
  data.frame(event_id = d_ri$event_id,
             downstream_exon_width = width(ri_de_s),
             downstream_exon_gc = letterFrequency(ri_de_s,
                                                  letters = "GC") |>
               as.numeric(),
             
             upstream_exon_width = width(ri_ue_s),
             upstream_exon_gc = letterFrequency(ri_ue_s,
                                                letters = "GC") |>
               as.numeric(),
             
             intron_width = width(ri_i_s),
             intron_gc = letterFrequency(ri_i_s,
                                         letters = "GC") |>
               as.numeric())
}


get_seq_se <- function(d_se){
  
  
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
  
  
  se_ui_seq <- gseq[se_ui]
  se_di_seq <- gseq[se_di]
  se_e_seq <- gseq[se_e]
  
  
  data.frame(event_id = d_se$event_id,
             upstream_intron_width = width(se_ui_seq),
             upstream_intron_gc = letterFrequency(se_ui_seq,
                                                  letters = "GC") |>
               as.numeric(),
             
             downstream_intron_width = width(se_di_seq),
             downstream_intron_gc = letterFrequency(se_di_seq,
                                                    letters = "GC") |>
               as.numeric(),
             
             exon_width = width(se_e_seq),
             exon_gc = letterFrequency(se_e_seq,
                                       letters = "GC") |>
               as.numeric())
  
}

get_seq_properties <- function(event_type, d){
    switch (event_type,
            A3 = get_seq_a3(d),
            A5 = get_seq_a5(d),
            AL = get_seq_al(d),
            AF = get_seq_af(d),
            MX = get_seq_mx(d),
            RI = get_seq_ri(d),
            SE = get_seq_se(d)
    )
  
}




# Main ----

all_seq_properties <- list()

for(ev_type in c("A5","A3","AF","AL","MX","RI","SE")){
  d_ev_type <- d_all[d_all$event_type == ev_type, "coords"][[1]][[1]] |>
    dplyr::select(-matches("^c[1-9]$"))
  
  all_seq_properties[[ev_type]] <- get_seq_properties(ev_type, d_ev_type)
}

qs::qsave(all_seq_properties, "intermediates/240814_dpsi/240814_seq_properties.qs")











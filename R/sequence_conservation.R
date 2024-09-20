

# with for a given granges and the bigwig, extract conservation score in region
get_score <- function(gr, bw_cons){
  
  problematic <- which(width(gr) == 1)
  
  res_problematic <- map_dbl(problematic,
                             ~ {
                               local_score <- subsetByOverlaps(bw_cons, gr[.x])$score
                               if(length(local_score) > 1) exit("problematic with length longer than 1")
                               if(length(local_score) == 0L) return(0)
                               local_score
                             })
  
  width(gr)[problematic] <- 2
  score <- genomation::ScoreMatrixBin(
    target = bw_cons,
    windows = gr,
    bin.num = 1, weight.col = "score"
  )
  res <- as.numeric(score@.Data)
  
  
  res[problematic] <- res_problematic
  res
}


#~ A3 ----



get_cons_a3 <- function(d_a3, bw_cons){
  a3_intron_gr <- GRanges(seqnames = d_a3$chr,
                          ranges = IRanges(start = d_a3$intron_start,
                                           end = d_a3$intron_end),
                          strand = d_a3$strand)
  
  
  a3_overhang_gr <- GRanges(seqnames = d_a3$chr,
                            ranges = IRanges(start = d_a3$overhang_start,
                                             end = d_a3$overhang_end),
                            strand = d_a3$strand)
  
  
  
  a3_intron_score <- get_score(a3_intron_gr, bw_cons)
  
  a3_overhang_score <- get_score(a3_overhang_gr, bw_cons)
  
  
  
  data.frame(event_id = d_a3$event_id,
             intron_conservation = as.numeric(a3_intron_score@.Data),
             overhang_conservation = as.numeric(a3_overhang_score@.Data))
  
}






# A5 ----
get_cons_a5 <- function(d_a5, bw_cons){
  a5_intron_gr <- GRanges(seqnames = d_a5$chr,
                          ranges = IRanges(start = d_a5$intron_start,
                                           end = d_a5$intron_end),
                          strand = d_a5$strand)
  
  
  a5_overhang_gr <- GRanges(seqnames = d_a5$chr,
                            ranges = IRanges(start = d_a5$overhang_start,
                                             end = d_a5$overhang_end),
                            strand = d_a5$strand)
  
  
  
  a5_intron_score <- get_score(a5_intron_gr, bw_cons)
  
  a5_overhang_score <- get_score(a5_overhang_gr, bw_cons)
  
  
  
  
  data.frame(event_id = d_a5$event_id,
             intron_conservation = as.numeric(a5_intron_score@.Data),
             overhang_conservation = as.numeric(a5_overhang_score@.Data))
  
}




# AF ----

get_cons_af <- function(d_af, bw_cons){
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
  
  
  
  
  af_prox_exon_score <- get_score(af_prox_exon_gr, bw_cons)
  
  af_prox_intron_score <- get_score(af_prox_intron_gr, bw_cons)
  
  af_dist_exon_score <- get_score(af_dist_exon_gr, bw_cons)
  
  af_dist_intron_score <- get_score(af_dist_intron_gr, bw_cons)
  
  
  
  
  data.frame(event_id = d_af$event_id,
             proximal_exon_conservation = af_prox_exon_score,
             proximal_intron_conservation = af_prox_intron_score,
             distal_exon_conservation = af_dist_exon_score,
             distal_intron_conservation = af_dist_intron_score
  )
}





# AL ----

get_cons_al <- function(d_al, bw_cons){
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
  
  
  
  al_p_ex_seq <- get_score(al_prox_ex_gr, bw_cons)
  al_p_i_seq <- get_score(al_prox_in_gr, bw_cons)
  al_d_e_seq <- get_score(al_dist_ex_gr, bw_cons)
  al_d_i_seq <- get_score(al_dist_in_gr, bw_cons)
  
  
  
  
  data.frame(event_id = d_al$event_id,
             proximal_exon_conservation = al_p_ex_seq,
             
             proximal_intron_conservation = al_p_i_seq,
             
             distal_exon_conservation = al_d_e_seq,
             
             distal_intron_conservation = al_d_i_seq
  )
}




# MX ----

get_cons_mx <- function(d_mx, bw_cons){
  mx_f_e <- GRanges(seqnames = d_mx$chr,
                    ranges = IRanges(start = d_mx$first_exon_start,
                                     end = d_mx$first_exon_end),
                    strand = d_mx$strand)
  
  mx_f_ui <- GRanges(seqnames = d_mx$chr,
                     ranges = IRanges(start = d_mx$first_upstream_intron_start,
                                      end = d_mx$first_upstream_intron_end),
                     strand = d_mx$strand)
  mx_f_di <- GRanges(seqnames = d_mx$chr,
                     ranges = IRanges(start = d_mx$first_downstream_intron_start,
                                      end = d_mx$first_downstream_intron_end),
                     strand = d_mx$strand)
  
  mx_s_e <- GRanges(seqnames = d_mx$chr,
                    ranges = IRanges(start = d_mx$second_exon_start,
                                     end = d_mx$second_exon_end),
                    strand = d_mx$strand)
  mx_s_ui <- GRanges(seqnames = d_mx$chr,
                     ranges = IRanges(start = d_mx$second_upstream_intron_start,
                                      end = d_mx$second_upstream_intron_end),
                     strand = d_mx$strand)
  mx_s_di <- GRanges(seqnames = d_mx$chr,
                     ranges = IRanges(start = d_mx$second_downstream_intron_start,
                                      end = d_mx$second_downstream_intron_end),
                     strand = d_mx$strand)
  
  mx_f_e_s <- get_score(mx_f_e, bw_cons)
  mx_f_ui_s <- get_score(mx_f_ui, bw_cons)
  mx_f_di_s <- get_score(mx_f_di, bw_cons)
  mx_s_e_s <- get_score(mx_s_e, bw_cons)
  mx_s_ui_s <- get_score(mx_s_ui, bw_cons)
  mx_s_di_s <- get_score(mx_s_di, bw_cons)
  
  
  
  data.frame(event_id = d_mx$event_id,
             first_exon_conservation = mx_f_e_s,
             first_upstream_intron_conservation = mx_f_ui_s,
             first_downstream_intron_conservation = mx_f_di_s,
             
             second_exon_conservation = mx_s_e_s,
             second_upstream_intron_conservation = mx_s_ui_s,
             second_downstream_intron_conservation = mx_s_di_s
  )
}




# RI ----


get_cons_ri <- function(d_ri, bw_cons){
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
  
  
  ri_de_s <- get_score(ri_ue_gr, bw_cons)
  ri_ue_s <- get_score(ri_de_gr, bw_cons)
  ri_i_s <- get_score(ri_intron_gr, bw_cons)
  
  
  
  data.frame(event_id = d_ri$event_id,
             downstream_exon_conservation = ri_de_s,
             
             upstream_exon_conservation = ri_ue_s,
             
             intron_conservation = ri_i_s
  )
}




# SE ----

get_cons_se <- function(d_se, bw_cons){
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


se_ui_seq <- get_score(se_ui, bw_cons)
se_di_seq <- get_score(se_di, bw_cons)
se_e_seq <- get_score(se_e, bw_cons)


data.frame(event_id = d_se$event_id,
           upstream_intron_conservation = se_ui_seq,

           downstream_intron_conservation = se_di_seq,

           exon_conservation = se_e_seq
             )
}


get_cons <- function(event_type, event_coordinates, bw_cons){
  switch (event_type,
          A3 = get_cons_a3(event_coordinates, bw_cons),
          A5 = get_cons_a5(event_coordinates, bw_cons),
          AL = get_cons_al(event_coordinates, bw_cons),
          AF = get_cons_af(event_coordinates, bw_cons),
          MX = get_cons_mx(event_coordinates, bw_cons),
          RI = get_cons_ri(event_coordinates, bw_cons),
          SE = get_cons_se(event_coordinates, bw_cons)
  )
}











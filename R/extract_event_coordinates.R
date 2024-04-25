


# For all events, calling c1, c2, ... the first, second,... coordinate given
# See here for correspondence to exon borders:
# https://github.com/comprna/SUPPA#generation-of-transcript-events-and-local-alternative-splicing-events
#
# here, schematics of exon borders: `###` constitutive exon, `---` alt exon, `...` intron





#dispatching function
extract_coords <- function(event_type, event_coordinates){
  switch (event_type,
    A3 = extract_coords_a3(event_coordinates),
    A5 = extract_coords_a5(event_coordinates),
    AL = extract_coords_al(event_coordinates),
    AF = extract_coords_af(event_coordinates),
    MX = extract_coords_mx(event_coordinates),
    RI = extract_coords_ri(event_coordinates),
    SE = extract_coords_se(event_coordinates)
  )
}


  
  

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



extract_coords_a3 <- function(event_coordinates){
  event_coordinates |>
    enframe(name = NULL,
            value = "event_coordinates") |>
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
}
  





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



extract_coords_a5 <- function(event_coordinates){
  event_coordinates |>
    enframe(name = NULL,
            value = "event_coordinates") |>
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
}








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



extract_coords_af <- function(event_coordinates){
  
  event_coordinates |>
    enframe(name = NULL,
            value = "event_coordinates") |>
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
}




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



extract_coords_al <- function(event_coordinates){
  
  event_coordinates |>
    enframe(name = NULL,
            value = "event_coordinates") |>
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
}







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


extract_coords_mx <- function(event_coordinates){
  
  event_coordinates |>
    enframe(name = NULL,
            value = "event_coordinates") |>
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
}




#~ RI ----

#'  +  #########------------------############
#     s1       e1                s2          e2
#     c1       c2                c3          c4
#
# Minus strand identical
#
#
#




extract_coords_ri <- function(event_coordinates){
  
  event_coordinates |>
    enframe(name = NULL,
            value = "event_coordinates") |>
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
}







#~ SE ----

#'  +  ########...........-----------.........############
#            e1          s2          e2       s3
#            c1          c2          c3       c4
#
# Minus strand identical
#
#
#



extract_coords_se <- function(event_coordinates){
  
  event_coordinates |>
    enframe(name = NULL,
            value = "event_coordinates") |>
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
}






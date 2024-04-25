# Inits ----

library(tidyverse)
source("R/utils.R")

library(printMat)
library(wbData)

tx2g <- wb_load_tx2gene(289)
gids <- wb_load_gene_ids(289)

# export_dir <- "data/outs/240305_fig"

outlier_samples <- c("AVKr113","RICr133","PVMr122", "ASIr154")




#~ Load ----

# filter from Alec's threshold
gene_expression_table <- read.delim("../majiq/data/2024-03-05_alec_integration/bsn12_subtracted_integrated_binarized_expression_withVDDD_FDR0.05_030424.tsv") |>
  as.data.frame() |>
  rownames_to_column("gene_id") |>
  pivot_longer(-gene_id,
               names_to = "neuron_id",
               values_to = "expressed") |>
  mutate(is_expressed = expressed == 1L) |> select(-expressed)


neurons_here <- unique(gene_expression_table$neuron_id)
genes_with_known_expr <- unique(gene_expression_table$gene_id)




psi_lg <- read.delim("data/240301b_psiPerEvent.psi") |>
  rownames_to_column("event_id") |>
  as_tibble() |>
  separate_wider_regex(event_id,
                       patterns = c(gene_id = "^WBGene[0-9]{8}", ";",
                                    event_type = "[SEA53MXRIFL]{2}", "\\:",
                                    event_coordinates = "[IXV]+\\:[0-9:\\-]+:[+-]$"),
                       cols_remove = FALSE) |>
  pivot_longer(-c(gene_id, event_type, event_coordinates, event_id),
               names_to = "sample_id",
               values_to = "PSI") |>
  filter(! sample_id %in% outlier_samples) |>
  mutate(neuron_id = str_match(sample_id, "^([A-Zef0-9]{2,4})r[0-9]{1,4}")[,2]) |>
  filter(neuron_id %in% neurons_here) |>
  left_join(gene_expression_table |>
              rename(expressed = is_expressed),
            by = c("gene_id", "neuron_id"))



#~ functions ----

rename_event <- function(nm){
  descriptor <- str_split_i(nm, ":",1)
  gene_name <- str_split_i(descriptor, ";",1) |> i2s(gids, warn_missing = TRUE)
  type <- str_split_i(descriptor, ";",2)
  paste0(type,"_",gene_name)
}


# among samples where they are both expressed, how often do they have the same inclusion
jaccard <- function(x,y){
  1 - sum(x == y,  na.rm = TRUE) / sum(! is.na(x&y))
}

# matrix of distances
get_jac_mat <- function(psi){
  
  deltapsi <- apply(psi,1, \(x) x - mean(x, na.rm = TRUE)) |> t()
  psi_bin <- (deltapsi > 0)
  
  res <- matrix(nrow = nrow(psi_bin), ncol = nrow(psi_bin))
  for(i in 1:nrow(psi_bin)){
    for(j in i:nrow(psi_bin)){
      res[i,j] <- jaccard(psi_bin[i, ],
                          psi_bin[j, ])
      res[j,i] <- res[i,j]
    }
  }
  dimnames(res) <- list(rownames(psi_bin),
                        rownames(psi_bin))
  res
}



# Clustering on PSI ----


#~ Pre-process ----


# by sample
psi_mat <- psi_lg |>
  filter(expressed) |>
  # filter(event_type == "AF") |>
  pivot_wider(id_cols = event_id,
              names_from = sample_id,
              values_from = PSI) |>
  column_to_rownames("event_id") |>
  as.matrix()


# filter
prefilt_psi <- psi_mat

xx <- rowSums(is.na(prefilt_psi))
hist(xx)
table(xx <= 100)
prefilt_psi <- prefilt_psi[xx < 100,]

xx <- matrixStats::rowSds(prefilt_psi, na.rm = TRUE)
hist(xx)
table(xx > .2)
prefilt_psi <- prefilt_psi[xx > .2,]







# by neuron
psi_neur <- psi_lg |>
  filter(expressed) |>
  # filter(event_type == "AF") |>
  summarize(PSI = mean(PSI, na.rm = TRUE),
            .by = c(event_id, neuron_id)) |>
  pivot_wider(id_cols = event_id,
              names_from = neuron_id,
              values_from = PSI) |>
  column_to_rownames("event_id") |>
  as.matrix()


prefilt_neur <- psi_neur[rownames(prefilt_psi),]

# more convenient rownames
rownames(prefilt_psi) <- paste0(1:nrow(prefilt_psi), "_",
                                rename_event(rownames(prefilt_psi)))
rownames(prefilt_neur) <- paste0(1:nrow(prefilt_neur), "_",
                                rename_event(rownames(prefilt_neur)))

all.equal(row.names(prefilt_psi), row.names(prefilt_neur))


dim(prefilt_psi)
dim(prefilt_neur)

matimage(prefilt_psi)




#~ Jaccard distances ----



# color all classes
color_coding <- tibble(sample_id = colnames(prefilt_psi),
       neuron_id = str_split_i(sample_id, "r", 1),
       color = rep(ggsci::pal_d3("category20")(20),3)[as.integer(factor(neuron_id))] ) |>
  column_to_rownames("sample_id")
# color some
color_coding <- tibble(sample_id = colnames(prefilt_psi),
                       neuron_id = str_split_i(sample_id, "r", 1)) |>
  filter(neuron_id %in% c("AVA","AVE","DVC", "OLL")) |>
  mutate(color = rep(ggsci::pal_d3("category20")(20),3)[as.integer(factor(neuron_id))] ) |>
  column_to_rownames("sample_id")


# color neur cluster
color_neur <- tibble(neur_id = colnames(prefilt_neur)) |>
  mutate(color = case_match(
    neur_id,
    c("I5","PHA","OLQ","PVC","LUA","OLL") ~ "#1F77B4FF",
    c("AVH","AVK") ~ "#FF7F0EFF",
    c("AVA","AVE") ~ "#2CA02CFF",
    c("DD","VD") ~ "#D62728FF",
    c("DA","VB") ~ "#9467BDFF",
    c("ASK", "ADF", "ASG", "ASEL", "ASER", "AFD", "BAG", "AWA", "AIY", "AWB") ~ "#8C564BFF"
  )) |>
  column_to_rownames("neur_id")



#~~ Distances ----

mjac_ev <- get_jac_mat(prefilt_psi)
djac_ev <- as.dist(mjac_ev)
plot(hclust(djac_ev), labels = FALSE)


mjac_sample <- get_jac_mat(t(prefilt_psi))
djac_sample <- as.dist(mjac_sample)
plot(hclust(djac_sample))

mjac_neur <- get_jac_mat(t(prefilt_neur))
djac_neur <- as.dist(mjac_neur)
plot(hclust(djac_neur))





par(cex = .7)
#~~ optimize hclust sample ----
hc_samp <- djac_sample |>
  hclust(method = "ward.D")
dend <- as.dendrogram(hc_samp)

dendextend::labels_colors(dend) <- color_coding[hc_samp$labels, "color"][hc_samp$order]
stats:::plot.dendrogram(dend)




par(cex=.8)
#~~ optimize hclust neur ----
for(method in c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")){
  hc_neur <- djac_neur |>
    hclust(method = method)
  dend <- as.dendrogram(hc_neur)
  
  dendextend::labels_colors(dend) <- color_neur[hc_neur$labels, "color"][hc_neur$order]
  # png(paste0("intermediates/240424_clust/", method,".png"), width = 650, height = 600)
  # stats:::plot.dendrogram(dend, main = method)
  # dev.off()
}

hc_neur <- djac_neur |>
  hclust(method = "mcquitty")
dend <- as.dendrogram(hc_neur)

dendextend::labels_colors(dend) <- color_neur[hc_neur$labels, "color"][hc_neur$order]
stats:::plot.dendrogram(dend)


#~~ How to optimize hc_events? ----


hc_events <- djac_ev |>
  hclust(method = "ward.D")

plot(hc_events, labels = FALSE)






pheatmap::pheatmap(prefilt_psi,
                   show_rownames = FALSE,
                   cluster_rows = hc_events,
                   cluster_cols = hc_samp,
                   cutree_rows = 10,
                   cutree_cols = 10)

pheatmap::pheatmap(prefilt_neur,
                   show_rownames = FALSE,
                   cluster_rows = hc_events,
                   cluster_cols = hc_neur,
                   cutree_rows = 10,
                   cutree_cols = 5)





#~ average by cluster and sample ----

ct_samp <- cutree(hc_samp, k=10)
table(ct_samp)

ct_ev <- cutree(hc_events, k = 8)
table(ct_ev)

# # reorder PSI to be highest in I5
# ordered_psi <- prefilt_psi
# to_invert <- which(ordered_psi[,"AVAr22"] < .5)
# ordered_psi[to_invert,] <- 1 - ordered_psi[to_invert,]

# matimage(ordered_psi)

averaged_psi <- matrix(-10, nrow = length(unique(ct_ev)), ncol = ncol(prefilt_psi))
dimnames(averaged_psi) <- list(unique(ct_ev),
                               colnames(prefilt_psi))
for(cl in unique(ct_ev)){
  submat <- prefilt_psi[ct_ev == cl,]
  averaged_psi[cl,] <- colMeans(submat, na.rm = TRUE)
}
# matimage(averaged_psi)
# head(averaged_psi)

pheatmap::pheatmap(averaged_psi,
                   cutree_cols = 10,fontsize = 5)





#~ average by cluster and neuron ----

ct_neur <- cutree(hc_neur, k=10)
table(ct_neur)

ct_ev <- cutree(hc_events, k = 8)
table(ct_ev)

# # reorder PSI to be highest in I5
# ordered_neur <- prefilt_neur
# to_invert <- which(prefilt_neur[,"I5"] < .5)
# ordered_neur[to_invert,] <- 1 - ordered_neur[to_invert,]

matimage(ordered_neur)

averaged_psi <- matrix(-10, nrow = length(unique(ct_ev)), ncol = ncol(prefilt_neur))
dimnames(averaged_psi) <- list(unique(ct_ev),
                               colnames(prefilt_neur))
for(cl in unique(ct_ev)){
  submat <- prefilt_neur[ct_ev == cl,]
  averaged_psi[cl,] <- colMeans(submat, na.rm = TRUE)
}
matimage(averaged_psi)
head(averaged_psi)

pheatmap::pheatmap(averaged_psi,
                   cutree_cols = 5)


#~ check some clusters ----

# cl9 in CAN vs all others

single_cl <- prefilt_psi[ct_ev == 2,]
# just_cl9 <- just_cl9[,colSums(is.na(just_cl9)) < 4]

rownames(single_cl) <- paste0(1:nrow(single_cl), "_",
                              rownames(single_cl) |> str_split_i(";",1) |> i2s(gids))
pheatmap::pheatmap(single_cl,
                   cluster_rows = FALSE,
                   cluster_cols = hc_samp,
                   cutree_cols = 10,
                   fontsize = 7)

annots <- tibble(sample_id = colnames(single_cl),
                neuron_id = str_split_i(sample_id, "r",1),
                is_this_cl = if_else(neuron_id %in% names(ct_neur)[ct_neur == 2],
                                     "yes","no")) |>
  select(sample_id, is_this_cl) |>
  column_to_rownames("sample_id")
pheatmap::pheatmap(single_cl,
                   show_rownames = FALSE,
                   cluster_rows = FALSE,
                   cluster_cols = hc_samp,
                   annotation_col = annots)

single_cl_neur <- prefilt_neur[ct_ev == 2,]
rownames(single_cl_neur) <- paste0(1:nrow(single_cl_neur), "_",
                              rownames(single_cl_neur) |> str_split_i(";",1) |> i2s(gids))
pheatmap::pheatmap(single_cl_neur,
                   show_rownames = TRUE,
                   cluster_rows = FALSE,
                   cluster_cols = hc_neur)










######


cc <- cor(t(prefilt_psi), use = "pairwise")

dim(cc)


xx <- matrixStats::rowSds(cc, na.rm = TRUE)
hist(xx)
table(xx > .5)
cc <- cc[xx > .5,xx > .5]


xx <- rowSums(is.na(cc))
hist(xx)
table(xx <= 2500)
cc <- cc[xx < 2500,xx < 2500]

cc_noname <- cc; dimnames(cc_noname) <- NULL

heatmap3::heatmap3(cc_noname,
                   Rowv = NA,
                   Colv = NA,
                   labRow = NA,
                   labCol = NA)

image(t(cc_noname), useRaster = TRUE)
image(cc_noname[, nrow(cc_noname):1], useRaster = TRUE, axes = FALSE)




cdist <- dist(cc)

hc <- hclust(cdist)
plot(hc,labels = FALSE)

ct <- kmeans(impute::impute.knn(cc)$data, centers = 10)$cluster

ct <- cutree(hc, k=15)

head(ct)
table(ct)


heatmap3::heatmap3(cc_noname[ct == 1, ct == 1],
                   Rowv = NA,
                   Colv = NA,
                   labRow = NA,
                   labCol = NA)
heatmap3::heatmap3(cc_noname[ct == 2, ct == 2],
                   Rowv = NA,
                   Colv = NA,
                   labRow = NA,
                   labCol = NA)

events_by_clust <- enframe(ct,
        name = "event_id",
        value = "cluster") |>
  separate_wider_regex(event_id,
                       patterns = c(gene_id = "^WBGene[0-9]{8}", ";",
                                    event_type = "[SEA53MXRIFL]{2}", "\\:",
                                    event_coordinates = "[IXV]+\\:[0-9:\\-]+:[+-]$"),
                       cols_remove = FALSE) |>
  mutate(gene_name = i2s(gene_id, gids, warn_missing = TRUE), .after = "gene_id") |>
  mutate(cluster = factor(cluster))


events_by_clust |>
  ggplot() +
  theme_classic() +
  geom_bar(aes(x = cluster, fill = event_type),
           position = "dodge")

events_by_clust |>
  # mutate(cluster = sample(cluster)) |>
  ggplot() +
  theme_classic() +
  geom_bar(aes(x = event_type, fill = cluster),
           position = "fill")

events_by_clust |>
  summarize(n = n(),
            .by = c(gene_id, cluster)) |>
  summarize(nb_per_clust = mean(n),
            nb_diff_clust = length(unique(cluster)),
            tot_events = sum(n),
            .by = gene_id) |>
  filter(tot_events > 1) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x = nb_diff_clust, y = nb_per_clust), alpha = .2)

events_by_clust |>
  mutate(cluster = sample(cluster)) |>
  summarize(n = n(),
            .by = c(gene_id, cluster)) |>
  summarize(nb_per_clust = mean(n),
            nb_diff_clust = length(unique(cluster)),
            tot_events = sum(n),
            .by = gene_id) |>
  filter(tot_events > 1) |>
  ggplot() +
  theme_classic() +
  ggbeeswarm::geom_quasirandom(aes(x = nb_diff_clust, y = nb_per_clust), alpha = .2)


annot <- psi_lg |> select(sample_id, neuron_id) |> distinct() |> column_to_rownames("sample_id")




neur_hc <- hclust(dist(cor(t(prefilt_psi), use = "pairwise"),
                       method = "maximum"))


pheatmap::pheatmap(t(prefilt_psi),
                   show_colnames = FALSE,
                   cluster_rows = neur_hc,
                   cluster_cols = hc,
                   annotation_row = annot,
                   annotation_legend = FALSE,
                   fontsize = 5)

submat <- function(ck){
  x <- t(prefilt_psi)[,ct == ck]
  x[rowSums(x, na.rm = TRUE) > 0,]
}

cl1 <- events_by_clust |> filter(cluster == 1) |> pull(event_id)
pheatmap::pheatmap(t(psi_mat)[,cl1],
                   show_colnames = FALSE,
                   cluster_rows = TRUE,
                   cluster_cols = FALSE)


pheatmap::pheatmap(submat(2),
                   show_colnames = FALSE,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE)

cl5 <- events_by_clust |> filter(cluster == 5) |> pull(event_id)
pheatmap::pheatmap(t(psi_mat)[,cl5],
                   show_colnames = FALSE,
                   cluster_cols = FALSE)

cl6 <- events_by_clust |> filter(cluster == 6) |> pull(event_id)
pheatmap::pheatmap(t(psi_mat)[,cl6],
                   show_colnames = FALSE,
                   cluster_cols = FALSE)




# ____________ ----
hclust(dist(cc, method = "euclidean")) |> plot()

pheatmap::pheatmap(cc,
                   show_rownames = FALSE,
                   show_colnames = FALSE)



image(psi_mat)



dd_euc <- dist(t(psi_mat), method = "euclidean") |> as.matrix()
dd_can <- dist(t(psi_mat), method = "canberra") |> as.matrix()

dd_euc[c("ASEL","ASER","AWA"),c("ASEL","ASER","AWA")]
dd_can[c("ASEL","ASER","AWA"),c("ASEL","ASER","AWA")]

dd_euc[c("AVA","AVE","DVB"),c("AVA","AVE","DVB")]
dd_can[c("AVA","AVE","DVB"),c("AVA","AVE","DVB")]



plot(psi_mat[,"ASEL"],
     psi_mat[,"ASER"])

plot(psi_mat[,"ASEL"],
     psi_mat[,"AWA"])

plot(psi_mat[,"AVA"],
     psi_mat[,"AVE"])

plot(psi_mat[,"AVA"],
     psi_mat[,"DVB"])

cc <- cor(psi_mat, use = "pairwise")
pheatmap::pheatmap(cc)


hclust(dist(t(psi_mat), method = "euclidean")) |> plot()
hclust(dist(t(psi_mat), method = "canberra")) |> plot()


plot(psi_mat[,"AIY"],
     psi_mat[,"AVG"])

plot(psi_mat[,"AIY"],
     psi_mat[,"DVB"])


km <- kmeans(psi_mat, centers = 10)

dd <- cluster::daisy(t(psi_mat))

as.matrix(dd)[1:3,1:3]

as.matrix(dd)[c("ASEL","ASER","AWA"),c("ASEL","ASER","AWA")]

as.matrix(dd)[c("AVA","AVE","DVB"),c("AVA","AVE","DVB")]

hclust(dd) |> plot()

cluster::agnes(dd, method = "ward") |> plot()

dd <- cluster::daisy(psi_mat)
cluster::agnes(dd, method = "ward") |> plot()


dbs <- dbscan::dbscan(dd_euc, eps = 9)
plot(dbs)
dbscan::hullplot(dbs)

dbscan::optics(psi_mat)



# remove events that are detectable in too few samples
xx <- rowSums(is.na(psi_mat))
hist(xx)
table(xx <= 30)
psi_mat <- psi_mat[which(rowSums(is.na(psi_mat)) <= 20),]

yy <- matrixStats::rowVars(psi_mat, na.rm = TRUE)
hist(yy)
table(yy > .05)

psi_mat <- psi_mat[which(yy > .05),]



hc_rows <- stats::hclust(dist(psi_mat, method = "canberra"), method = "complete")
hc_cols <- hclust(dist(t(psi_mat), method = "canberra"), method = "ward.D2")

pheatmap::pheatmap(psi_mat,
                   show_rownames = FALSE,
                   cluster_rows = FALSE,
                   # cluster_cols = hc_cols,
                   cutree_rows = 3)



plot(psi_mat[,"ASEL"],
     psi_mat[,"ASER"])

plot(psi_mat[,"ASEL"],
     psi_mat[,"AWA"])

plot(psi_mat[,"AVA"],
     psi_mat[,"AVE"])


#~ test Hellinger correlation ----
# mat <- t(psi_mat)
# 
# Hcor <- matrix(nrow = ncol(mat), ncol = ncol(mat))
# for(i in seq_len(ncol(mat))){
#   for(j in seq_len(ncol(mat))){
#     x <- mat[,i]
#     y <- mat[,j]
#     obs <- which(!is.na(x) & !is.na(y))
#     Hcor[i,j] <- HellCor::HellCor(x[obs], y[obs])$Hcor
#   }
# }
# 
# dimnames(Hcor) <- list(colnames(mat),
#                        colnames(mat))
# 


spearcor <- cor(mat, use = "pairwise", method = "spearman")


hc_hcor <- stats::hclust(as.dist(1 - spearcor), method = "complete")

plot(hc_hcor,labels = FALSE)

pheatmap::pheatmap(t(mat),
                   show_rownames = FALSE,
                   cluster_rows = hc_hcor,
                   cluster_cols = hc_cols,
                   annotation_col = annot_df)

annot_df <- read.csv("../majiq/data/neuron_properties.csv") |>
  select(-include) |>
  column_to_rownames("Neuron_type")



# MDS ----
cor_mat <- cor(t(psi_mat), use = "pairwise.complete.obs")
cor_mat[1:3,1:3]

diag(cor_mat) <- 0
mds_ds <- cmdscale(cor_mat, eig = TRUE)

ggplot(tibble(eig = mds_ds$eig[seq_len(100)], k = seq(along = eig)),
       aes(x = k, y = eig)) + theme_minimal() +
  scale_x_discrete("k", limits = as.factor(seq_len(100))) + 
  geom_bar(stat = "identity", width = 0.5, fill = "#ffd700", col = "#0057b7")

mds_ds_df <- mds_ds$points |>
  as.data.frame() |>
  rownames_to_column("event_id") |>
  mutate(gene_id = str_split_i(event_id, ";", 1),
         gene_name = i2s(gene_id, gids, warn_missing = TRUE))

ggplot(mds_ds_df,
       aes(x = V1, y = V2, label = gene_name)) +
  geom_point() +
  ggrepel::geom_text_repel(col = "#0057b7") + coord_fixed() 



#~ Cluster by PCA ----

mat_psi_imp <- impute::impute.knn(psi_mat)$data

pca <- prcomp(mat_psi_imp)

eigs <- pca$sdev^2
plot(cumsum(eigs)/sum(eigs), ylim = c(0,1))





km <- kmeans(mat_psi_imp, 5)

pca$x[,1:2] |>
  as.data.frame() |>
  rownames_to_column("event_id") |>
  mutate(gene_id = str_split_i(event_id, ";", 1),
         gene_name = i2s(gene_id, gids, warn_missing = TRUE),
         cluster = factor(km$cluster)) |>
  ggplot(aes(x = PC1, y = PC2, label = gene_name, color = cluster)) +
  geom_point() +
  ggrepel::geom_text_repel(col = "#0057b7") + coord_fixed() 




# WGCNA ----

GS1 <- as.numeric(cor_mat)

pvals <- WGCNA::corPvalueFisher(GS1, nSamples = ncol(psi_mat))

table(is.na(pvals))
qvals <- WGCNA::qvalue(pvals)$qvalues

StandardGeneScreeningResults <- data.frame(event_id = rownames(cor_mat),
                                           PearsonCorrelation = GS1,
                                           pvals,
                                           qvals)

table(qvals < .05)

# soft thresholding
adj1 <- abs(cor(t(psi_mat), use = "pairwise"))^5

k <- WGCNA::softConnectivity(datE=t(psi_mat),power=5)


par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

datExpr <- t(psi_mat)
dissADJ <- 1 - adj1

dissTOM <- TOMdist(adj1)

hierADJ <- hclust(as.dist(dissADJ),
                  method="average" )

plotDendroAndColors(hierADJ,
                    colors = rep("black", 772),
                    dendroLabels = FALSE, hang = 0.03,
                    main = "Gene hierarchical clustering dendrogram" )


colorStaticADJ <- as.character(cutreeStaticColor(hierADJ,
                                                 cutHeight=.99,
                                                 minSize=10))


plotDendroAndColors(hierADJ, colors = data.frame(colorStaticADJ),
                    dendroLabels = FALSE, abHeight = 0.99,
                    main = "Gene dendrogram and module colors")


branch.number=cutreeDynamic(hierADJ,method="tree")
# This function transforms the branch numbers into colors
colorDynamicADJ=labels2colors(branch.number )




powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

sft = pickSoftThreshold(
  t(psi_mat),             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)



par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")



# run!

cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(psi_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = 5,                # <= power here
                          networkType = "unsigned",
                          
                          # == Tree and Block Options ==
                          deepSplit = 4,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

library(dplyr)

matrixFromTreeXY <- function(wg_stat_df) {
  all_taxa <- unique(c(wg_stat_df$pop1, wg_stat_df$pop2))
  distances_df <- data.frame(matrix(ncol = length(seq(1, length(all_taxa))), nrow = length(seq(1, length(all_taxa)))))
  rownames(distances_df) <- seq(1, length(all_taxa))
  colnames(distances_df) <- seq(1, length(all_taxa))
  distances_df[is.na(distances_df)] <- 0
  for(seq in rownames(distances_df)) {
    for(seq2 in rownames(distances_df)){
      if(seq != seq2){
        distances_df[seq,seq2] <- wg_stat_df[wg_stat_df$pop1 == seq & wg_stat_df$pop2 == seq2 | wg_stat_df$pop1 == seq2 & wg_stat_df$pop2 == seq,c(1)]
      }
    }
  }
  return(distances_df)
}

#####################
#### WG dXY tree ####
#####################

# path to the directory containing the output CSV files from the treeXY analysis
comp_dir <- "/path/to/treeXY/output/"

# read treeXY files
# make sure to change the pattern matched by dir if your treeXY filenames are different
read_scaffs_list <- list()
for(window_file in dir(comp_dir, pattern = "Chr[1,2,3,4,5,6,7,8]_treeXY_filtered_15_200_2_2_treeXY.csv")) {
  read_scaffs_list <- append(read_scaffs_list, list(read.csv(paste(comp_dir, window_file, sep = "/"), header=T)))
}

in_treeXY <- bind_rows(read_scaffs_list, .id = "column_label")

WG_dXY <- in_treeXY[,grep("dXY", colnames(in_treeXY))]
WG_dXY <- as.data.frame(colMeans(WG_dXY))
WG_dXY$pop1 <- unlist(lapply(rownames(WG_dXY), function(x) strsplit(x, "_")[[1]][2]))
WG_dXY$pop2 <- unlist(lapply(rownames(WG_dXY), function(x) strsplit(x, "_")[[1]][3]))

distances_df <- matrixFromTreeXY(WG_dXY)
dend_mean_dXY  <- as.dendrogram(hclust(d = as.dist(distances_df), method = "average"))

plot(dend_mean_dXY)

###################
#### WG D tree ####
###################

WG_D <- in_treeXY[,grep("D", colnames(in_treeXY))]
WG_D <- as.data.frame(colMeans(WG_D))

WG_D$pop1 <- lapply(rownames(WG_D), function(x) strsplit(x, "_")[[1]][2])
WG_D$pop2 <- lapply(rownames(WG_D), function(x) strsplit(x, "_")[[1]][3])

distances_df <- matrixFromTreeXY(WG_D)
dend_mean_D  <- as.dendrogram(hclust(d = as.dist(distances_df), method = "average"))

plot(dend_mean_D)


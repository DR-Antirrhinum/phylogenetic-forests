library(dplyr)
library(cluster)

options(scipen=999)

###################
#### Functions ####
###################

# measure shortest root branch, and populate a table of statistics for each tree
getClusterSummaries <- function(distanceDF, agnesList) {
  clusterSummaries <- as.data.frame(matrix(nrow = nrow(na.omit(distanceDF)), ncol = 6))
  rownames(clusterSummaries) <- rownames(na.omit(distanceDF))
  colnames(clusterSummaries) <- c("cluster_one_size", "cluster_two_size", "tree_height", "longest_root_branch_length", "shortest_root_branch_length", "cophenetic_correlation")
  for(window in names(agnesList)) {
    agnesDxy <- agnesList[[window]][[1]]
    cophCorDxy <- cor(agnesDxy[[5]], cophenetic(agnesDxy))
    cophMatDxy <- as.matrix(cophenetic(agnesDxy))
    heightDxy <- max(cophMatDxy)
    # Find first cluster, and count number of leaves
    clusterTwo <- max(cophMatDxy[cophMatDxy != max(cophMatDxy)])
    clusterTwoPops <- names(apply(cophMatDxy, 1, function(x) clusterTwo %in% x)
                            [apply(cophMatDxy, 1, function(x) clusterTwo %in% x) == TRUE])
    clusterOnePops <- rownames(cophMatDxy)[!rownames(cophMatDxy) %in% clusterTwoPops]
    clusterOne <- max(cophMatDxy[rownames(cophMatDxy) %in% clusterOnePops,colnames(cophMatDxy) %in% clusterOnePops])
    # N.B. these might be the wrong way round...
    longestRoot <- heightDxy - clusterOne
    shortestRoot <- heightDxy - clusterTwo
    # Find second cluster, and count number of leaves
    clusterTwoPops <- rownames(cophMatDxy)[!rownames(cophMatDxy) %in% clusterOnePops]
    
    # Collect summary stats
    clusterSummaries[window,] <- c(length(clusterOnePops), 
                                   length(clusterTwoPops), 
                                   heightDxy, 
                                   longestRoot,
                                   shortestRoot,
                                   cophCorDxy)
  }
  return(clusterSummaries)
}

##############
#### Main ####
##############

# chromosome details
lg_stats <- read.table("Chromosome_offset_tab.txt", header=T, sep="\t")
# calculate chromosome offsets for WG plotting
lg_offsets <- cumsum(as.numeric(lg_stats$scaffold_length))
lg_offsets <- c(0, lg_offsets[-c(length(lg_offsets))])
lg_stats$offset_whole_genome <- lg_offsets

# treeXY files
# make sure to edit the dir pattern if your treeXY output files are not in the same directory as this script
read_scaffs_list <- list()
for(window_file in dir(pattern = "*_15_200_2_18_50000_25000_treeXY.csv")) {
  read_scaffs_list <- append(read_scaffs_list, list(read.csv(paste(window_file, sep = "/"), header=T)))
}

in_treeXY <- bind_rows(read_scaffs_list, .id = "column_label")
Filtered_treeXY <- in_treeXY[in_treeXY$n_window_sites > 5000,]
# Filtered_treeXY <- in_treeXY[in_treeXY$n_window_sites > 1000,]
treeXY_off <- merge(Filtered_treeXY, lg_stats, by.x = "scaffold", by.y = "scaffold", all.x=T)
treeXY_off$window_midkb <- ((treeXY_off$window_start + treeXY_off$window_end) / 2000)
treeXY_off$window_midMb<-((treeXY_off$window_start + treeXY_off$window_end) / 2000000)
treeXY_off$linkage_midkb <- treeXY_off$window_midkb + (treeXY_off$offset_LG/1000)
treeXY_off$linkage_midMb <- treeXY_off$linkage_midkb/1000
treeXY_off$chromosome_midkb <- treeXY_off$window_midkb + (treeXY_off$offset_whole_genome/1000)
treeXY_off$chromosome_midMb <- treeXY_off$chromosome_midkb/1000

# population / individual details
# change these names if you're working with a different set of populations
pop_names <- c("UNA", "BED", "LU", "AXA", "MIJ", "MON", "PER", "BOU", "VIL", "ARS", "THU", "BAN", "ARL", "CIN", "YP4", "YP1", "MP4", "MP11")
pop_species <- c("ps", "ps", "st", "st", "st", "st", "ps", "st", "st", "ps", "st", "ps", "ps", "ps", "st", "st", "ps", "ps")
pop_tab <- cbind.data.frame(pop_names, pop_species)

long_piw <- treeXY_off[,grep("piw", colnames(treeXY_off))]
rownames(long_piw) <- treeXY_off$chromosome_midMb
numbered_pop_names <- pop_names
names(numbered_pop_names) <- seq(1, length(pop_names))
colnames(long_piw) <- unlist(lapply(strsplit(colnames(long_piw), split = "_"), function(x) paste(x[1], numbered_pop_names[[x[2]]], sep = "_")))

long_piT <- treeXY_off[,grep("piT", colnames(treeXY_off))]
rownames(long_piT) <- treeXY_off$chromosome_midMb
numbered_pop_names <- pop_names
names(numbered_pop_names) <- seq(1, length(pop_names))
colnames(long_piT) <- unlist(lapply(strsplit(colnames(long_piT), split = "_"), function(x) paste(x[1], numbered_pop_names[[x[2]]], numbered_pop_names[[x[[3]]]], sep = "_")))

# make dxy trees for each row
long_dXY <- treeXY_off[,grep("dXY", colnames(treeXY_off))]
rownames(long_dXY) <- treeXY_off$chromosome_midMb
numbered_pop_names <- pop_names
names(numbered_pop_names) <- seq(1, length(pop_names))
colnames(long_dXY) <- unlist(lapply(strsplit(colnames(long_dXY), split = "_"), function(x) paste(x[1], numbered_pop_names[[x[2]]], numbered_pop_names[[x[[3]]]], sep = "_")))

WG_dXY_variance_t <- as.data.frame(t(long_dXY))
splitPops <- strsplit(as.character(rownames(WG_dXY_variance_t)), c("_"), fixed = TRUE)
tSplitPops <- t(as.data.frame(splitPops))
tSplitPops <- as.data.frame(tSplitPops[,-c(1)])
colnames(tSplitPops) <- c("Pop1", "Pop2")
UniquePops <- unique(c(paste(tSplitPops[,1]), paste(tSplitPops[,2])))

all_dXY_edges <- cbind(as.data.frame(tSplitPops), as.data.frame(WG_dXY_variance_t))

all_DMs_dXY <- list()
all_hclusts_dXY <- list()
all_agnes_dXY <- list()
for (i in rownames(na.omit(long_dXY))){
  locus <- paste("\\b^", i, "$\\b", sep = "")
  edges_for_bimodal_test <- cbind.data.frame(all_dXY_edges[,c(1,2)], all_dXY_edges[,-c(1,2)][,grep(locus, colnames(all_dXY_edges)[-c(1,2)])[1]])
  window_df <- data.frame(matrix(nrow = length(UniquePops), ncol = length(UniquePops)))
  rownames(window_df) <- UniquePops
  colnames(window_df) <- UniquePops
  window_df[is.na(window_df)] <- 0
  for(seq in rownames(window_df)){
    for(seq2 in rownames(window_df)){
      if(seq != seq2){
        window_df[seq,seq2] <- edges_for_bimodal_test[edges_for_bimodal_test$Pop1 == seq & edges_for_bimodal_test$Pop2 == seq2 | edges_for_bimodal_test$Pop1 == seq2 & edges_for_bimodal_test$Pop2 == seq,c(3)]
      }
    }
  }
  # min_mean_piw <- min(colMeans(WG_off[WG_off$chromosome_midMb == i,grep("piAdj", colnames(WG_off))], na.rm = T))
  # window_df <- window_df - min_mean_piw
  all_DMs_dXY[[i]] <- window_df
  all_hclusts_dXY[[i]] <- list(hclust(as.dist(window_df), method = "average"))
  all_agnes_dXY[[i]] <- list(agnes(as.dist(window_df), method = "average"))
}
# remove trees with 0 height
all_hclusts_dXY <- all_hclusts_dXY[unlist(lapply(all_hclusts_dXY, function(x) max(x[[1]]$height) > 0))]
all_agnes_dXY <- all_agnes_dXY[unlist(lapply(all_agnes_dXY, function(x) max(x[[1]]$height) > 0))]

# get table of summary stats for all trees
cluster_summaries <- getClusterSummaries(long_dXY, all_agnes_dXY)

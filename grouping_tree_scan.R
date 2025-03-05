###################
#### Functions ####
###################

# Run clustering
runClustering <- function(agList, dXY, seeds, cutoff) {
  corDF <- data.frame(matrix(ncol = nrow(na.omit(dXY)), 
                             nrow = nrow(na.omit(dXY))))
  rownames(corDF) <- rownames(na.omit(dXY))
  colnames(corDF) <- rownames(na.omit(dXY))
  bootList <- list()
  currBin <- c()
  allBinned <- c()
  currClust <- c()
  clustList <- list()
  for(j in seeds) {
    if(!is.null(ncol(corDF))) {
      if (!j %in% allBinned){
        corDF[j,] <- lapply(agList, function(x) cor(cophenetic(agList[[j]][[1]]), 
                                                    cophenetic(x[[1]])))
        currBin <- c(j, colnames(corDF[j,])[corDF[j,] > cutoff])
        currClust <- dXY[rownames(dXY) %in% currBin,]
        clustList <- append(clustList, list(currClust))
        agList <- agList[!names(agList) %in% currBin]
        corDF <- corDF[!rownames(corDF) %in% currBin,!colnames(corDF) %in% currBin]
        allBinned <- c(allBinned, currBin)
      }
    }
  }
  return(clustList)
}

# before proceeding, make sure you've loaded generate_dXY_trees.R

############################################################
#### Grouping tree scan based on cophenetic correlation ####
############################################################

# generate a list of seed regions for the cluster analysis
clustering_seeds <- sample(rownames(na.omit(long_dXY)))

# run the clustering algorithm, to generate forests based on cophenetic correlation coefficient
dXY_forests <- runClustering(all_agnes_dXY, long_dXY, clustering_seeds, 0.5)
dXY_filtered_forests <- dXY_forests[as.logical(lapply(dXY_forests, function(x) !is.null(nrow(x))))]

# generate summary statistics for all forests
dXY_forests_size <- unlist(lapply(dXY_filtered_forests, function(x) length(rownames(x))))
dXY_forests_height <- unlist(lapply(dXY_filtered_forests, function(x) mean(cluster_summaries[rownames(x),"tree_height"], na.rm = T)))
dXY_forests_longest_RB <- unlist(lapply(dXY_filtered_forests, function(x) mean(cluster_summaries[rownames(x),"longest_root_branch_length"], na.rm = T)))
dXY_forests_shortest_RB <- unlist(lapply(dXY_filtered_forests, function(x) mean(cluster_summaries[rownames(x),"shortest_root_branch_length"], na.rm = T)))

dXY_forests_stats <- cbind.data.frame(dXY_forests_size, dXY_forests_height, dXY_forests_longest_RB, dXY_forests_shortest_RB)

# generate mean trees for all forests
dXY_mean_forests <- list()
iter_names <- c()
for (i in seq(1, length(dXY_filtered_forests))) {
  iter_names <- c(iter_names, i)
  temp_dXY <- cbind.data.frame(all_dXY_edges[,c(1,2)], colMeans(dXY_filtered_forests[[i]]))
  dXY_df <- data.frame(matrix(nrow = length(UniquePops), ncol = length(UniquePops)))
  rownames(dXY_df) <- UniquePops
  colnames(dXY_df) <- UniquePops
  dXY_df[is.na(dXY_df)] <- 0
  for(pop in rownames(dXY_df)){
    for(pop2 in rownames(dXY_df)){
      if(pop != pop2){
        dXY_df[pop,pop2] <- temp_dXY[temp_dXY$Pop1 == pop & temp_dXY$Pop2 == pop2 | temp_dXY$Pop1 == pop2 & temp_dXY$Pop2 == pop,c(3)]
      }
    }
  }
  dend_dXY <- as.dendrogram(hclust(as.dist(dXY_df)))
  dXY_mean_forests[[i]] <- list(dend_dXY)
}

# forest summary plot - forest height against forest mean Shortest Root Branch (SRB)
# png(file="size_SFB_plot_V4.png", width = 1200, height = 1200, res = 200)
par(mfrow=c(1,1))
par(mar = c(5.1, 5.1, 0, 2.1))
plot(dXY_forests_stats$dXY_forests_shortest_RB, dXY_forests_stats$dXY_forests_size, log = "y", pch = 19,
     ylab = "Forest size", xlab = "Mean SRB", cex.lab = 2)
text(dXY_forests_stats$dXY_forests_shortest_RB[order(dXY_forests_stats$dXY_forests_shortest_RB, decreasing = T)][seq(1,3)], 
     dXY_forests_stats$dXY_forests_size[order(dXY_forests_stats$dXY_forests_shortest_RB, decreasing = T)][seq(1,3)], 
     adj = c(0.5,-1), labels = c("(I)", "(II)", "(III)"), col = "red")
text(dXY_forests_stats$dXY_forests_shortest_RB[order(dXY_forests_stats$dXY_forests_size, decreasing = T)][1], 
     dXY_forests_stats$dXY_forests_size[order(dXY_forests_stats$dXY_forests_size, decreasing = T)][1], 
     adj = c(0.5,--1.5), labels = c("(A)"), col = "blue")
text(dXY_forests_stats$dXY_forests_shortest_RB[order(dXY_forests_stats$dXY_forests_size, decreasing = T)][seq(2,3)], 
     dXY_forests_stats$dXY_forests_size[order(dXY_forests_stats$dXY_forests_size, decreasing = T)][seq(2,3)], 
     adj = c(0.5,-1), labels = c("(B)", "(C)"), col = "blue")
# dev.off()

############################################################
#### Grouping tree scan based on root division analysis ####
############################################################

# specify populations belonging to one clade of the tree topology you wish to scan for
clade_to_match <- pop_tab[pop_tab$pop_species == "ps","pop_names"]

# for all subgenomic trees, use the cutree function to split into the outermost two clades
# check if either of the two clades matches clade_to_match
rda_1 <- unlist(lapply(all_hclusts_dXY, function(x) all(clade_to_match %in% names(cutree(x[[1]], 2))[as.logical(cutree(x[[1]], 2) == 1)]) & sum(as.logical(cutree(x[[1]], 2) == 1)) == length(clade_to_match)))
rda_2 <- unlist(lapply(all_hclusts_dXY, function(x) all(clade_to_match %in% names(cutree(x[[1]], 2))[as.logical(cutree(x[[1]], 2) == 2)]) & sum(as.logical(cutree(x[[1]], 2) == 1)) == length(clade_to_match)))

# all subgenomic trees matching your specific root division
matching_trees <- c(names(which(rda_1)), names(which(rda_2)))

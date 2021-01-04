library(dplyr)
library(spaa)
library(reshape2)


# Get list of *.results files
files <- list.files()[grepl(".results",list.files())]
# Read them all into a list of data.frames()
results <- lapply(files,read.table,col.names=c("genome","gene","allele"))

# Create empty data.frame, fill it, convert to "distance matrix"
shared_mat <- d_mat <- diag(0,length(results))
row.names(shared_mat) <- row.names(d_mat) <- abbreviate(unlist(lapply(results,FUN=function(x){as.character(x$genome[1])})), minlength=15)
colnames(shared_mat) <- colnames(d_mat) <- abbreviate(unlist(lapply(results,FUN=function(x){as.character(x$genome[1])})), minlength=15)

for( i in 1:(length(results)-1) ) {
  for( j in (i+1):length(results) ) {
    tmp <- left_join(results[[i]],results[[j]],by="gene") %>%
      mutate(diff=allele.x!=allele.y,
             shared=!(is.na(allele.x)|is.na(allele.y)))
    diff <- sum(tmp$diff,na.rm=T)
    shared <- sum(tmp$shared,na.rm=T)

        # # Straight number of differing alleles
    d_mat[i,j] <- d_mat[j,i] <- diff
    shared_mat[i,j] <- shared_mat[j,i] <- shared
  }
}

# Transform the distance matrix into a pairwise comparison and apply a maximum allelic difference cut-off (default = 20)
d_pair <- subset(melt(d_mat), value<=20)

# Output the pairwise file with cut-off
write.table(d_pair, file="/PATH/TO/OUTPUT/FILE/d_pair.txt", sep="\t")

# Output the matrix of shared genes (out of 1,748) between each pair of genomes
write.table(shared_mat, file="PATH/TO/OUTPUT/FILE/shared_mat.txt", sep="\t")

# Create a matrix of proportion of shared alleles between each pair of genomes
p_mat <- d_mat/shared_mat

# Output the matrix of proportion of shared alleles between each pair of genomes
write.table(p_mat, file="PATH/TO/OUTPUT/FILE/p_mat.txt", sep="\t")

# Output the matrix of absolute distances (number of allelic differences) between each pair of genomes
write.table(d_mat, file="PATH/TO/OUTPUT/FILE/d_mat.txt", sep="\t")




# Store distance matrix as dataframe d_df
raw_d_df <- as.matrix(d_mat) %>%
  data.frame()
d_df <- as.matrix(p_mat) %>%
  data.frame()


# dendrogram from hclust object using complete linkage (default)
p_mat_dist <- as.dist(p_mat)
hc <- hclust(p_mat_dist,method="complete")
plot(hc)

# agglomerative clustering from cluster::agnes also using complete linkage (default)
library(cluster)
ag <- agnes(d_mat,diss=TRUE,method="complete")
summary(ag)
plot(ag)


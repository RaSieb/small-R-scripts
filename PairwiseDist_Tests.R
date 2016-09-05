#!/usr/bin/Rscript
# This script is based on the packages "adegenet" and "ape" and 
# provides some usefull tools to describe groups (read from a .csv) 
# in a phylogenetic (NJ-)tree created from from a fasta alignment of SNPs.

require(adegenet)
# load data
snp_matr <- fasta2genlight("./test_fasta.fasta", saveNbAlleles=T, n.cores = 8)
epi_data <- read.csv("./test_grps.csv", sep = ";", header = T)

# Assign groups to sequence file
pop(snp_matr) <- epi_data[match(indNames(snp_matr) , epi_data$name), "group"]
# Overview group sizes
table(pop(snp_matr))

#Distance matrix of the whole dataset
D <- dist(as.matrix(snp_matr), method = "manhattan")

# plot a NJ-tree
require(ape)
quartz(width = 6, height = 6) # Open a plot window (skip and use default plot window if not working)
plot(nj(D), type="unrooted", tip.color=fac2col(pop(snp_matr), col.pal=funky, na.col="transparent", seed=NULL),.6, cex=0.5)

## Plot the distances within groups
D_grp <- seppop(snp_matr)
quartz(width = 8, height = 8) # Open a plot window (skip and use default plot window if not working)
par(mfrow=c(ceiling(sqrt(length(levels(pop(snp_matr))))), ceiling(sqrt(length(levels(pop(snp_matr)))))))
dm <- lapply(seq_along(D_grp), function(x) {
  d <- dist(as.matrix(D_grp[[x]]), method = "manhattan")
  hist(d, main=paste(names(D_grp)[x]))
  abline(v=mean(d), col="firebrick", lwd=3)
  as.vector(d)
}) 

## Is there a difference in the distances within groups? (p-values of pairwise t-tests for all groups)
m <- matrix(rep(0, length(dm)^2), nrow = length(dm))
dimnames(m) <- list(names(D_grp), names(D_grp))
for (i in 1:(length(dm))) {
  for (j in i:length(dm)){
    m[i,j] <- m[j,i] <- t.test(dm[[j]], dm[[i]])$p.value
  }
};m

### Distance between group members and non-members
# Which group do you want to look at?
test_group = "3"

gr <- ifelse(epi_data$group==test_group, "in", "out")[match(attr(D, which = "Labels"), epi_data$name)]
cl2rest <- as.matrix(D)[gr=="out", gr=="in"]
clWithin <- as.matrix(D)[gr=="in", gr=="in"]

# Plot histogram of distances and means
quartz(width = 6, height = 4) # Open a plot window (skip and use default plot window if not working)
d <- density(clWithin)
e <- density(cl2rest)
plot(NULL, type="n", main="Distribution of genetic distances", ylab="proportion", xlab="Distance"
     , xlim=c(min(d$x, e$x),max(d$x, e$x))
     , ylim=c(min(d$y, e$y),max(d$y, e$y)+0.2))
polygon(e, col="#ff606080", border=NA); rug(cl2rest, side = 1, col="#ff606080", lwd = 2, ticksize = 0.1)
abline(v = mean(cl2rest), lwd = 2, col = "#ff6060")
polygon(d, col="#6060ff80", border=NA); rug(clWithin, side = 3, col="#6060ff80", lwd = 2, ticksize = 0.1)
abline(v = mean(clWithin), lwd = 2, col = "#6060ff")

# T-test for difference in means (be carefull, the test of means over the two matrices can be missleading)
t.test(clWithin, cl2rest)
# Minimum distance from group to outgroup vs. maximum distance between group members
min(cl2rest)
max(clWithin)

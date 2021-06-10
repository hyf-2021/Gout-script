argv <- commandArgs(T)
if(length(argv) != 3) { 
	stop ("Rscript argv[0] <abd.tab> <distance_method> <outfile>
	abd.tab: Abundance matrix of gene/genus/species etc. Rows represent features (species for example), columns represent observations (samples).
	dist_method: Distance method to calculate sample-sample distance, could be any distance method used in 'vegdist'. The accepted distances are:  manhattan, euclidean, canberra, bray, kulczynski, jaccard, gower, altGower, morisita, horn, mountford, raup, binomial, chao, cao or mahalanobis.
	outfile: ouput file")
}

if (grepl("\\.gz$",argv[1]))
{
	dat <- read.table(gzfile(argv[1]), sep = "\t", header = T, row.names = 1,check.names=F)
} else {
	dat <- read.table(argv[1], sep = "\t", header = T, row.names = 1,check.names=F)
}

dist_method <- argv[2]

outfile <- argv[3]
dist_m <- matrix(0, ncol(dat), ncol(dat))
rownames(dist_m) <- colnames(dat)
colnames(dist_m) <- colnames(dat)

library("vegan")
for(i in 1:(ncol(dat)-1)) {
	for(j in (i + 1):ncol(dat)) {
		dist_m[i, j] <- vegdist(t(dat[, c(i, j)]), method = dist_method,binary=TRUE)
		dist_m[j, i] <- vegdist(t(dat[, c(i, j)]), method = dist_method,binary=TRUE)
	}
}

write.table(dist_m, outfile, quote = F, sep = "\t", row.names = T, col.names = NA)

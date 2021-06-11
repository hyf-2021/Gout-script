argv <- commandArgs(T)
if(length(argv) != 4){stop("Rscripts pam.r [input.prof] [sample.config] [MAXIMUM] [OUTDIR]
			input.prof: profile table, features in row while samples in column, row with identifier of \"Unclassified\" will be removed for JSD calculation
			sample.config: group information of samples, column1 is sample ID and column2 is group info.
			MAXIMUM: maximun number of clusters to be evaluated before to determine the best one
			OUTDIR: directory where to store the output results
#LOW: remove the low abundance or not.Y/N
			")}

ShareTable <- argv[1]
Config <- argv[2]
ROUNDS <- as.numeric(argv[3])
dir.create(argv[4],recursive = TRUE,mode = "0755")
OUTDIR <- argv[4]
genus <- read.table(ShareTable, head=T, check.names=F, row.names=1, dec=".")
if(ROUNDS > ncol(genus)) ROUNDS <- ncol(genus) - 1
col_name <- colnames(genus)
group <- read.table(Config, head=T)
rownames(group) <- group[, 1]
group <- group[col_name, ]  ## sorted group according to the sample order of genus profile

colnames(group) <- c("SampleID","Group")

#----- For plot -----
pchs <- c(21, 24, 15, 21, 23, 24, 25, 0, 7:12)  ## indicates samples group information, support at most 14 subgroups
group_pch <- pchs[as.factor(group[, 2])]
legend_pch <- pchs[1:length(group[, 2])]
colors <- c("#4DAF4A","#E41A1C","#377EB8","purple","orange", "darkorchid", "cornflowerblue", "brown1", "deeppink", "cyan","black") ## indicates clusters, support at most 11 clusters F4A"

#----- Preprocessing -----
cat("\n***Preprocessing:\n")
if(sum(colSums(genus)) > (ncol(genus) + 1)) {
	cat("Normalizing data to the sample's total reads count\n")
		genus <- sweep(genus, 2, colSums(genus), "/")
} else {}
cat("Removing \'Unclassified\' feature without further normalization\n***\n")
genus <- genus[rownames(genus)!="Unclassified", ] ## remove unclassified row
genus <- genus[rownames(genus)!="unknown", ]
genus <- genus[rownames(genus)!="Unknown", ]

#----- delete for JSD calculation -----
genus <- as.data.frame(genus)

############################ Subfunctions ####################################
JSD <- function(InMat, pseudocount = 0.0000000000001){
	kld <- function(x, y){ sum(x * log(x/y)) }
	jsd <- function(x,y){ sqrt(0.5*kld(x, (x+y)/2) + 0.5*kld(y, (x+y)/2)) }
	ncol1 <- length(colnames(InMat))
	colname <- colnames(InMat)
	resultMatrix <- matrix(0, ncol1, ncol1)
	InMat <- apply(InMat, 1:2, function(x) ifelse (x==0, pseudocount, x))
	for(i in 1:ncol1){
		for(j in 1:ncol1){
			resultMatrix[i, j] <- jsd(as.vector(InMat[, i]), as.vector(InMat[, j]))
			}
	}
	colname -> colnames(resultMatrix) -> rownames(resultMatrix)
	resultMatrix <- as.dist(resultMatrix)
	attr(resultMatrix, "method") <- "dist"
	return(resultMatrix)
}

#----- denoise only for BCA -----
Noise.remove <- function(dataframe, percent=0.01, top=NULL){ ## cutoff 0.01%
	Matrix <- dataframe
	big.ones <- rowSums(Matrix)*100 / (sum(rowSums(Matrix))) > percent
	Matrix.1 <- Matrix[big.ones, ]
	return(Matrix.1)
}

#----- PAM ----- 
pam.clustering <- function(x, k) { # x is a distance matrix and k the number of clusters
	library("cluster",lib.loc="/hwfssz4/BC_PUB/Software/03.Soft_ALL/R-3.5.1/library/")
	cluster <- as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
	return(cluster)
}

############################ Subfunctions ####################################
genus.dist <- JSD(genus)  ## whithout denoising 
write.table(as.matrix(genus.dist), paste(OUTDIR, "/JSD_distance.txt", sep = ""), quote = F, sep = "\t")

library("clusterSim",lib.loc="/hwfssz4/BC_PUB/Software/03.Soft_ALL/R-3.5.1/library/")
num.cluster <- NULL
max.ch <- 0
for(k in 1:ROUNDS){
	if(k==1){ 
		num.cluster[k]=0 
	} else{
		genus.cluster_temp <- pam.clustering(genus.dist, k)
		num.cluster[k] <- index.G1(t(genus), genus.cluster_temp, d=genus.dist, centrotypes="medoids")
		if (num.cluster[k] > max.ch)
		{
			max.ch <- num.cluster[k]
			CLUSTERS <- k
		}
	}
}

cat("Line 99 CLUSTER", CLUSTERS, "\n", sep=" ")
genus.cluster <- pam.clustering(genus.dist, k=CLUSTERS)  
cluster.result <- cbind(col_name,genus.cluster)
write.table(cluster.result, paste(OUTDIR, "/cluster.result.xls", sep=""), row.names = F, col.names = F, sep = "\t", quote = F)

#----- to make visualization -----
library("ggplot2",lib.loc="/hwfssz4/BC_PUB/Software/03.Soft_ALL/R-3.5.1/library/")
pdf(paste(OUTDIR,"/CH.Index.plot.pdf",sep=""),width=4,height=3) 	
chplotdata <- data.frame(vec=num.cluster)
ggplot(chplotdata,aes(x=1:ROUNDS,y=vec)) + 
	geom_line() +
	labs(x="K clusters",y="CH Index") +
	theme (axis.text.x = element_text(size=7),
			axis.text.y = element_text(size=7),
			plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "in"),
			panel.background = element_blank(),
			panel.border = element_rect(fill=NA,colour="black")
	) +
	scale_x_continuous(breaks = round(seq(1, ROUNDS, by = 1),1)) +
	annotate("rect",xmin=CLUSTERS-0.3,xmax=CLUSTERS+0.3,ymin=-Inf,ymax=Inf,alpha=.4)
dev.off()


#----- for driver species plot -----
library("ade4",lib.loc="/hwfssz4/BC_PUB/Software/03.Soft_ALL/R-3.5.1/library/")
library("reshape2",lib.loc="/hwfssz4/BC_PUB/Software/03.Soft_ALL/R-3.5.1/library/")

cluster.matrix <- as.matrix(paste("Enterotype",genus.cluster,sep=""))
rownames(cluster.matrix) <- col_name

#----- plot rate of Enterotype in each group -----
EnterotypeRate <- cbind(group,Type=(cluster.matrix[match(rownames(group),rownames(cluster.matrix))]))
rownames(EnterotypeRate) <- NULL
EnterotypeRate <- table(EnterotypeRate[,2:3])
EnterotypeRate <- melt(sweep(EnterotypeRate, 1, rowSums(EnterotypeRate), FUN="/"))

pdf(paste(OUTDIR,"/Enterotype_rate.pdf",sep=""),width=3,height=2.5)
ggplot(EnterotypeRate,aes(x=Group,y=value,fill=Type)) + 
	geom_bar(stat = "identity",position="fill",colour="white") +
	labs(x="",y="",fill="") +
	scale_fill_manual(values = colors) +
	coord_equal(1/0.25) +
	scale_y_continuous(expand = c(0,0)) +
	theme(axis.text.y = element_blank(),
		axis.text.x = element_text(size=8),
		axis.ticks.y = element_blank(),
		axis.ticks.x = element_blank(),
		legend.text = element_text(size=7),
		legend.position="top",
		legend.key.width = unit(0.15, "in"),
		legend.key.height = unit(0.15, "in"),
		legend.margin=margin(0,0,-8,0),
		panel.background = element_blank(),
		plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "in")
	)
dev.off()

#----- Run PCA-----
obs.pca <- dudi.pca(data.frame(t(sqrt(genus))), scale=F,scannf=F, nf=10)
	
#----- ouput driver items, top 20 -----
drivers <- head(obs.pca$c1[order(abs(obs.pca$c1[, 1]) + abs(obs.pca$c1[, 2]), decreasing=T),], n=10)
write.table(drivers, paste(OUTDIR, "/PCA.drivers.txt", sep = ""), sep = "\t", quote = F, row.names = T, col.names = F)

#----- boxplot of drivers -----
drivers.abd <- genus[rownames(drivers), ]
drivers.abd <- log10(drivers.abd + 1e-12) ## introduce a small number for log
	
data.temp <- as.data.frame(cbind(Type=(cluster.matrix[match(rownames(t(drivers.abd)),rownames(cluster.matrix))]),t(drivers.abd)))
rownames(data.temp) <- NULL
write.table(data.temp,paste(OUTDIR, "/PCA.drivers_temp.txt", sep = ""), sep = "\t", quote = F, row.names = F, col.names = T)

data.temp <- read.table(paste(OUTDIR, "/PCA.drivers_temp.txt", sep = ""),head=T)
data.melt <- melt(data.temp,id.vars="Type")
data.melt$Type <- factor(data.melt$Type)
	
pdf(paste(OUTDIR, "/PCA.drivers.boxplot.pdf", sep=""),width=4,height=3)
ggplot(data.melt,aes(x=variable,y=value)) +
	geom_boxplot(aes(fill=factor(Type)),outlier.size = 0.5,position=position_dodge(0.9)) +
	geom_point(aes(fill=factor(Type)),size=0.1,position = position_jitterdodge())+
	labs(x="",y=expression(Relative~abundance~(log["10"])),fill="") +
	scale_fill_manual(values = colors) +
	theme(
		axis.text.x = element_text(colour="black",size=7,face="italic", hjust=1, angle=45),
		axis.line = element_line(color="black"),
		legend.position=c(0,0),
		legend.justification=c(0,0),
		legend.key = element_blank(),
		legend.text = element_text(size=9),
		legend.key.width = unit(0.15, "in"),
		legend.key.height = unit(0.2, "in"),
		legend.background = element_blank(),
		panel.background = element_blank(),
		plot.margin = unit(c(0.2, 0.2, 0.1, 0.2), "in")
	)
dev.off()

#----- Output a PCA plot when CLUSTERS == 2 -----
cat("\n***\nThe number of cluster is", CLUSTERS, "which is less than 3, thus the BCA cannot work, instead PCA is performed.\n***\n")

pca.pc1 <- sprintf("%0.2f",100*obs.pca$eig[1]/sum(obs.pca$eig))
pca.pc2 <- sprintf("%0.2f",100*obs.pca$eig[2]/sum(obs.pca$eig))
			
pca_figure <- paste(OUTDIR, "/PCA.cluster.plot.pdf", sep = "")
pdf(pca_figure,width=3.2, height=3.4)
plot(obs.pca$li[,1:2],type="n",frame.plot=T,xlab=paste("PC1"," (",pca.pc1,"%",")",sep=""),ylab=paste("PC2"," (",pca.pc2,"%",")",sep=""),tcl=-0.3,mgp=c(1.5,0.2,0),las=2,xlim=c(1.2*min(obs.pca$li[,1]),1.2*max(obs.pca$li[,1])),ylim=c(1.2*min(obs.pca$li[,2]),1.2*max(obs.pca$li[,2])),cex.axis=0.8,cex.lab=0.8,main="Principle component analysis",cex.main=1.2,xaxt="n",yaxt="n")
axis(side=1,tcl=-0.2,lwd.ticks=1.4,mgp=c(1.5,0.2,0),cex.axis=0.8)
axis(side=2,tcl=-0.2,lwd.ticks=1.4,mgp=c(1.5,0.2,0),cex.axis=0.8)
box(lwd=1.4)
abline(h=0,v=0,col=rgb(0,0,0,0.4),lwd=1.4)
s.class(obs.pca$li[,1:2], fac = as.factor(genus.cluster), grid = F, col =colors, pch=group_pch, csub=0.7,cpoint=0.6,clabel=0,add.plot=T,cstar=0,cellipse=1.5,addaxes=F,axesell=F)
for(i in 1:CLUSTERS){
		text(0.4*drivers[i,1],0.4*drivers[i,2],labels=rownames(drivers)[i],cex=0.7)
	}
legend("bottomleft",pch=legend_pch, col="black", legend=levels(as.factor(group[, 2])),cex=0.7, pt.cex=0.7,xpd=T,bty="n", y.intersp=1)
legend("topleft",pch=16,col=colors,legend=levels(as.factor(cluster.matrix)),cex=0.7, pt.cex=0.7,xpd=T,bty="n", y.intersp=1)
dev.off()

#----- output samples' coordinates for further visualization -----
colnames(obs.pca$li) <- paste("PC", 1:ncol(obs.pca$li), sep = "")
write.table(obs.pca$li, paste(OUTDIR, "/PCA.coordinates.txt", sep = ""), quote = F, sep = "\t")

#----- ouput eigen value of PCA -----
eig1 <- cbind(obs.pca$eig, obs.pca$eig/sum(obs.pca$eig))
colnames(eig1) <- c("eigenValue", "%")
rownames(eig1) <- paste("PC", 1:nrow(eig1), sep="")
write.table(eig1, paste(OUTDIR, "/PCA.eigenValue.txt", sep = ""), quote = F, sep = "\t")

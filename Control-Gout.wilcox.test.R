
infile <- file("Control-Gout.Genus.profile.xls","r")
outfile <- file("Control-Gout.wilcox.test.xls","w")

line <- readLines(infile, n=1)
samp <- as.vector(unlist(strsplit(line, split="	")))

phen <- read.table("Control-Gout.group.info",head=T)
phen_samp <- cbind("samp"=as.vector(phen[,1]), "state"=as.vector(phen[,2]))
groupname <- as.matrix(sort(unique(phen$groupname)))
title <- "ID"
for (i in 1:nrow(groupname))
{
	mean <- paste("mean(",groupname[i,1],")",sep="")
	sd <- paste("sd(",groupname[i,1],")",sep="")
	rank <- paste("mean-rank(",groupname[i,1],")",sep="")
	occ <- paste("occ-rate(",groupname[i,1],")",sep="")
	title <- paste(title,mean,sd,rank,occ,sep="	")
}
outmatrixname <- paste(title,"enriched","pvalue",sep="\t")
writeLines(outmatrixname, con=outfile, sep="\n")

while(length(line <- readLines(infile, n=1))) {
	line <- as.vector(unlist(strsplit(line, split="	")))
	tmp <- cbind(samp=samp[-1], abund=line[-1])
	dat <- merge(phen_samp, tmp, by="samp")
	dat$abund <- as.numeric(as.vector((dat$abund)))
	dat <- dat[complete.cases(dat),]
	mean.sd <- as.matrix(aggregate(dat$abund,by=list(dat$state),FUN=function(x)c(mean=sprintf("%0.9f",mean(x)),sd=sprintf("%0.9f",sd(x)))))
	testresult <- wilcox.test(dat$abund ~ dat$state,alternative ="two.sided")
	Rank <- rank(dat$abund)
	Rank_mean <- as.matrix(aggregate(Rank,by=list(dat$state),FUN=function(x)c(mean=mean(x))))
	Occ <- as.matrix(aggregate(dat$abund,by=list(dat$state),FUN=function(x)sum(x != 0)/length(x)))
	enriched <- as.character(Rank_mean[which.max(Rank_mean[,2]),1])
	output <- paste(line[1],mean.sd[1,2],mean.sd[1,3],format(as.numeric(Rank_mean[1,2]),digits=0),format(as.numeric(Occ[1,2]),digits=3),mean.sd[2,2],mean.sd[2,3],format(as.numeric(Rank_mean[2,2]),digits=0),format(as.numeric(Occ[2,2]),digits=3),enriched,testresult$p.value,sep="\t")
	writeLines(output, con=outfile, sep="\n")
}
close(infile)
close(outfile)

data <- read.table("Control-Gout.wilcox.test.xls",head=T,check.names=F,sep="	")
p_adjusted <- p.adjust(data$pvalue, method="BH")
data$qvalue <- p_adjusted
write.table(format(data,scientific=FALSE,digits=9),file="Control-Gout.wilcox.test.xls",row.names=F, col.names=T, quote=F, sep="	")

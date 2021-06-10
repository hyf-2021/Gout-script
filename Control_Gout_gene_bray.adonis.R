library("vegan")
beta<- read.table("Gene_beta_bray.xls",header=T,row.names=1,sep="	")

title <- c("Phenotype","Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)","PermanovaQvalue")
write.table(t(title),file="./Control_Gout_gene.adonis.result.xls",quote = F, sep = "	", row.names = F, col.names = F)
beta.names <- colnames(beta)
beta <- as.matrix(beta)
rowtitle <- factor(c("SamplingTime","Age","BMI","CRP","eGFR","ESR","SCr","SUA","Urea"))
outmatrix <- matrix(0, nrow = length(rowtitle), ncol = length(title))
for (i in 1:length(rowtitle))
{
	set.seed(999)
	info.file <- paste("./",rowtitle[i],".info_for_adonis.xls",sep="")
	data.info <-  read.table(info.file,header=T,row.names=1,sep="	")

	data.info.names <- rownames(data.info)
	intersection <- intersect(data.info.names,beta.names)
	beta_subsample <- beta[pmatch(intersection,beta.names),pmatch(intersection,beta.names)]

	run1 <- paste("permanova.result <- adonis(formula=as.dist(beta_subsample)~",rowtitle[i],",data=data.info,permutations=9999)",sep="")
	eval(parse(text = run1))
	result.table <- permanova.result$aov.tab
	
	out <- t(c(rownames(result.table)[1],result.table$Df[1],result.table$SumsOfSqs[1],result.table$MeanSqs[1],result.table$F.Model[1],result.table$R2[1],format(result.table$"Pr(>F)"[1],digits=9)))
	outmatrix[i, 1:length(out)] <- out
}
outmatrix[,ncol(outmatrix)] <- p.adjust(outmatrix[, ncol(outmatrix)-1],method="BH")

write.table(outmatrix,file="./Control_Gout_gene.adonis.result.xls",row.names=F,col.names=F,quote=F,sep ="	",append=T)

library('randomForest')
library('ggplot2')
library('reshape2')
library('pROC')
library('dplyr')
library('robustbase')
library('cowplot')

dat1 <- read.table('Gout_train_test_Control-Gout_Train_Sample_profile.xls', head=T,sep='	',row.names=1)
conf1 <- read.table('Gout_train_test_Control-Gout_Train_Sample_group.info',head=T,row.names=1,sep='	')

outcome = conf1$state
outcome <- sub('Control','0',outcome)
outcome <- sub('Gout','1',outcome)
outcome <- as.factor(outcome)
X <- as.data.frame(t(dat1))
X$outcome = outcome

#----- crossvalidation -----
pdf('Gout_train_test_Control-Gout_ROC.pdf',width=9,height=2)
#par(mfrow = c(1,3),mai=c(0.54,0.5,0.29,0.32),tcl=-0.3,lwd=1.4,mgp=c(2.3,0.5,0),las=1)

set.seed(999)
source('randomforest.crossvalidation_v2.2.r')
result <- replicate(5, rfcv1(X[,-ncol(X)], X$outcome, cv.fold=10,step=0.9), simplify=FALSE)
error.cv <- sapply(result, '[[', 'error.cv')
error.cv.cbm <- as.data.frame(cbind(rowMeans(error.cv), error.cv))
cutoff <- min(error.cv.cbm[,1])+sd(error.cv.cbm[,1])
cutoff.num <- nrow(error.cv.cbm[error.cv.cbm[,1]<=cutoff,])
optimal.set.feature.num <- as.numeric(rownames(error.cv.cbm[error.cv.cbm[,1]<=cutoff,])[cutoff.num]) # obtain marker number

##### plot figure 1 #####
error.cv.data <- cbind( error.cv,rowMeans(error.cv))
colnames(error.cv.data) <- c('E1','E2','E3','E4','E5','Mean') 
error.cv.data <- melt(error.cv.data)

number_ticks <- function(n) {function(limits) pretty(limits, n)} ## for axis tick number
P1 <- ggplot(error.cv.data,aes(x=Var1,y=value,color=Var2)) +
	geom_line() +
	geom_vline(xintercept = optimal.set.feature.num, color='red') +
	coord_trans(x="log2") +
	scale_x_continuous(breaks=c(1,2,5,10,20,50),labels=c(c(1,2,5,10,20,50))) +
	labs(x='Number of variables',y='CV Error') +
	scale_color_manual(values= c('gray','gray','gray','gray','gray','black')) +
	theme(axis.text.x = element_text(color='black',size=7),
		axis.text.y = element_text(color='black',size=7,angle=90,hjust=0.5),
		axis.title = element_text(color='black',size=7),
		axis.line = element_line(color='black'),
		axis.ticks = element_line(color='black'),
		plot.margin = unit(c(0.351, 0.05, 0.05, 0.051), 'in'),
		legend.position = 'none',
		panel.border = element_rect(colour = 'black', fill=NA),
		panel.grid = element_blank(),
		panel.background = element_blank()
	)

#----- pick out markers from the crossvalidation result -----
k = 1
b <- matrix(0,ncol=100,nrow=50)	## ncol is gene, genus or mlg number, nrow = 5*10
accuracy <- matrix(0,ncol=100,nrow=50) ## for plot feature importance
gini <- matrix(0,ncol=100,nrow=50)

for(i in 1:5)
{
	for(j in 1:10)
	{
		b[k,] <- result[[i]]$res[[j]]
		accuracy[k,] <- result[[i]]$acc[[j]]
		gini[k,] <- result[[i]]$gini[[j]]
		k = k+1
	}
}

mlg.list <- b[,1:10]
list <- c()
k = 1
for(i in 1:10)
{
	for(j in 1:50)
	{
		list[k] <- mlg.list[j,i]
		k = k+1
	}
}

mlg.sort <- as.matrix(table(list))
mlg.sort <- mlg.sort[rev(order(mlg.sort[,1])),]
pick <- as.numeric(names(head(mlg.sort,optimal.set.feature.num)))
tmp = X[,-ncol(X)]
mlg.pick <- colnames(tmp)[pick]
write.table(mlg.pick,'Gout_train_test_Control-Gout.cross_validation_pick.txt',sep='	',quote=F)


colnames(accuracy) <- rownames(dat1)
colnames(gini) <- rownames(dat1)

median_acc <- as.matrix(rowMedians(t(accuracy)))
median_gini <- as.matrix(rowMedians(t(gini)))

importance.data <- cbind(median_acc,median_gini)
colnames(importance.data) <- c('MeanDecreaseAccuracy','MeanDecreaseGini')
write.table(importance.data,'Gout_train_test_Control-Gout.feature.importance.xls',sep='	',quote=F,col.names=NA)

#----- predict result of train.set by markers -----
train1 <- X[,c(pick,100+1)] ##
set.seed(999)
train1 <- data.frame(train1)
train1.rf <- randomForest(outcome~., data = train1,importance = TRUE)
train1.pre <- predict(train1.rf,type='prob')
p.train <- train1.pre[,2]
combine <- as.data.frame(cbind(predict.value=as.matrix(p.train)[match(rownames(as.matrix(p.train)),rownames(as.matrix(conf1)))],as.matrix(conf1)))
write.table(combine,'Gout_train_test_Control-Gout.cross_validation.marker.predict.in.train.txt',sep='	',quote=F)

###### plot figure 2 #####
temp.data <- read.table('Gout_train_test_Control-Gout.cross_validation.marker.predict.in.train.txt',head=T)
temp.data$state <- factor(temp.data$state,levels=c('Control','Gout'))
colors <- c('#4DAF4A','#E41A1C')

P2 <- ggplot(data=temp.data,aes(x=state,y=predict.value)) +
	stat_boxplot(geom = 'errorbar',width=0.3,color=colors) +
	geom_boxplot(aes(color=state),lwd=0.5,outlier.shape=1,notch=T) + 
	geom_hline(yintercept=0.5,linetype='dashed',color='grey') +
	labs(x='',y='Probability of Disease',fill= '') +
	scale_color_manual(values= colors) +
	scale_y_continuous(limits = c(0, 1),breaks=seq(0, 1, 0.5)) + 
	theme(axis.text.x = element_text(color='black',size=7),
		axis.text.y = element_text(color='black',size=7,angle=90,hjust=0.5),
		axis.title = element_text(color='black',size=7),
		axis.line = element_line(color='black'),
		axis.ticks = element_line(color='black'),
		plot.margin = unit(c(0.351, 0.05, 0.05, 0.051), 'in'),
		legend.position = 'none',
		panel.border = element_rect(colour = 'black', fill=NA),
		panel.grid = element_blank(),
		panel.background = element_blank()
	)

#----- ROC in train -----
##### plot figure 3 #####
train.roc <- roc(outcome,p.train,percent=FALSE,partial.auc.correct=TRUE,ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,plot=F,auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE)
sens.ci <- as.data.frame(ci.se(train.roc, specificities=seq(0, 1, 0.05)))
sens.ci <- setNames(cbind(rownames(sens.ci),sens.ci,row.names=NULL),c('sp','se.low','se.median','se.high'))
write.table(sens.ci,'Gout_train_test_Control-Gout_train_ci_result.txt',sep='	',quote=F,col.names=T,row.names=F)
train.ci <- read.table('Gout_train_test_Control-Gout_train_ci_result.txt',head=T)
P3 <- ggroc(train.roc,color='red',legacy.axes=TRUE)+
	annotate(geom='text',x=0.2,y=0.15,label=paste('AUC =',round(train.roc$ci[2],4)),size=2.5,hjust=0) +
	annotate(geom='text',x=0.2,y=0.05,label=paste('95% CI:',round(train.roc$ci[1],4),'-',round(train.roc$ci[3],4)),size=2.5,hjust=0) +
	geom_abline(intercept = 0, slope = 1, color = 'gray') +
	geom_ribbon(data=train.ci,aes(x=1-sp,ymin=se.low,ymax=se.high,fill='red'),fill='red',alpha=0.2) +
	scale_x_continuous(breaks=seq(0, 1, 0.5)) +
	scale_y_continuous(breaks=seq(0, 1, 0.5)) +
	theme(axis.text.x = element_text(color='black',size=7),
		axis.text.y = element_text(color='black',size=7,angle=90,hjust=0.5),
		axis.title = element_text(color='black',size=7),
		axis.line = element_line(color='black'),
		axis.ticks = element_line(color='black'),
		legend.position = 'none',
		plot.margin = unit(c(0.351, 0.05, 0.05, 0.051), 'in'),
		panel.border = element_rect(colour = 'black', fill=NA),
		panel.grid = element_blank(),
		panel.background = element_blank()
	)

#----- predice resule of test.set by markers -----
dat2 <- read.table('Gout_train_test_Control-Gout_Test_Sample_profile.xls',header=T,sep='	',row.names=1)
conf2 <- read.table('Gout_train_test_Control-Gout_Test_Sample_group.info',header=T,sep='	',row.names=1)
dat2 <- data.frame(t(dat2))

set.seed(999)
p.test<-predict(train1.rf, dat2,type='prob')
pre.test <- as.data.frame(cbind(pre.value=as.matrix(p.test[,2])[match(rownames(as.matrix(p.test)),rownames(as.matrix(conf2)))],as.matrix(conf2)))
write.table(pre.test,'Gout_train_test_Control-Gout.cross_validation.marker.predict.in.test.txt',sep='	',quote=F)

##### plot figure 4 #####
pre.data <- read.table('Gout_train_test_Control-Gout.cross_validation.marker.predict.in.test.txt',head=T)
rownames(pre.data) <- NULL
pre.data$state <- factor(pre.data$state,levels=c('Control','Gout'))
colors <- c('#4DAF4A','#E41A1C')
pre.data <- pre.data[order(pre.data$pre.value),]

P4 <- ggplot(pre.data) +
	geom_point(aes(as.numeric(reorder(as.numeric(row.names(pre.data)),pre.value)),pre.value,color=state),size=0.4) +
	geom_hline(yintercept=0.5,linetype='dashed',color='grey') +
	labs(x='Samples',y='Probability of Disease',color='') +
	scale_color_manual(values= colors) +
	scale_y_continuous(limits = c(0, 1),breaks=seq(0,1,0.5)) +
	scale_x_continuous(breaks=number_ticks(5)) +
	theme(axis.text.x = element_text(color='black',size=7),
		axis.text.y = element_text(color='black',angle=90,hjust=0.5,size=7),
		axis.title = element_text(color='black',size=7),
		axis.line = element_line(color='black'),
		axis.ticks = element_line(color='black'),
		legend.position = c(1,0),
		legend.justification = c(1,0),
		legend.text = element_text(size=7),
		legend.key = element_blank(),
		legend.background = element_blank(),
		legend.key.size = unit(0.1,'in'),
		plot.margin = unit(c(0.351, 0.05, 0.05, 0.051), 'in'),
		panel.border = element_rect(colour = 'black', fill=NA),
		panel.grid = element_blank(),
		panel.background = element_blank()
	)

#----- test.ROC -----
outcome.test = conf2$state
outcome.test <- sub('Control','0',outcome.test)
outcome.test <- sub('Gout','1',outcome.test)
 
##### plot figure 5 #####
test.roc <- roc(outcome.test,p.test[,2],percent=FALSE,partial.auc.correct=TRUE,ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,plot=F,auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE)
sens.ci <- as.data.frame(ci.se(test.roc, specificities=seq(0, 1, 0.05)))
sens.ci <- setNames(cbind(rownames(sens.ci),sens.ci,row.names=NULL),c('sp','se.low','se.median','se.high'))
write.table(sens.ci,'Gout_train_test_Control-Gout_test_ci_result.txt',sep='	',quote=F,col.names=T,row.names=F)
test.ci <- read.table('Gout_train_test_Control-Gout_test_ci_result.txt',head=T)
P5 <- ggroc(test.roc,color='red',legacy.axes=TRUE)+
	annotate(geom='text',x=0.2,y=0.15,label=paste('AUC =',round(test.roc$ci[2],4)),size=2.5,hjust=0) +
	annotate(geom='text',x=0.2,y=0.05,label=paste('95% CI:',round(test.roc$ci[1],4),'-',round(test.roc$ci[3],4)),size=2.5,hjust=0) +
	geom_abline(intercept = 0, slope = 1, color = 'gray') +
	geom_ribbon(data=test.ci,aes(x=1-sp,ymin=se.low,ymax=se.high),fill='red',alpha=0.2) +
	scale_x_continuous(breaks=seq(0, 1, 0.5)) +
	scale_y_continuous(breaks=seq(0, 1, 0.5)) +
	theme(axis.text.x = element_text(color='black',size=7),
		axis.text.y = element_text(color='black',size=7,angle=90,hjust=0.5),
		axis.title = element_text(color='black',size=7),
		axis.line = element_line(color='black'),
		axis.ticks = element_line(color='black'),
		legend.position = 'none',
		plot.margin = unit(c(0.351, 0.05, 0.05, 0.051), 'in'),
		panel.border = element_rect(colour = 'black', fill=NA),
		panel.grid = element_blank(),
		panel.background = element_blank()
	)
plot_grid(P1,P2,P3,P4,P5,ncol=5,align='hv')

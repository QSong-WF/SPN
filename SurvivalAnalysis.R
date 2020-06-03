setwd("D:/Breast_Cancer/Datasets") # set path for running the following program
install.packages(“survival”)
library("survival")  
# install and invoke the “survival” package, which is used for survival analysis

gen<-"significant subnework markers(SPNs)" # input the SPNs identifed after the three tests
NEW<-read.csv(file="GeneExpression.csv",row.names=1)# the gene expression data of breast cancer patients
NEW2<-NEW[colnames(NEW)%in% training[,1]]  # the gene expression matrix for training patients (same for the test set)

number<-0
atra<-matrix(nrow=length(gen),ncol=ncol(NEW2)) # the matrix of activity scores based on SPNs in training patients (same for the test set)
for (i in 1:length(gen) )
{
	if (length(gen[[i]])>0)
	{
	number=number+1
	tempp<-NEW2[rownames(NEW2) %in% gen[[i]],]
	tem<-(tempp-mean(as.matrix(tempp)))/sd(as.numeric(as.matrix(tempp)))
	atra[number,]<-colSums(tem)/sqrt(nrow(tem))
	}
}
rownames(atra)<-paste("sub",1:length(gen),sep=" ")
colnames(atra)<-colnames(NEW2)

#### survival analysis using hclust classification in training patients

#hclust on training patients(same for the test set); patients are classifed to two subgroups.
pp<-hclust(dist(t(atra)))
qq<-2
pat<- cutree(pp,k=qq)
pat1<-pat[pat==1]
length(pat1)
pat2<-pat[pat==2]
length(pat2)

ttime<-dmfs[rownames(dmfs) %in% training[,1],2:3]    #DMFS in training patients (same for the test set)
p1<-cbind(ttime[rownames(ttime) %in% names(pat1),],1)
colnames(p1)<-c("time","event","type")
p2<-cbind(ttime[rownames(ttime) %in% names(pat2),],2)
colnames(p2)<-c("time","event","type")
surv<-rbind(p1,p2)
MY_DMFSdata<-surv

my.fit <- survfit(Surv(MY_DMFSdata[,1],MY_DMFSdata[,2])~ unlist(MY_DMFSdata[,3]), conf.type="none") # The function "survfit" creats survival curves from the survival object created by the formula "Surv".
survdiff(Surv(MY_DMFSdata[,1],MY_DMFSdata[,2])~ unlist(MY_DMFSdata[,3])) # The function "survdiff" tests if there is a difference between the two survival curves using the log-rank test.

# plot the survival curves 
X11()
par(lwd=4)
plot(my.fit,xlab="Time",col=c(2,4,3),ylab="Probability of DMFS",pch ="s",font=4,font.lab=2,cex.lab=1.3,cex.axis=1.5,lwd=5)

foo<- c("(Ph) high","(Ph) intermediate",
		"(Ph) low")
nn<-c(sum(pat==1),sum(pat==2),sum(pat==3))
foo1 <- nn
foo2 <- paste(foo,foo1,sep= "    n=")
legend(4.8,0.23, paste(foo2),col=c(2,4,3),pch=15, lty=1,cex=1,text.font=2,box.col="white", pt.bg="white")
legend(0.2,0.2,  box.col="white", paste("p= 0"),cex=1.5,text.font=2)

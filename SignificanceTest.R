setwd("D:/Breast_Cancer/Datasets") # set path for running the following program

GEN <- read.csv(file="Subnetwork.csv",check.names=F,header=F,sep=NULL,stringsAsFactors=F) # input the identified subnetworks initiated from each gene
SCORE<-read.csv(file="Score.csv",check.names=F,header=F,stringsAsFactors=F)  # input the corresponding score of each subnetwork

###for the 1-th test

all1<-matrix(nrow=nrow(GEN),ncol=100)
rownames(all1)<-paste("Score",1:nrow(GEN),sep="")
j<-0
for ( i in 1:100 )
{
a<-read.csv(file=paste("Test1/SCORE",i,".csv",sep=""),check.names=F,header=F,stringsAsFactors=F)  # input the first 100 random permution scores as null distribution
b<-a[match(rownames(all1),a[,1]),2]
if(length(b)>0){
j<-j+1
all1[,j]<- b
}}
colnames(all1)<-paste("trial",1:100)

foo1<-quantile(all1, probs = seq(0, 1, 0.05),na.rm=T) # select the significant subnetworks during the 1-th test from the null distribution
fo1<-as.numeric(foo1[names(foo1)=="95%"])
select_test1<-SCORE[SCORE[,2]>=fo1,1]

###for the 2-th test
all2<-matrix(nrow=nrow(GEN),ncol=100)
rownames(all2)<-paste("Score",1:nrow(GEN),sep="")
j<-0
for ( i in 1:100 )
{
a<-read.csv(file=paste("Test2/SCORE",i,".csv",sep=""),check.names=F,header=F,stringsAsFactors=F) # input the second 100 random permution scores as null distribution
b<-a[match(rownames(all2),a[,1]),2]
if(length(b)>0){
j<-j+1
all2[,j]<- b
}}
colnames(all2)<-paste("trial",1:100)

foo2<-quantile(all2, probs = seq(0, 1, 0.05),na.rm=T)   # select the significant subnetworks during the t-th test from the null distribution
fo2<-as.numeric(foo2[names(foo2)=="95%"])
select_test2<-SCORE[SCORE[,2]>=fo3,1]

###for the 3-th test

all3<-matrix(nrow=nrow(GEN),ncol=100)
rownames(all3)<-paste("Score",1:nrow(GEN),sep="")
j<-0
for ( i in 1:100 )
{
a<-read.csv(file=paste("Test3/SCORE",i,".csv",sep=""),check.names=F,header=F,stringsAsFactors=F) # input the third 100 random permution scores as null distribution
b<-a[match(rownames(all3),a[,1]),2]
if(length(b)>0){
j<-j+1
all1[,j]<- b
}}
colnames(all3)<-paste("trial",1:100)

foo3<-quantile(all3, probs = seq(0, 1, 0.05),na.rm=T) # select the significant subnetworks during the 3-th test from the null distribution
fo3<-as.numeric(foo3[names(foo3)=="95%"])
select_test3<-SCORE[SCORE[,2]>=fo3,1]

# the final significant subnetwork markers through the above three tests 
select3<-gsub("Score","Gene",select_test3)
select3subnet<-GEN[match(select3,GEN[,1]),]
str(select3subnet)
unig <- select3subnet  # the final significant subnetwork markers

HighGen<-vector("list")
n<-0
for ( i in 1:nrow(unig))
  {
          temm<-as.character(unig[i,])
          temp<-strsplit(temm,'\"')
          tem<-unlist(temp)
          GE<-tem[-which(tem==",")]
          #ge<-GE[-1]
          if ( length(GE)>0 )
            {
              n<-n+1
              HighGen[[n]]<-GE
            }
        }

# invoke the function for writing list to a file
fnlist <- function(x, fil){ z <- deparse(substitute(x))
                            cat(z, "\n", file=fil)
                            nams=names(x) 
                            for (i in 1:length(x) ){ cat(nams[i], "\t",  x[[i]], "\n", file=fil, append=TRUE) }}      
                            
                                                                               
fnlist(HighGen, "SignificantMarkers.txt") # save the SPN to a file "SignificantMarkers.txt"


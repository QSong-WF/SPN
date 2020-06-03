setwd("D:/Breast_Cancer/Datasets") # set path for running the following program
install.packages(“survival”)
library("survival")  
# install and invoke the “survival” package, which is used for survival analysis

FI<-read.delim(file="FIdatabase.txt",sep="",header=F,stringsAsFactors=F)
# the human functional protein interaction networks database
FII<-FI
FI1<-sapply(FII,as.character)
geee<-FI1[(FI1[,1]%in%rownames(ALL)) &(FI1[,2]%in%rownames(ALL)),]
gene<-geee[geee[,2]%in%rownames(ALL),]
# organize the FI database to be a matrix with two columns

ALL<-read.csv(file="GeneExpression.csv",row.names=1)
# the gene expression data of breast cancer patients
dmfs<-read.csv(file="DMFSdata.csv",sep=",",row.names=1,check.names=F)
# the distant metastasis-free survival time for breast cancer patients
high <-read.csv(file="P_high.csv",header=F,stringsAsFactors=F)
# P-high tertile of breast cancer patients 
training <- high[,1] # the name of P-high patients

NEW1<-ALL[,colnames(ALL)%in% training]
# the gene expression data of P-high patients
DMFS<-dmfs[,2:3]
DMF<-DMFS[rownames(DMFS)%in% training,]
# the DMFS of P-high patients

begin<-unique(c(gene[,1],gene[,2])
# the total genes for the process below
N<-length(begin)  # number of genes

NETWORK<-vector("list")  # identified network during searching process
Score<-vector("list")  # the score of corresponding identified network
ptm <- proc.time() # for calculating the elapsed time

for (i in 1:N)    # the loop of greedy search algorithm
  {
    cat ("i=",i,"\n")  # print the i-th gene that serves as initiate gene for searching
    AcceptN<-1   # AcceptN represents the number of networks
    why<-1   # “why” represent the number of loops when adding and deleting genes 
    gene1<-begin[i]   # the initiate gene
    range<-unique(c(gene[which(gene[,1] %in% gene1),2],gene[which(gene[,2]%in% gene1),1]))  # the neighbors of gene in the FI database
    N2 <-length(range)
    my1<-NEW2[rownames(NEW2) %in% gene1,]  # the gene expression matrix across the patients
    myy1<-(my1-mean(as.matrix(my1)))/sd(as.numeric(as.matrix(my1))) # normalization of the matrix
    myda1<-colSums(myy1)/sqrt(nrow(myy1))  # calculation of activity score
    mydata1<-cbind(myda1[match(rownames(DMF),names(myda1))],DMF) #combine the activity value with DMFS, and the recurrence event
    colnames(mydata1)<-c("value","time","event")  # define the colnames
    DATAcox1<-coxph(Surv(time,event) ~ value,data=as.data.frame(mydata1)) # coxph function invoked in the survival package, which is used for fitting proportional hazards regression model
    MAXcox1=as.numeric(-log((summary(DATAcox1))$logtest[3])) # –log(P-value); P-value is obtained by the log-rank test
    GEN<-vector("list")
    MAXc<-vector("list")
    GEN[[1]]<-gene1
    MAXc[[1]]<-MAXcox1
    V<-matrix(0,nrow=30,ncol=30) # the position of the neighbor of genes
    rownames(V)=paste("row",1:30,"")
    colnames(V)=paste("col",1:30,"")
    E<-matrix(0,nrow=1,ncol=60)
    colnames(E)=paste("col",1:60,"")
    V[1,1]<-1
    E[1,1]<-gene1
    while ( why<30 )    
      {
			mm<-0
        MM<-0
        MAXcox2<-vector("list")
			# The loop for adding genes into current network; when we add a gene, we try to delete one of genes in current network to ensure higher score
        for ( j in 1:N2 )
          {
            cat(" j=",j,"\n")
            gene2<-unique(c(range[[j]],GEN[[AcceptN]]))
            my2<-NEW2[rownames(NEW2) %in% gene2,]
            myy2<-(my2-mean(as.matrix(my2)))/sd(as.numeric(as.matrix(my2)))
            myda2<-colSums(myy2)/sqrt(nrow(myy2))
            mydata2<-cbind(myda2[match(rownames(DMF),names(myda2))],DMF)
            colnames(mydata2)<-c("value","time","event")
            DATAcox2<-coxph(Surv(time,event) ~ value,data=as.data.frame(mydata2))
            temp2<-as.numeric(summary(DATAcox2)$logtest[3])
            if ( temp2 !=0 )
              {
                MAXcox2[[j]]=-log(temp2)
              }else {
                MAXcox2[[j]]=0
              }
          }   
        MAXCOX2<-max(unlist(MAXcox2))
        Ee<-range[which(MAXcox2 %in% MAXCOX2)]
        v1<- which(GEN[[AcceptN]]%in% gene[which(gene[,1] %in% Ee),2])
        v2<- which(GEN[[AcceptN]]%in% gene[which(gene[,2]%in% Ee),1])
        gene2<-unique(c(Ee,GEN[[AcceptN]]))
        gene2<-gene2[!is.na(gene2)]
        if (  MAXCOX2 > MAXc[[AcceptN]] && (MAXCOX2- MAXc[[AcceptN]])/MAXc[[AcceptN]] >0.1 )
          {
            AcceptN <- AcceptN + 1
            GEN[[AcceptN]]<- gene2
            MAXc[[AcceptN]]<-MAXCOX2
            E[1,AcceptN]<-Ee[1]
            if(length(v1)!=0){
              V[v1,AcceptN]<-1
              V[AcceptN,v1]<-1
            }
            if(length(v2)!=0){
              V[AcceptN,v2]<-1
              V[v2,AcceptN]<-1
            }
            MM<-AcceptN
          }        
			#deletion process
        G <- length(GEN[[AcceptN]])
        delete <- 1
        g <- 1
        km <- vector("list")
        km[[g]] <- 0
        while ( delete<10 ){
          cat(" delete=",delete,"\n")
          if( G -1 >0 )
            {
              MAXcox3<-vector("list")
              for (k in 1:G)
                {
                gene3<- GEN[[AcceptN]][-k]
                pos<-c(1:length(gene3))[-k]
                if ( (length(pos)>1&&sum(colSums(V[pos,pos])!=0)==length(pos)) |
                    (length(pos)==1&& V[pos,pos]!=0) )
                  {
                    my3<-NEW2[rownames(NEW2) %in% gene3,]
                    myy3<-(my3-mean(as.matrix(my3)))/sd(as.numeric(as.matrix(my3)))
                    myda3<-colSums(myy3)/sqrt(nrow(myy3))
                    mydata3<-cbind(myda3[match(rownames(DMF),names(myda3))],DMF)
                    colnames(mydata3)<-c("value","time","event")
                    DATAcox3<-coxph(Surv(time,event) ~ value,data=as.data.frame(mydata3))
                    temp3<-as.numeric(summary(DATAcox3)$logtest[3])
                    if ( temp3 != 0 )
                      {
                        MAXcox3[[k]]=-log(temp3)
                      }else {
                        MAXcox3[[k]]=0
                      }
                    } else {
                      MAXcox3[[k]]= 0
                    }                  # for if condition
                cat(" k=",k)
              } 
              if ( sum(MAXcox3!=0)>0 ){
                foo <- max(unlist(MAXcox3))
                foopos <- which(MAXcox3 %in% foo)
                gen3 <- GEN[[AcceptN]][-foopos]
              if ( ( foo > MAXc[[AcceptN]] ) && ( (foo- MAXc[[AcceptN]])/MAXc[[AcceptN]] > 0.1 ) )
                {
                  g <- g+1
                  AcceptN <- AcceptN+1
                  GEN[[AcceptN]] <- gen3
                  MAXc[[AcceptN]] <- foo
                  V <- V[-foopos,-foopos]
                  E[1,foopos] <- 0
                  G <- length(GEN[[AcceptN]])
                  km[[g]] <- G
                  mm <- mm+1
                }
              }
              if ( (g==1)|(g>2 &&  km[[g]]==km[[g-1]] ) ) {
                break
              }
            } else {      
              break
            }
          delete <- delete+1
        } 
            range<-unique(c(gene[which(gene[,1] %in% GEN[[AcceptN]]),2],
                            gene[which(gene[,2]%in% GEN[[AcceptN]]),1]))
            N2 <-length(range)
            gene1<-GEN[[AcceptN]]
            why<-why+1
        if ( (MM==0)  && (mm ==0) )
          {
            break
          }
      }  
# through the adding process and deleting process, we get the network and corresponding score initiated from the start gene
    NETWORK[[i]]<- GEN[[length(GEN)]]
    Score[[i]] <- MAXc[[length(MAXc)]]
    write.table(t(NETWORK[[i]]), file ="Subnetwork.csv",
                row.names = paste("Network",27*(id-1)+i,sep=""),append = T, col.names = FALSE, sep =",")
# save the identified network 
    write.table(t(Score[[i]]), file ="Score.csv",
                row.names = paste("Score",27*(id-1)+i,sep=""),append = T, col.names = FALSE, sep =",")
# save the corresponding score of network
}
proc.time() – ptm
# show the elapsed time for running the above program




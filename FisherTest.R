setwd("D:/Breast_Cancer/Datasets") # set path for running the following program

subnetwork<- "identified subnetwork markers"
fileName="BP.txt" # read the file "BP.txt" 
#fileName="KEGG.txt"  # read the file "KEGG.txt" 
conn=file(fileName,open="r")
linn=readLines(conn)
close(conn)
temp <- grep("\t",linn)
result_bp<- lapply(linn[temp],function(x) unlist(strsplit(x,"\t")))
result <- result_bp # extract the biological process sets from the file; same for the KEGG pathway sets
name <- unlist(lapply(result,function(x) x[1])) # the name of each biological process; same for the KEGG pathway sets

net <- subnetwork  # input the identified network markers 
score <- matrix(nrow=length(net),ncol=length(result))
rownames(score) <-paste("net",seq(1:length(net)))
colnames(score) <- name

# Fisher's Exact test applied to each subnetwork marker
for (j in 1:length(net)){
    for ( i in 1:length(result)){
        dd <- result[[i]][3:length(result[[i]])]
        n1 <- intersect(net[[j]],dd)
        n2 <- setdiff(dd,n1)
        qq <- setdiff(rownames(NEW),result[[i]])
        n3 <- intersect(net[[j]],unlist(qq))
        n4 <- setdiff(qq,n3)
        foa <- matrix(c(length(n1),length(n2),length(n3),length(n4)),nrow=2)#table for fisher test(what havemisunderstood before)
        rownames(foa) <- c(paste("net",j),paste("non-net",j))
        colnames(foa) <- c("BP","non-BP")
        tt <- fisher.test(foa)
        score[j,i] <- as.numeric(tt$p.value)
    }
}

# the procedure for ajusting the p.value 
phighkegg <- score
phighkegg1 <- as.matrix(phighkegg)
adjusthighkegg1 <- t(apply(phighkegg1, 1, function(x)   p.adjust(x,method="fdr",n=length(x)))) # use the fdr method to adjust the P-value

#select the significantly enriched subnetwork in certain biological processes or KEGG pathways 
try <- apply(adjusthighkegg1,1,function(x) lapply(x,function(y) if (y<0.05) {y <-1} else {y <- 0}))
res <- do.call(rbind,try)
heathigh_kegg <- res
mode(heathigh_kegg) <- "numeric"
colnames(heathigh_kegg) <- unlist(lapply(colnames(heathigh_kegg),function(x) unlist(strsplit(x,split="KEGG_"))[2]))
n <- apply(heathigh_kegg,1,function(x) sum(x))
m <- apply(heathigh_kegg,2,function(x) sum(x))
aa <- as.numeric(which(n!=0))
bb <- as.numeric(which(m!=0))

imagehigh <- t(heathigh_kegg[aa,bb])
rownames(imagehigh) <- colnames(res)[bb]
colnames(imagehigh)<-aa

source("MyImageFunction.R") # invoke the function "MyImageFunction.R"
myImagePlot(imagehigh)  


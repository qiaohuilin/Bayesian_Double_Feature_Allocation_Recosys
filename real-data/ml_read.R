rm(list=ls())
tb=read.table("Movielens/ml-1m/ratings.dat",sep=":")
tb=as.matrix(tb[,c(1,3,5)])
colnames(tb) <- NULL
dim1=max(tb[,1])
dim2=max(tb[,2])
len=dim(tb)[1]


Rcpp::sourceCpp("Movielens/ml-readin.cpp")
R=matrix(0, nrow=dim1, ncol=dim2)
R=readintb(tb,R,dim1,dim2,0,len)

R1=R[,order(colSums(R!=0),decreasing = T)][,1:200]
R1=R1[rowSums((R1!=0))>=3,]

save(R1,file="movielensrating.RData")



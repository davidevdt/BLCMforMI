#MISSINGNESS CREATION
N<-500
dset2<-dset
for(i in 1:N){
	cat("dset=",i,"\n")
	attach(dset2[[i]])
	R<-matrix(0,nrow(dset2[[i]]),ncol(dset2[[i]]))
	
	
	R[which(X1==0 & X4==0),3]<-sample(c(1,0),length(which(X1==0 & X4==0)),rep=T,prob=c(.1,.9))
	R[which(X1==0 & X4==1),3]<-sample(c(1,0),length(which(X1==0 & X4==1)),rep=T,prob=c(.025,.975))
	R[which(X1==0 & X4==2),3]<-sample(c(1,0),length(which(X1==0 & X4==2)),rep=T,prob=c(.125,.875))
	R[which(X1==1 & X4==0),3]<-sample(c(1,0),length(which(X1==1 & X4==0)),rep=T,prob=c(.15,.85))
	R[which(X1==1 & X4==1),3]<-sample(c(1,0),length(which(X1==1 & X4==1)),rep=T,prob=c(.075,.925))
	R[which(X1==1 & X4==2),3]<-sample(c(1,0),length(which(X1==1 & X4==2)),rep=T,prob=c(.05,.95))
	R[which(X1==2 & X4==0),3]<-sample(c(1,0),length(which(X1==2 & X4==0)),rep=T,prob=c(.125,.875))
	R[which(X1==2 & X4==1),3]<-sample(c(1,0),length(which(X1==2 & X4==1)),rep=T,prob=c(.2,.8))
	R[which(X1==2 & X4==2),3]<-sample(c(1,0),length(which(X1==2 & X4==2)),rep=T,prob=c(.15,.85))
	
	
	
	R[which(Y==0 & X5==0),4]<-sample(c(1,0),length(which(Y==0 & X5==0)),rep=T,prob=c(.125,.875))
	R[which(Y==0 & X5==1),4]<-sample(c(1,0),length(which(Y==0 & X5==1)),rep=T,prob=c(.1,.9))
	R[which(Y==0 & X5==2),4]<-sample(c(1,0),length(which(Y==0 & X5==2)),rep=T,prob=c(.15,.85))
	R[which(Y==1 & X5==0),4]<-sample(c(1,0),length(which(Y==1 & X5==0)),rep=T,prob=c(.075,.925))
	R[which(Y==1 & X5==1),4]<-sample(c(1,0),length(which(Y==1 & X5==1)),rep=T,prob=c(.15,.85))
	R[which(Y==1 & X5==2),4]<-sample(c(1,0),length(which(Y==1 & X5==2)),rep=T,prob=c(.05,.95))
	R[which(Y==2 & X5==0),4]<-sample(c(1,0),length(which(Y==2 & X5==0)),rep=T,prob=c(.1,.9))
	R[which(Y==2 & X5==1),4]<-sample(c(1,0),length(which(Y==2 & X5==1)),rep=T,prob=c(.175,.825))
	R[which(Y==2 & X5==2),4]<-sample(c(1,0),length(which(Y==2 & X5==2)),rep=T,prob=c(.125,.925))

	
	
	
	detach(dset2[[i]])
	dset2[[i]][which(R[,3]==1),3]<-NA
	dset2[[i]][which(R[,4]==1),4]<-NA

	
}



#Check
J<-ncol(dset2[[1]])
totmiss<-0
unitmiss<-0
varmiss<-rep(0,J)
for(h in 1:N){

	Yy<-dset2[[h]]
	J<-ncol(Yy)
	n<-nrow(Yy)
	Ra<-ifelse(is.na(Yy),1,0)
	totmiss<-totmiss+(sum(Ra)/(n*J))
	unitmiss<-unitmiss+((n-nrow(na.omit(Yy)))/(n))
	varmiss<-varmiss+((apply(Ra,2,sum))/(n))
	
}
#Total missingness
totmiss/N
# % (average) Missing Units 
unitmiss/N
# % (average) Missing Variables
varmiss/N
rm(Yy)


dset<-dset2
rm(dset2)
rm(list=setdiff(ls(),c("dset","N","n","b","d","par1","par2","ests1","ests2")))



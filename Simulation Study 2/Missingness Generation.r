#MISSINGNESS CREATION
N<-200
dset2<-dset
for(i in 1:N){
	cat("dset=",i,"\n")
	attach(dset2[[i]])
	R<-matrix(0,nrow(dset2[[i]]),ncol(dset2[[i]]))
	
	
	R[which(X3==0 & X4==0),2]<-sample(c(1,0),length(which(X3==0 & X4==0)),rep=T,prob=c(.15,.85))
	R[which(X3==0 & X4==1),2]<-sample(c(1,0),length(which(X3==0 & X4==1)),rep=T,prob=c(.05,.95))
	R[which(X3==1 & X4==0),2]<-sample(c(1,0),length(which(X3==1 & X4==0)),rep=T,prob=c(.25,.75))
	R[which(X3==1 & X4==1),2]<-sample(c(1,0),length(which(X3==1 & X4==1)),rep=T,prob=c(.3,.7))
	
	
	
	R[which(X5==0 & Y==0),7]<-sample(c(1,0),length(which(X5==0 & Y==0)),rep=T,prob=c(.3,.7))
	R[which(X5==0 & Y==1),7]<-sample(c(1,0),length(which(X5==0 & Y==1)),rep=T,prob=c(.2,.8))
	R[which(X5==1 & Y==0),7]<-sample(c(1,0),length(which(X5==1 & Y==0)),rep=T,prob=c(.1,.9))
	R[which(X5==1 & Y==1),7]<-sample(c(1,0),length(which(X5==1 & Y==1)),rep=T,prob=c(.35,.65))
	
	
	R[which(X9==0 & X10==0),17]<-sample(c(1,0),length(which(X9==0 & X10==0)),rep=T,prob=c(.3,.7))
	R[which(X9==0 & X10==1),17]<-sample(c(1,0),length(which(X9==0 & X10==1)),rep=T,prob=c(.25,.75))
	R[which(X9==1 & X10==0),17]<-sample(c(1,0),length(which(X9==1 & X10==0)),rep=T,prob=c(.1,.9))
	R[which(X9==1 & X10==1),17]<-sample(c(1,0),length(which(X9==1 & X10==1)),rep=T,prob=c(.4,.6))

	
	
	R[which(X14==0 & X15==0),21]<-sample(c(1,0),length(which(X14==0 & X15==0)),rep=T,prob=c(.35,.65))
	R[which(X14==0 & X15==1),21]<-sample(c(1,0),length(which(X14==0 & X15==1)),rep=T,prob=c(.1,.9))
	R[which(X14==1 & X15==0),21]<-sample(c(1,0),length(which(X14==1 & X15==0)),rep=T,prob=c(.1,.9))
	R[which(X14==1 & X15==1),21]<-sample(c(1,0),length(which(X14==1 & X15==1)),rep=T,prob=c(.45,.55))


	
		
	detach(dset2[[i]])
	dset2[[i]][which(R[,2]==1),2]<-NA
	dset2[[i]][which(R[,7]==1),7]<-NA
	dset2[[i]][which(R[,17]==1),17]<-NA
	dset2[[i]][which(R[,21]==1),21]<-NA
	
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
totmiss/N
unitmiss/N
varmiss/N
rm(Yy)


dset<-dset2
rm(dset2)
rm(list=setdiff(ls(),c("dset","N","n","b","d","par1","par2","ests1","ests2")))



library(mice)
library(foreach)
library(doParallel)


#Function for combining results obtained in different PC cores
comb <- function(x,...) {
  lapply(seq_along(x),
    function(j) c(x[[j]], lapply(list(...), function(y) y[[j]])))
}



no_cores <- detectCores() 
cl<-makeCluster(no_cores)
registerDoParallel(cl)




#Simulation Preparation
B<-N
npar<-length(b)
vobs<-n-(npar*2)
mm<-20
bp<-b
dp<-d
MAX_IT=20
CRa1<-CRa2<-rep(0,npar)
BIASa1<-BIASa2<-rep(0,npar)



#Parallel Simulations
clusterExport(cl,c("dset","npar","vobs","mm","MAX_IT"))
ptm <- proc.time()
results<-foreach(h=1:B,.combine=comb,.multicombine=TRUE,.init=list(list(),list(),list(),list(),list(),list(),list(),list()),.packages=c('mice',"nnet")) %dopar% {


	para1<-para2<-matrix(0,mm,npar)
	tsea1<-tsea2<-matrix(0,mm,npar)
	ra1<-ra2<-rep(0,npar)
	ua1<-ua2<-rep(0,npar)
	Ba1<-Ba2<-rep(0,npar)
	Ta1<-Ta2<-matrix(0,npar)
	lamdaa1<-lamdaa2<-rep(0,npar)
	nia1<-nia2<-rep(0,npar)
	DOFa1<-DOFa2<-rep(0,npar)
	LOWa1<-LOWa2<-rep(0,npar)		
	UPPa1<-UPPa2<-rep(0,npar)	

	J<-ncol(dset[[h]])
	datt<-dset[[h]]
	for(i in 1:J){
		if(any(is.na(dset[[1]][,i]))){
			datt[,i]<-as.factor(datt[,i])
		}
	}
	
	bb<-mice(datt,m=20,method="polyreg",maxit=MAX_IT,print=TRUE)
	
	dataset<-vector("list",20)
	
	for (i in 1:20) {
		dataset[[i]]<-dset[[h]]
		for(j in 1:J){	
			if(!is.null(bb$imp[[j]])){		
				dataset[[i]][as.numeric(row.names(bb$imp[[j]])),j]<-as.numeric(bb$imp[[j]][,i])		
			}
		}
		daat<-dataset[[i]]	
		daat<-as.data.frame(daat)
		colnames(daat)<-colnames(dset[[h]])
		for(hh in 1:J){
			daat[,hh]<-as.numeric(daat[,hh])
		}
		moda<-multinom(as.factor(Y)~X1+X2+X3+X4+X5+X2:X5+X3:X4,dat=daat)		
		
		para1[i,]<-coefficients(moda)[1,]
		para2[i,]<-coefficients(moda)[2,]		

		tsea1[i,]<-(summary(moda)[[30]][1,])^2
		tsea2[i,]<-(summary(moda)[[30]][2,])^2	
		
	
	}

	
	esta1<-apply(para1,2,mean)
	esta2<-apply(para2,2,mean)	

	ua1<-apply(tsea1,2,mean)
	ua2<-apply(tsea2,2,mean)	

	Ba1<-(apply((t(t(para1)-esta1))^2,2,sum))/(mm-1)			
	Ba2<-(apply((t(t(para2)-esta2))^2,2,sum))/(mm-1)			

	Ta1<-sqrt(ua1+((1+(1/mm))*Ba1))
	Ta2<-sqrt(ua2+((1+(1/mm))*Ba2))	

	lambdaa1<-((1+(1/mm))*Ba1)/(Ta1^2)	
	lambdaa2<-((1+(1/mm))*Ba2)/(Ta2^2)		

	ra1<-(mm-1)/(lambdaa1)^2
	ra2<-(mm-1)/(lambdaa2)^2	

	nia1<-((vobs+1)/(vobs+3))*vobs*(1-lambdaa1)
	nia2<-((vobs+1)/(vobs+3))*vobs*(1-lambdaa2)	

	DOFa1<-((1/ra1)+(1/nia1))^(-1)
	DOFa2<-((1/ra2)+(1/nia2))^(-1)	

	LOWa1<-esta1-(qt(0.975,DOFa1)*Ta1)

	LOWa2<-esta2-(qt(0.975,DOFa2)*Ta2)

	
	
	UPPa1<-esta1+(qt(0.975,DOFa1)*Ta1)
	UPPa2<-esta2+(qt(0.975,DOFa2)*Ta2)

	
	list(esta1,esta2,Ta1,Ta2,LOWa1,LOWa2,UPPa1,UPPa2)
		
	


}
ptmf<-proc.time()-ptm

stopCluster(cl)

#Unlist Simulation Results
esta1<-matrix(unlist(results[[1]]),B,npar,byrow=TRUE)
esta2<-matrix(unlist(results[[2]]),B,npar,byrow=TRUE)
Ta1<-matrix(unlist(results[[3]]),B,npar,byrow=TRUE)
Ta2<-matrix(unlist(results[[4]]),B,npar,byrow=TRUE)
LOWa1<-matrix(unlist(results[[5]]),B,npar,byrow=TRUE)
LOWa2<-matrix(unlist(results[[6]]),B,npar,byrow=TRUE)
UPPa1<-matrix(unlist(results[[7]]),B,npar,byrow=TRUE)
UPPa2<-matrix(unlist(results[[8]]),B,npar,byrow=TRUE)


##################

BIASa1<-(apply(esta1,2,mean))-bp
BIASa2<-(apply(esta2,2,mean))-dp

ASEa1<-apply(Ta1,2,mean)
ASEa2<-apply(Ta2,2,mean)

for(i in 1:B){
	for(j in 1:npar){	
		if(bp[j]>=LOWa1[i,j] & bp[j]<=UPPa1[i,j]){
			CRa1[j] = CRa1[j]+1
		}
	}
}
for(i in 1:B){
	for(j in 1:npar){	
		if(dp[j]>=LOWa2[i,j] & dp[j]<=UPPa2[i,j]){
			CRa2[j] = CRa2[j]+1
		}		
	}
}

CRa1<-CRa1/B
CRa2<-CRa2/B



#Final Results
round(rbind(BIASa1),3)
round(rbind(BIASa2),3)
round(rbind(BIASa1/b),3)
round(rbind(BIASa2/d),3)
round(rbind(ASEa1),3)
round(rbind(ASEa2),3)
rbind(CRa1)
rbind(CRa2)
ptmf








library(mice)
library(foreach)
library(doParallel)



no_cores <- detectCores() 
cl<-makeCluster(no_cores)
registerDoParallel(cl)


#Function to combine results from different PC cores
comb <- function(x,...) {
  lapply(seq_along(x),
    function(j) c(x[[j]], lapply(list(...), function(y) y[[j]])))
}


#Simulation options
B<-N
npar<-length(b)
vobs<-n-(npar)
mm<-20
bp<-b
MAX_IT=20
CRa1<-rep(0,npar)



#Simulations, in parallel
clusterExport(cl,c("dset","npar","vobs","mm","MAX_IT"))
ptm <- proc.time()
results<-foreach(h=1:B,.combine=comb,.multicombine=TRUE,.init=list(list(),list(),list(),list()),.packages=c('mice')) %dopar% {


	para1<-matrix(0,mm,npar)
	tsea1<-matrix(0,mm,npar)
	ra1<-rep(0,npar)	
	ua1<-rep(0,npar)
	Ba1<-rep(0,npar)
	Ta1<-matrix(0,npar)
	lamdaa1<-rep(0,npar)
	nia1<-rep(0,npar)
	DOFa1<-rep(0,npar)
	LOWa1<-rep(0,npar)		
	UPPa1<-rep(0,npar)	

	J<-ncol(dset[[h]])
	datt<-dset[[h]]
	for(i in 1:J){
		if(any(is.na(dset[[h]][,i]))){
			datt[,i]<-as.factor(datt[,i])
		}
	}
	
	bb<-mice(datt,m=20,defaultMethod=rep("polyreg",4),maxit=MAX_IT,print=TRUE)
	
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
		moda<-glm(as.factor(Y)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X1:X5+X1:X17+X1:X5:X17,family='binomial',data=daat)		
		
		para1[i,]<-coefficients(moda)
		
		tsea1[i,]<-(summary(moda)[[12]][,2])^2	
		
	
	}

	
	esta1<-apply(para1,2,mean)

	ua1<-apply(tsea1,2,mean)

	Ba1<-(apply((t(t(para1)-esta1))^2,2,sum))/(mm-1)			

	Ta1<-sqrt(ua1+((1+(1/mm))*Ba1))

	lambdaa1<-((1+(1/mm))*Ba1)/(Ta1^2)		

	ra1<-(mm-1)/(lambdaa1)^2

	nia1<-((vobs+1)/(vobs+3))*vobs*(1-lambdaa1)	

	DOFa1<-((1/ra1)+(1/nia1))^(-1)

	LOWa1<-esta1-(qt(0.975,DOFa1)*Ta1)

	
	
	UPPa1<-esta1+(qt(0.975,DOFa1)*Ta1)

	
	list(esta1,Ta1,LOWa1,UPPa1)
		
	


}
ptmf<-proc.time()-ptm


stopCluster(cl)



#Unlist Results
esta1<-matrix(unlist(results[[1]]),B,npar,byrow=TRUE)
Ta1<-matrix(unlist(results[[2]]),B,npar,byrow=TRUE)
LOWa1<-matrix(unlist(results[[3]]),B,npar,byrow=TRUE)
UPPa1<-matrix(unlist(results[[4]]),B,npar,byrow=TRUE)

##################

BIASa1<-(apply(esta1,2,mean))-bp


ASEa1<-apply(Ta1,2,mean)

for(i in 1:B){
	for(j in 1:npar){	
		if(bp[j]>=LOWa1[i,j] & bp[j]<=UPPa1[i,j]){
			CRa1[j] = CRa1[j]+1
		}
	}
}


CRa1<-CRa1/B



#Final Results
round(rbind(BIASa1),3)
round(rbind(BIASa1/b),3)
round(rbind(ASEa1),3)
round(apply(esta1,2,sd),3)
round(rbind(CRa1),3)
ptmf







#Note: before running the following code, perform model selection from the script "R Parallel - Bayesian Model Selection.r"

library(MMLCimpute)
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



#Simulation Options
B<-200	
classes<-classes	
npar<-length(b)
vobs<-n-(npar)
mm<-20
bp<-b
CRa1<-rep(0,npar)
BIASa1<-rep(0,npar)

IT1<-1000
IT2<-mm*200
IT3<-1
ITPRINT<-250
PRI1<-80				#Simulation Conditions (re-run simulations by modifying this value): {1,80}
PRIRESP<-0.01			


#Simulations, in parallel
clusterExport(cl,c("dset","classes","npar","vobs","mm","IT1","IT2","IT3","ITPRINT","PRI1","PRIRESP"))
ptm <- proc.time()
results<-foreach(i=1:B,.combine=comb,.multicombine=TRUE,.init=list(list(),list(),list(),list()),.packages=c('MMLCimpute')) %dopar% {

	cat("Dataset number ",i,"\n")	
	
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
	
	a<-multilevelLCMI(dat=dset[[i]],GID=NULL,UID=NULL,var2=NULL,L=1,K=classes[i],it1=IT1,it2=IT2,it3=IT3,it.print=ITPRINT,v=10,I=mm,pri2=1.0,pri1=PRI1,priresp=PRIRESP,priresp2=1.0,random=TRUE,num=TRUE,estimates=FALSE,count=FALSE,plot.loglik=FALSE,prec=4,restrict=FALSE)
	
	for(j in 1:mm){
		moda<-glm(as.factor(Y)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X1:X5+X1:X17+X1:X5:X17,family='binomial',data=a[[3]][[j]])
		
		para1[j,]<-coefficients(moda)
		tsea1[j,]<-(summary(moda)[[12]][,2])^2	
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

#Unlist the results
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

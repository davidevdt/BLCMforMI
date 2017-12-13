#Note: use the file containing the number of classes selected with the BLC model!

#LG set-up
library(foreach)
library(doParallel)
library(nnet)

no_cores <- detectCores() 


#Function to combine results from different PC cores
comb <- function(x,...) {
  lapply(seq_along(x),
    function(j) c(x[[j]], lapply(list(...), function(y) y[[j]])))
}



#LatentGOLD 5.1 file folder and Syntax function
LG<-'...//lg51.exe'				#LatentGOLD folder
			
	
makeNewSyntax = function(in_file,out_file,M,K){
paste("//LG5.1//
version = 5.1
infile '",in_file,"'

model
options
	 algorithm 
 tolerance=1e-008 emtolerance=0.01 emiterations=5000 nriterations=0;
 startvalues
 seed=0 sets=100 tolerance=1e-005 iterations=250;
 bayes
 categorical=1 variances=1 latent=1 poisson=1;
 missing  includeall;
 output profile;
 outfile '",out_file,"' imputation= ",M," ;
 variables
  dependent Y,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20 nominal;
  latent
   Z nominal ",K," ;
equations
	Z  <- 1;
	Y  <- 1+Z;
	X1 <- 1+Z;
	X2 <- 1+Z;
	X3 <- 1+Z;
	X4 <- 1+Z;
	X5 <- 1+Z;
	X6 <- 1+Z;
	X7 <- 1+Z;
	X8 <- 1+Z;
	X9 <- 1+Z;
	X10 <- 1+Z;
	X11 <- 1+Z;
	X12 <- 1+Z;
	X13 <- 1+Z;
	X14 <- 1+Z;
	X15 <- 1+Z;
	X16 <- 1+Z;
	X17 <- 1+Z;
	X18 <- 1+Z;
	X19 <- 1+Z;
	X20 <- 1+Z;
end model
",sep="")

}




#HERE: SELECT FROM R-CONSOLE THE FOLDER WHERE YOU WANT TO STORE AND READ THE LG-FILES 



#Simulation Options

B<-N
classes<-classes+0			#{0,20}

npar<-length(b)
vobs<-n-(npar)
mm<-20
bp<-b
CRa1<-rep(0,npar)
BIASa1<-rep(0,npar)



cl<-makeCluster(no_cores)
registerDoParallel(cl)
#Simulations, in parallel
clusterExport(cl,c("dset","classes","npar","vobs","mm","makeNewSyntax"))

#Running the simulations
results<-foreach(i=1:B,.combine=comb,.multicombine=TRUE,.init=list(list(),list(),list(),list())) %dopar% {


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
	
	
	
	thrID<-Sys.getpid()
	in_file<-paste("parallel",thrID,".txt",sep="")
	imp_dat_file<-paste("imputed_data",thrID,".dat",sep="")
	outfile3<-paste("lc_imp",thrID,".lgs",sep="")
	
	write.table(dset[[i]],in_file,na=".",sep=" ",row.names=FALSE,quote=FALSE)
	write.table(makeNewSyntax(in_file,imp_dat_file,mm,classes[i]),outfile3,row.names=FALSE,quote=FALSE,col.names=FALSE)
	
	T1<-proc.time()
	
	shell(paste(LG,outfile3,"/b"))
	

	imp_dat<-read.table(imp_dat_file,sep="",header=TRUE)
	
	for(j in 1:mm){
		
		tmp = imp_dat[which(imp_dat[,22]==j),-22]
		moda<-glm(as.factor(Y)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X1:X5+X1:X17+X1:X5:X17,family='binomial',data=tmp)	
		
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
	
	T2<-proc.time()-T1
	
	list(esta1,Ta1,LOWa1,UPPa1,T2[[3]])
	
	

	
}


stopCluster(cl)

#Unlist the results
esta1<-matrix(unlist(results[[1]]),B,npar,byrow=TRUE)
Ta1<-matrix(unlist(results[[2]]),B,npar,byrow=TRUE)
LOWa1<-matrix(unlist(results[[3]]),B,npar,byrow=TRUE)
UPPa1<-matrix(unlist(results[[4]]),B,npar,byrow=TRUE)
Times<-unlist(results[[5]])

#####################################################################



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
sum(Times)






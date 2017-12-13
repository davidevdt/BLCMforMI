#FIND K
library(MMLCimpute)
library(foreach)
library(doParallel)
B<-N


no_cores <- detectCores() 
cl<-makeCluster(no_cores)
registerDoParallel(cl)



#Model Selection, in parallel
clusterExport(cl,c("dset"))
#Find Number of classes
ptm <- proc.time()
classes<-foreach(i=1:B,.combine='c',.packages='MMLCimpute',.inorder=TRUE) %dopar% {
	a0<-multilevelLCMI(dat=dset[[i]],GID=NULL,UID=NULL,var2=NULL,L=1,K=50,it1=1000,it2=2000,it3=1,it.print=250,v=10,I=0,pri1=1/50,random=TRUE,num=TRUE,estimates=FALSE,count=TRUE,plot.loglik=FALSE,prec=4,restrict=FALSE,scale=0.01)
	
	max(which(a0[[14]]!=0))
	
}
ptmf_mod_selection<-proc.time()-ptm
stopCluster(cl)

#Detected Classes
classes




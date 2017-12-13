library(nnet)


#Bootstrap replications and Sample Size
N=500
n=5000 	#{300,1000}

#INDEPENDENT VARIABLES - loglinear Model
a<-c(0,1,2)
IND<-expand.grid(a,a,a,a,a)
colnames(IND)<-c("X1","X2","X3","X4","X5")
attach(IND)
logp0<-0.5*(apply(IND,1,sum))-1*(X1*X2+X1*X3+X1*X4+X1*X5+X2*X3+X2*X4+X2*X5+X3*X4+X3*X5+X4*X5)-0.2*(X1*X3*X5)+0.5*(X2*X4*X5)
b0<--log(sum(exp(logp0)))
lp<-b0+logp0
rm(logp0)
detach(IND)



#CREATION OF THE DATASETS
S=5
dset<-vector("list",N)

for(i in 1:N){
	dset[[i]]<-matrix(0,n,S)
	ind<-sample(1:dim(IND)[1],n,rep=T,prob=exp(lp))
	dset[[i]]<-IND[ind,]
	rownames(dset[[i]])<-1:n
	dset[[i]]<-as.data.frame(dset[[i]])
	colnames(dset[[i]])<-c("X1","X2","X3","X4","X5")
}
rm(IND)


#MULTI LOGISTIC MODEL PARAMETERS
b0<--0.1
b1<-1
b2<--1.7
b3<-1.5
b4<--0.6
b5<-.5
b25<--0.25
b34<-0.1


d0<--0.6
d1<-1.8
d2<--1.25
d3<-1
d4<-1
d5<--0.5
d25<--0.5
d34<-0.2


#Creation of the Outcome
for(i in 1:N){
	attach(dset[[i]])
	eta1<-b0+b1*X1+b2*X2+b3*X3+b4*X4+b5*X5+b25*X2*X5+b34*X3*X4
	eta2<-d0+d1*X1+d2*X2+d3*X3+d4*X4+d5*X5+d25*X2*X5+d34*X3*X4
	detach(dset[[i]])
	p1<-round((exp(eta1)/(1+exp(eta1)+exp(eta2))),5)
	p2<-round((exp(eta2)/(1+exp(eta1)+exp(eta2))),5)
	p0<-1-(p1+p2)
	Y<-apply(cbind(p0,p1,p2),1,function(x) sample(c(0,1,2),1,rep=T,prob=x))
	dset[[i]]<-as.data.frame(cbind(Y,dset[[i]]))
	colnames(dset[[i]])<-c("Y","X1","X2","X3","X4","X5")
	

}
rm(Y)




#Complete-data simulation results 
b<-c(b0,b1,b2,b3,b4,b5,b25,b34)
d<-c(d0,d1,d2,d3,d4,d5,d25,d34)
par1<-matrix(0,N,length(b))
par2<-matrix(0,N,length(d))
spar1<-matrix(0,N,length(b))
spar2<-matrix(0,N,length(d))
dof<-n-(length(b)*2)
covR1<-rep(0,length(b))
covR2<-rep(0,length(d))
for(i in 1:N){
	attach(dset[[i]])
	cat("dset=",i,"\n")
	mod<-suppressMessages(multinom(as.factor(Y)~X1+X2+X3+X4+X5+X2:X5+X3:X4,dat=dset[[i]]))
	par1[i,] = coefficients(mod)[1,]
	par2[i,] = coefficients(mod)[2,]
	spar1[i,]<-summary(mod)[[30]][1,]
	spar2[i,]<-summary(mod)[[30]][2,]
	for(k in 1:length(b)){
		if(b[k]>par1[i,k]-(qt(0.975,dof)*spar1[i,k])&b[k]<par1[i,k]+(qt(0.975,dof)*spar1[i,k])){
			covR1[k]<-covR1[k]+1
		}
		
		if(d[k]>par2[i,k]-(qt(0.975,dof)*spar2[i,k])&d[k]<par2[i,k]+(qt(0.975,dof)*spar2[i,k])){
			covR2[k]<-covR2[k]+1
		}
	}
	detach(dset[[i]])	
}


#Results
ests1<-apply(par1,2,mean)
ests2<-apply(par2,2,mean)
ests1-b
ests2-d
apply(spar1,2,mean)
apply(spar2,2,mean)
covR1/N
covR2/N

ind<-which(b==0)
rbb<-round((ests1-b)/b,2)
rbd<-round((ests2-d)/d,2)


rm(list=setdiff(ls(),c("dset","ests1","ests2","spar1","spar2","covR1","covR2","ind","rbb","rbd","abb","abd","N","n","b","d","par1","par2")))

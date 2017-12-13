#Bootstrap replications and sample size
N=200
n=2000 	



#INDEPENDENT VARIABLES - loglinear Model
a<-c(0,1,2)
z<-c(0,1)
IND<-expand.grid(z,z,z,z,z,z,z,z,z,z,z,z,z,z,z)
colnames(IND)<-c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15")
attach(IND)
logp0<--0.15*(apply(IND,1,sum))+0.5*(X1*X2+X1*X3+X1*X4+X1*X5+X2*X3+X2*X4+X2*X5+X3*X4+X3*X5+X4*X5)-
0.1*(X6*X7+X6*X8+X6*X9+X6*X10+X6*X11+X7*X8+X7*X9+X7*X10+X7*X11+X8*X9+X8*X10+X8*X11+X9*X10+X9*X11+X10*X11)+0.15*
(X12*X13+X12*X14+X12*X15+X13*X14+X13*X15+X14*X15)+0.3*(X1*X2*X7)+0.6*(X3*X4*X8)-0.4*(X6*X9*X10)
b0<--log(sum(exp(logp0)))
lp<-b0+logp0
sum(exp(lp))
boxplot(exp(lp))
hist(sample(1:dim(IND)[1],(500*300),rep=T,prob=exp(lp)),nclass=dim(IND)[1])
sort(table(sample(1:dim(IND)[1],(200*300),rep=T,prob=exp(lp))),decreasing=T)
length(sort(table(sample(1:dim(IND)[1],(500*300),rep=T,prob=exp(lp))),decreasing=T))
rm(logp0)
detach(IND)

#INDEPENDENT VARIABLES - indipendently created predictors
creat_ind<-function(J,p,n){
	mat<-matrix(-1,n,J)
	for(j in 1:J){
		mat[,j]<-sample(c(0,1),n,rep=T,prob=p[j,])
	}
	mat
}


ps<-matrix(c(
0.3,0.7,
0.4,0.6,
0.5,0.6,
0.4,0.6,
0.3,0.7
),5,2,byrow=TRUE)



#CREATION OF THE DATASETS
S=15
dset<-vector("list",N)

for(i in 1:N){
	dset[[i]]<-matrix(0,n,S)
	ind<-sample(1:dim(IND)[1],n,rep=T,prob=exp(lp))
	dset[[i]]<-IND[ind,]
	dset[[i]]<-cbind(dset[[i]],creat_ind(5,ps,n))
	rownames(dset[[i]])<-1:n
	dset[[i]]<-as.data.frame(dset[[i]])
	colnames(dset[[i]])<-c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20")
}
rm(IND)

	

#Multinomial Logistic coefficients
b0<--1.9
b1<-0.8			
b2<-1.8				
b3<--0.95			
b4<--0.9		
b5<-.8						
b6<-1.1
b7<--0.5
b8<-0.6			
b9<-1			
b10<-0.55
b11<--0.6
b12<-0.75
b13<--1.2
b14<-0
b15<-0
b16<--0.45		
b17<--.85						
b18<-0.55
b19<-0
b20<-0
b1.5<-1.3				
b1.17<--.85				
b1.5.17<-.45	



#Outcome Creation
for(i in 1:N){
	attach(dset[[i]])
	eta1<-b0+b1*X1+b2*X2+b3*X3+b4*X4+b5*X5+b6*X6+b7*X7+b8*X8+b9*X9+b10*X10+b11*X11+b12*X12+b13*X13+b14*X14+b15*X15+b16*X16+b17*X17+b18*X18+b19*X19+b20*X20+(b1.5*X1*X5)+(b1.17*X1*X17)+(b1.5.17*X1*X5*X17)
	detach(dset[[i]])
	p1<-round((exp(eta1)/(1+exp(eta1))),5)
	p0<-1-p1
	Y<-apply(cbind(p0,p1),1,function(x) sample(c(0,1),1,rep=T,prob=x))
	dset[[i]]<-as.data.frame(cbind(Y,dset[[i]]))
	colnames(dset[[i]])<-c("Y","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20")
	

}
rm(Y)





#Simulation on full dataset
b<-c(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b1.5,b1.17,b1.5.17)
par1<-matrix(0,N,length(b))
spar1<-matrix(0,N,length(b))
dof<-n-(length(b))
covR1<-rep(0,length(b))
for(i in 1:N){
	attach(dset[[i]])
	cat("dset=",i,"\n")
	mod<-glm(as.factor(Y)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X1:X5+X1:X17+X1:X5:X17,family='binomial',data=dset[[i]])
	par1[i,] = coefficients(mod)
	spar1[i,]<-summary(mod)[[12]][,2]
	for(k in 1:length(b)){
		if(b[k]>par1[i,k]-(qt(0.975,dof)*spar1[i,k])&b[k]<par1[i,k]+(qt(0.975,dof)*spar1[i,k])){
			covR1[k]<-covR1[k]+1
		}
	}
	detach(dset[[i]])	
}
ests1<-apply(par1,2,mean)
ests1-b
apply(spar1,2,mean)
covR1/N

#Results
ind<-which(b==0)
rbb<-round((ests1[-ind]-b[-ind])/b[-ind],2)
abb<-round(ests1[ind],2)

rm(list=setdiff(ls(),c("dset","ests1","ests2","spar1","spar2","covR1","covR2","ind","rbb","rbd","abb","abd","N","n","b","d","par1","par2")))	

#Final report
setwd('/Users/mac/Documents/IM')
library(data.table)
library(nloptr)
#part1:MVO vs Risk parity
win=108
Assets=fread('asset_china.csv')
Assets
Assets=Assets[(nrow(Assets)-win+1):nrow(Assets),-1,with=FALSE]
AssetCov=cov(Assets)
AssetRtn=apply(Assets,2,mean)
cor(Assets)
w0=rep(1/14,14)
ret=w0%*%AssetRtn#calculate the original portfolio return
ret

#MVO
Vol_tgt=0.02
ObjFun=function(w)
{
  -sum(AssetRtn*w)
}
ObjFunGradient=function(w)
{
  -AssetRtn
}
Equal=function(w)
{
  sum(w)-1
}
EqualGradient=function(w)
{
  rep(1,length(w))
}
Inequal=function(w)
{
  t(w) %*% AssetCov %*% w-Vol_tgt^2 
}
InequalGradient=function(w)
{
  t(w) %*% (AssetCov+t(AssetCov))
}
w0=rep(1/14,14)
w0=w0/sum(w0)
opt=nloptr(x0=w0,eval_f=ObjFun,eval_grad_f=ObjFunGradient,
           eval_g_eq=Equal,eval_jac_g_eq=EqualGradient,
           eval_g_ineq=Inequal,eval_jac_g_ineq=InequalGradient,
           lb=rep(-0.025,14),ub=rep(1,14),
           opts=list('algorithm'='NLOPT_LD_SLSQP','xtol_abs'=1.0e-8,'maxeval'=100000000))
w_MVO=opt$solution
sum(AssetRtn*w_MVO)
sqrt(t(w_MVO) %*% AssetCov %*% w_MVO)
MCR=AssetCov %*% w_MVO / sqrt(t(w_MVO) %*% AssetCov %*% w_MVO)[1,1]
data.frame(w_MVO,w_MVO*MCR)

#risk parity
Assets1=Assets[,-'SHIBOR3m']
AssetRtn1=apply(Assets1,2,mean)
AssetCov1=cov(Assets1)
ObjFun=function(w)
{
  -sum(log(abs(w)))
}
ObjFunGradient=function(w)
{
  -(1/w)
}
Inequal=function(w)
{
  t(w) %*% AssetCov1 %*% w-Vol_tgt^2 
}
InequalGradient=function(w)
{
  t(w) %*% (AssetCov1+t(AssetCov1))
}
w0=rep(1/13,13)
opt=nloptr(x0=w0,eval_f=ObjFun,eval_grad_f=ObjFunGradient,
           eval_g_ineq=Inequal,eval_jac_g_ineq=InequalGradient,
           lb=rep(-1,13),ub=rep(1,13),
           opts=list('algorithm'='NLOPT_LD_SLSQP','xtol_abs'=1.0e-8,'maxeval'=100000000))
w_RP=opt$solution
sum(w_RP)
sum(AssetRtn*c(w_RP,1-sum(w_RP)))
sqrt(t(w_RP) %*% AssetCov1 %*% w_RP)
MCR=AssetCov1 %*% w_RP / sqrt(t(w_RP) %*% AssetCov1 %*% w_RP)[1,1]
data.frame(w_RP,w_RP*MCR)

#factor based asset allocation
#PCA
Ast=fread('asset_pca.csv')
Ast
Ast[,length(unique(ticket))]
Ast=dcast(Ast,date~ticket,value.var = 'return')
Ast
pz=prcomp(Ast[,-1,with=F])
summary(pz)
pz$rotation

#factor model
Factors=fread('factors_china.csv')
Factors=Factors
Factors=Factors[(nrow(Factors)-win+1):nrow(Factors),-1,with=FALSE]
FacCov=cov(Factors)
FacRtn=apply(Factors,2,mean)
cor(Factors)
AssetExpo=matrix(0,ncol(Assets),ncol(Factors))
rownames(AssetExpo)=colnames(Assets)
colnames(AssetExpo)=colnames(Factors)

Spe=matrix(0,nrow(Factors),ncol(Assets))
rownames(Spe)=rownames(Factors)
colnames(Spe)=colnames(Assets)

Mapping=function(y)
{
  dat=data.table(Y=y,Factors)
  full=lm(Y~.,data = dat)
  null=lm(Y~1,data = dat)
  best=step(null,scope = list(upper=full,lower=null),direction = 'both')
  list(coeff=best$coefficients[-1],res=summary(best)$residuals)
}
for (i in 1:ncol(Assets)) {
  a=Mapping(data.frame(Assets[,i,with=FALSE])[,1])
  AssetExpo[i,names(a$coeff)]=a$coeff
  Spe[,i]= a$res
}
AssetExpo
FacCov=cov(Factors)
SpeCov=cov(Spe)

TargetExpo=w_MVO %*% AssetExpo
TargetExpo
lambda=0.99
ObjFun=function(w)
{
  M=w %*% AssetExpo-TargetExpo
  ((1-lambda) * M %*% t(M) + lambda * M %*% FacCov %*% t(M) +lambda * w %*% SpeCov %*% w)[1,1]
    }
ObjFunGradient=function(w)
{
  M=w %*% AssetExpo-TargetExpo
  2*(1-lambda) * M %*% t(AssetExpo) + 2*lambda * M %*% FacCov %*% t(AssetExpo) +2*lambda * w %*% SpeCov
}
fullyInvested=function(w)
{
  sum(w)-1
}
fullyInvestedGradient=function(w)
{
  rep(1,length(w))
}
#case1:weight bound(-1,1)
w0=rep(1/14,14)
opt=nloptr(x0=w0,eval_f=ObjFun,eval_grad_f=ObjFunGradient,
           eval_g_eq=fullyInvested,eval_jac_g_eq=fullyInvestedGradient,
           lb=rep(-1,14),ub=rep(1,14),
           opts=list('algorithm'='NLOPT_LD_SLSQP','xtol_abs'=1.0e-4,'maxeval'=100000000))
w=opt$solution
w
w %*% AssetExpo
TargetExpo
sum(AssetRtn*w)
#case2:weight bound(-0.025,1)
w0=rep(1/14,14)
opt=nloptr(x0=w0,eval_f=ObjFun,eval_grad_f=ObjFunGradient,
           eval_g_eq=fullyInvested,eval_jac_g_eq=fullyInvestedGradient,
           lb=rep(-0.025,14),ub=rep(1,14),
           opts=list('algorithm'='NLOPT_LD_SLSQP','xtol_abs'=1.0e-4,'maxeval'=100000000))
w=opt$solution
w
w %*% AssetExpo
TargetExpo
sum(AssetRtn*w)


#case3:factor risk-parity(target exposure calculated by MVO)
Equal=function(w)
  {
    return(rbind(sum(w)-1,w%*%(AssetExpo[,1]-AssetExpo[,2]),w%*%(AssetExpo[,1]-AssetExpo[,3]),w%*%(AssetExpo[,1]-AssetExpo[,4]),w%*%(AssetExpo[,1]-AssetExpo[,5]),w%*%(AssetExpo[,1]-AssetExpo[,6]),w%*%(AssetExpo[,1]-AssetExpo[,7])))
  }
EqualGradient=function(w)
  {
    return(rbind(rep(1,length(w)),AssetExpo[,1]-AssetExpo[,2],AssetExpo[,1]-AssetExpo[,3],AssetExpo[,1]-AssetExpo[,4],AssetExpo[,1]-AssetExpo[,5],AssetExpo[,1]-AssetExpo[,6],AssetExpo[,1]-AssetExpo[,7]))
  }
w0=rep(1/14,14)
opt=nloptr(x0=w0,eval_f=ObjFun,eval_grad_f=ObjFunGradient,
           eval_g_eq=Equal,eval_jac_g_eq=EqualGradient,
           lb=rep(-1,14),ub=rep(1,14),
           opts=list('algorithm'='NLOPT_LD_SLSQP','xtol_abs'=1.0e-4,'maxeval'=100000000))
w=opt$solution
w
w %*% AssetExpo
sum(AssetRtn*w)
#case4:factor risk-parity(target exposure calculated by Risk-parity)
ObjFun1=function(w)
{
  -sum(log(abs(w)))
}
ObjFunGradient1=function(w)
{
  -(1/w)
}
Inequal=function(w)
{
  t(w) %*% AssetCov1 %*% w-Vol_tgt^2 
}
InequalGradient=function(w)
{
  t(w) %*% (AssetCov1+t(AssetCov1))
}
w0=rep(1/13,13)
opt=nloptr(x0=w0,eval_f=ObjFun1,eval_grad_f=ObjFunGradient1,
           eval_g_ineq=Inequal,eval_jac_g_ineq=InequalGradient,
           lb=rep(-0.025,13),ub=rep(1,13),
           opts=list('algorithm'='NLOPT_LD_SLSQP','xtol_abs'=1.0e-8,'maxeval'=100000000))
w_RP=opt$solution
sum(AssetRtn*c(w_RP,1-sum(w_RP)))
sqrt(t(w_RP) %*% AssetCov1 %*% w_RP)

AssetExpo=matrix(0,ncol(Assets1),ncol(Factors))
rownames(AssetExpo)=colnames(Assets1)
colnames(AssetExpo)=colnames(Factors)

Spe=matrix(0,nrow(Factors),ncol(Assets1))
rownames(Spe)=rownames(Factors)
colnames(Spe)=colnames(Assets1)

for (i in 1:ncol(Assets1)) {
  a=Mapping(data.frame(Assets[,i,with=FALSE])[,1])
  AssetExpo[i,names(a$coeff)]=a$coeff
  Spe[,i]= a$res
}
AssetExpo
FacCov=cov(Factors)
SpeCov=cov(Spe)
TargetExpo=w_RP %*% AssetExpo
TargetExpo

Equal=function(w)
{
  return(rbind(sum(w)-1,w%*%(AssetExpo[,1]-AssetExpo[,2]),w%*%(AssetExpo[,1]-AssetExpo[,3]),w%*%(AssetExpo[,1]-AssetExpo[,4]),w%*%(AssetExpo[,1]-AssetExpo[,5]),w%*%(AssetExpo[,1]-AssetExpo[,6]),w%*%(AssetExpo[,1]-AssetExpo[,7])))
}
EqualGradient=function(w)
{
  return(rbind(rep(1,length(w)),AssetExpo[,1]-AssetExpo[,2],AssetExpo[,1]-AssetExpo[,3],AssetExpo[,1]-AssetExpo[,4],AssetExpo[,1]-AssetExpo[,5],AssetExpo[,1]-AssetExpo[,6],AssetExpo[,1]-AssetExpo[,7]))
}
w0=rep(1/13,13)
opt=nloptr(x0=w0,eval_f=ObjFun,eval_grad_f=ObjFunGradient,
           eval_g_eq=Equal,eval_jac_g_eq=EqualGradient,
           lb=rep(-1,13),ub=rep(1,13),
           opts=list('algorithm'='NLOPT_LD_SLSQP','xtol_abs'=1.0e-4,'maxeval'=100000000))
w=opt$solution
w
w %*% AssetExpo
sum(w)
sum(AssetRtn1*w)


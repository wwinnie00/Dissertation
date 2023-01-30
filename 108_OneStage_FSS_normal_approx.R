###################################################################
### Dunnett's One-Stage for t=1 or 2
### Z_i = sqrt(4*n)*arcsin(sqrt(X_i/n))
### If Z_{k-t+1} >= Z_0+sqrt(2)*c, select the t largest Z_i
### Otherwise, select only Z_i >= Z_0+c*sqrt(1/(2n))
###################################################################

source(file="~/Dropbox/105_OneStage_FSS_PCS_exact.R")
source(file="~/Dropbox/107_OneStage_FSS_search value.R")

### P(CS0), any k and t 
pcs0.norm <- function(p0,n,k,c,delta1){
  A<-function(x){x+sqrt(4*n)*(asin(sqrt(p0))-asin(sqrt(p0+delta1)))+sqrt(2)*c}
  
  integrand <- function(x){dnorm(x)*pnorm(A(x))^k}
  integrate(integrand,-Inf,Inf)$value
}

#pcs0.norm(p0=0.1,n=46,k=3,c=0.4,delta1=0.05)

# pcs0.norm.alt <- function(p0,n,k,c,delta1){
#   ## no integral needed:
#   ## P(CS0)=pnorm(c-sqrt(2*n)*(asin(sqrt(p0+delta1))-asin(sqrt(p0))))
#   ## result slightly different from the original function
#   stat<-c-sqrt(2*n)*(asin(sqrt(p0+delta1))-asin(sqrt(p0)))
#   pnorm(stat)^k
# }
# n<-46
# h<-1.7335
# p0=0.1
# delta1=0.05
# c<-h+sqrt(2*n)*(asin(sqrt(p0+delta1))-asin(sqrt(p0)))
# pcs0.norm(p0,n,k=5,c,delta1)
# pcs0.norm.alt(p0,n,k=3,c,delta1)
# 
# pcs0.norm.hval <- function(k,h){
#   A<-function(x){x+sqrt(2)*h}
#   integrand <- function(x){dnorm(x)*pnorm(A(x))^k}
#   integrate(integrand,-Inf,Inf)$value
# }
# pcs0.norm(3,1.7335)

### P(CS1), t=1
pcs1.norm.t1 <- function(p0,n,k,c,delta1,delta2){
  A<-function(x){x+sqrt(4*n)*(asin(sqrt(p0+delta2))-asin(sqrt(p0)))-sqrt(2)*c}
  B<-function(x){x+sqrt(4*n)*(asin(sqrt(p0+delta2))-asin(sqrt(p0+delta1)))}
  
  integrand<-function(x) {dnorm(x)*pnorm(A(x))*pnorm(B(x))^(k-1)}
  integrate(integrand, -Inf,Inf)$value
}

#pcs1.norm.t1(p0=0.1,n=46,k=3,c=0.4,delta1=0.05,delta2=0.2)

### P(CS1), t=2
pcs1.norm.t2 <- function(p0,n,k,c,delta1,delta2){
  A1<-function(x){x+sqrt(4*n)*(asin(sqrt(p0+delta2))-asin(sqrt(p0)))-sqrt(2)*c}  #integrand function when x<=y
  B1<-function(x){x+sqrt(4*n)*(asin(sqrt(p0+delta2))-asin(sqrt(p0+delta1)))}   #integrand function when x<=y
  A2<-function(y){y+sqrt(4*n)*(asin(sqrt(p0+delta2))-asin(sqrt(p0)))-sqrt(2)*c}  #integrand function when x>y
  B2<-function(y){y+sqrt(4*n)*(asin(sqrt(p0+delta2))-asin(sqrt(p0+delta1)))}   #integrand function when x>y
  
  InnerIntegral1<- function(y){
    InnerFunc<-function(x) {dnorm(x)*dnorm(y)*pnorm(A1(x))*pnorm(B1(x))^(k-2)}
    integrate(InnerFunc, -Inf,y)$value
  }
  InnerIntegral2<- function(y){
    InnerFunc<-function(x) {dnorm(x)*dnorm(y)*pnorm(A2(y))*pnorm(B2(y))^(k-2)}
    integrate(InnerFunc, y+0.00000001,Inf)$value
  }
  integrate(Vectorize(InnerIntegral1) , -Inf, Inf)$value +
    integrate(Vectorize(InnerIntegral2) , -Inf, Inf)$value
}

pcs1.norm.t2(p0=0.2,n=109,k=3,c=2.063,delta1=0,delta2=0.2)


# Search values, P(CS0) and P(CS1) requirements

dun.norm.fvalue <- function(k,t,p0star,p1star,p0, delta1=0.05, delta2=0.2,
                            c_L,c_U,c_step=0.001,n_L,n_U,n_step=1){
  done<-FALSE

  for (c in seq(c_L,c_U,by=c_step)){
    for (n in seq(n_L,n_U,by=n_step)){
      
    if(t==1){
      pcs1 <- pcs1.norm.t1(p0,n,k,c,delta1,delta2)
    }
    else if (t==2){
      pcs1 <- pcs1.norm.t2(p0,n,k,c,delta1,delta2)
    }
    else {stop('Undefined value of t!')}
    if (pcs1>=p1star){
      pcs0 <- pcs0.norm(p0,n,k,c,delta1)
      if(pcs0>=p0star){
        N_val <- (k+1)*n
            
        out<-data.frame(k,t,p0,delta1,delta2,c,n,N_val,pcs0,pcs1)
        print(out)

        done<-TRUE
        break
        }
      else 
          break  ## P(CS0) will only be smaller when increasing n while c is fixed
      }
    }
    if (done){break} ## after a set of (n,c) found, if continue to increase c, the appropriate n will only be larger
  }
  out
}

###################################
#### k=2, p=0.2, P0*=0.95, P1*=0.7
appr.out<-dun.norm.fvalue(k=2,t=1,p0star=0.95,p1star=0.7,p0=0.2,delta1=0,delta2=0.2,
                          c_L=1,c_U=3,c_step=0.001,n_L=40,n_U=70) #n=62,c=1.917

## Compared with exact binomial:
out<-NULL
exact.out <- onestage.design.out(k=2,t=1,n=62,c=10,p0=0.2,delta1=0,delta2=0.2) 
optimal <- RCW1.fvalue(k=2,t=1,p0star=0.95,p1star=0.7,p0=0.2,delta1=0,delta2=0.2,c_L=8,c_U=10,n_L=50,n_U=62,n_step=1) 
tmp<-data.frame(appr.out[1,1:7],P0star=0.95,P1star=0.70,c.bin=exact.out[1,5],PCS0.bin=exact.out[1,7],PCS1.bin=exact.out[1,9],
                opt.c=optimal[1,6],opt.n=optimal[1,7])
out<-rbind(out,tmp)
write.csv(out,file="R_CW1_NormApp_vs_Exact.csv")

#### k=2, p=0.2, P0*=0.95, P1*=0.8
appr.out<-dun.norm.fvalue(k=2,t=1,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,
               c_L=1,c_U=3,c_step=0.001,n_L=70,n_U=90) #n=78,c=1.917

## Compared with exact binomial:
out<-read.csv(file="R_CW1_NormApp_vs_Exact.csv")[,-1]
exact.out <- onestage.design.out(k=2,t=1,n=78,c=11,p0=0.2,delta1=0,delta2=0.2) 
optimal <- RCW1.fvalue(k=2,t=1,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,c_L=8,c_U=11,n_L=55,n_U=70,n_step=1) 
tmp<-data.frame(appr.out[1,1:7],P0star=0.95,P1star=0.80,c.bin=exact.out[1,5],PCS0.bin=exact.out[1,7],PCS1.bin=exact.out[1,9],
                opt.c=optimal[1,6],opt.n=optimal[1,7])
out<-rbind(out,tmp)
write.csv(out,file="R_CW1_NormApp_vs_Exact.csv")

# #### k=2, p=0.1, P0*=0.95, P1*=0.9
# appr.out<-dun.norm.fvalue(k=2,t=1,p0star=0.95,p1star=0.9,p0=0.1,delta1=0,delta2=0.2,
#                           c_L=1,c_U=3,c_step=0.001,n_L=60,n_U=80) #n=77,c=1.917
# 
# ## Compared with exact binomial:
# out<-read.csv(file="R_CW1_NormApp_vs_Exact.csv")[,-1]
# exact.out <- onestage.design.out(k=2,t=1,n=77,c=8,p0=0.1,delta1=0,delta2=0.2) 
# tmp<-data.frame(appr.out[1,1:7],P0star=0.95,P1star=0.90,c.bin=exact.out[1,5],PCS0.bin=exact.out[1,7],PCS1.bin=exact.out[1,9])
# out<-rbind(out,tmp)
# write.csv(out,file="R_CW1_NormApp_vs_Exact.csv")

#### k=3, t=1, p=0.2, P0*=0.95, P1*=0.7
appr.out<-dun.norm.fvalue(k=3,t=1,p0star=0.95,p1star=0.7,p0=0.2,delta1=0,delta2=0.2,
                          c_L=1,c_U=3,c_step=0.001,n_L=60,n_U=120) #n=69,c=2.063

## Compared with exact binomial:
out<-read.csv(file="R_CW1_NormApp_vs_Exact.csv")[,-1]
exact.out <- onestage.design.out(k=3,t=1,n=69,c=11,p0=0.2,delta1=0,delta2=0.2) 
optimal <- RCW1.fvalue(k=3,t=1,p0star=0.95,p1star=0.7,p0=0.2,delta1=0,delta2=0.2,c_L=8,c_U=10,n_L=60,n_U=80,n_step=1) 
tmp<-data.frame(appr.out[1,1:7],P0star=0.95,P1star=0.70,c.bin=exact.out[1,5],PCS0.bin=exact.out[1,7],PCS1.bin=exact.out[1,9],
                opt.c=optimal[1,6],opt.n=optimal[1,7])
out<-rbind(out,tmp)
write.csv(out,file="R_CW1_NormApp_vs_Exact.csv")

#### k=3, t=1, p=0.2, P0*=0.95, P1*=0.8
appr.out<-dun.norm.fvalue(k=3,t=1,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,
                          c_L=1,c_U=3,c_step=0.001,n_L=60,n_U=120) #n=87,c=2.063

## Compared with exact binomial:
out<-read.csv(file="R_CW1_NormApp_vs_Exact.csv")[,-1]
exact.out <- onestage.design.out(k=3,t=1,n=87,c=12,p0=0.2,delta1=0,delta2=0.2) 
optimal <- RCW1.fvalue(k=3,t=1,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,c_L=8,c_U=12,n_L=70,n_U=90,n_step=1) 
tmp<-data.frame(appr.out[1,1:7],P0star=0.95,P1star=0.80,c.bin=exact.out[1,5],PCS0.bin=exact.out[1,7],PCS1.bin=exact.out[1,9],
                opt.c=optimal[1,6],opt.n=optimal[1,7])
out<-rbind(out,tmp)
write.csv(out,file="R_CW1_NormApp_vs_Exact.csv")

#### k=4, t=1, p=0.2, P0*=0.95, P1*=0.7
appr.out<-dun.norm.fvalue(k=4,t=1,p0star=0.95,p1star=0.7,p0=0.2,delta1=0,delta2=0.2,
                          c_L=1,c_U=3,c_step=0.001,n_L=60,n_U=120) #n=75,c=2.161

## Compared with exact binomial:
out<-read.csv(file="R_CW1_NormApp_vs_Exact.csv")[,-1]
exact.out <- onestage.design.out(k=3,t=1,n=75,c=11,p0=0.2,delta1=0,delta2=0.2) 
optimal <- RCW1.fvalue(k=4,t=1,p0star=0.95,p1star=0.7,p0=0.2,delta1=0,delta2=0.2,c_L=9,c_U=12,n_L=60,n_U=90,n_step=1) 
tmp<-data.frame(appr.out[1,1:7],P0star=0.95,P1star=0.70,c.bin=exact.out[1,5],PCS0.bin=exact.out[1,7],PCS1.bin=exact.out[1,9],
                opt.c=optimal[1,6],opt.n=optimal[1,7])
out<-rbind(out,tmp)
write.csv(out,file="R_CW1_NormApp_vs_Exact.csv")

#### k=4, t=1, p=0.2, P0*=0.95, P1*=0.8
appr.out<-dun.norm.fvalue(k=4,t=1,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,
                          c_L=2,c_U=3,c_step=0.001,n_L=60,n_U=120) #n=93,c=2.161

## Compared with exact binomial:
out<-read.csv(file="R_CW1_NormApp_vs_Exact.csv")[,-1]
exact.out <- onestage.design.out(k=3,t=1,n=93,c=12,p0=0.2,delta1=0,delta2=0.2) 
optimal <- RCW1.fvalue(k=4,t=1,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,c_L=9,c_U=12,n_L=60,n_U=90,n_step=1) 
tmp<-data.frame(appr.out[1,1:7],P0star=0.95,P1star=0.80,c.bin=exact.out[1,5],PCS0.bin=exact.out[1,7],PCS1.bin=exact.out[1,9],
                opt.c=optimal[1,6],opt.n=optimal[1,7])
out<-rbind(out,tmp)
write.csv(out,file="R_CW1_NormApp_vs_Exact.csv")

#### k=3, t=2, p=0.2, P0*=0.95, P1*=0.7
appr.out<-dun.norm.fvalue(k=3,t=2,p0star=0.95,p1star=0.7,p0=0.2,delta1=0,delta2=0.2,
                          c_L=2,c_U=3,c_step=0.001,n_L=80,n_U=110) #n=89,c=2.063

## Compared with exact binomial:
out<-read.csv(file="R_CW1_NormApp_vs_Exact.csv")[,-1]
exact.out <- onestage.design.out(k=3,t=2,n=89,c=12,p0=0.2,delta1=0,delta2=0.2) 
optimal <- RCW1.fvalue(k=3,t=2,p0star=0.95,p1star=0.7,p0=0.2,delta1=0,delta2=0.2,c_L=11,c_U=13,n_L=70,n_U=105,n_step=1) 
tmp<-data.frame(appr.out[1,1:7],P0star=0.95,P1star=0.70,c.bin=exact.out[1,5],PCS0.bin=exact.out[1,7],PCS1.bin=exact.out[1,9],
                opt.c=optimal[1,6],opt.n=optimal[1,7])
out<-rbind(out,tmp)
write.csv(out,file="R_CW1_NormApp_vs_Exact.csv")

#### k=3, t=2, p=0.2, P0*=0.95, P1*=0.8
appr.out<-dun.norm.fvalue(k=3,t=2,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,
                          c_L=2,c_U=3,c_step=0.001,n_L=90,n_U=110) #n=107,c=2.063

## Compared with exact binomial:
out<-read.csv(file="R_CW1_NormApp_vs_Exact.csv")[,-1]
exact.out <- onestage.design.out(k=3,t=2,n=107,c=13,p0=0.2,delta1=0,delta2=0.2) 
optimal <- RCW1.fvalue(k=3,t=2,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,c_L=12,c_U=14,n_L=90,n_U=100,n_step=1) 
tmp<-data.frame(appr.out[1,1:7],P0star=0.95,P1star=0.80,c.bin=exact.out[1,5],PCS0.bin=exact.out[1,7],PCS1.bin=exact.out[1,9],
                opt.c=optimal[1,6],opt.n=optimal[1,7])
out<-rbind(out,tmp)
write.csv(out,file="R_CW1_NormApp_vs_Exact.csv")

#### k=4, t=2, p=0.2, P0*=0.95, P1*=0.7
appr.out<-dun.norm.fvalue(k=4,t=2,p0star=0.95,p1star=0.7,p0=0.2,delta1=0,delta2=0.2,
                          c_L=2,c_U=3,c_step=0.001,n_L=80,n_U=100) #n=95,c=2.161

## Compared with exact binomial:
out<-read.csv(file="R_CW1_NormApp_vs_Exact.csv")[,-1]
exact.out <- onestage.design.out(k=4,t=2,n=95,c=13,p0=0.2,delta1=0,delta2=0.2) 
optimal <- RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.7,p0=0.2,delta1=0,delta2=0.2,c_L=11,c_U=13,n_L=80,n_U=100,n_step=1) 
tmp<-data.frame(appr.out[1,1:7],P0star=0.95,P1star=0.70,c.bin=exact.out[1,5],PCS0.bin=exact.out[1,7],PCS1.bin=exact.out[1,9],
                opt.c=optimal[1,6],opt.n=optimal[1,7])
out<-rbind(out,tmp)
write.csv(out,file="R_CW1_NormApp_vs_Exact.csv")

#### k=4, t=2, p=0.2, P0*=0.95, P1*=0.8
appr.out<-dun.norm.fvalue(k=4,t=2,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,
                          c_L=2.1,c_U=3,c_step=0.001,n_L=110,n_U=125) #n=114,c=2.161

## Compared with exact binomial:
out<-read.csv(file="R_CW1_NormApp_vs_Exact.csv")[,-1]
exact.out <- onestage.design.out(k=4,t=2,n=114,c=14,p0=0.2,delta1=0,delta2=0.2) 
optimal <- RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,c_L=13,c_U=14,n_L=100,n_U=110,n_step=1) 
tmp<-data.frame(appr.out[1,1:7],P0star=0.95,P1star=0.80,c.bin=exact.out[1,5],PCS0.bin=exact.out[1,7],PCS1.bin=exact.out[1,9],
                opt.c=optimal[1,6],opt.n=optimal[1,7])
out<-rbind(out,tmp)
write.csv(out,file="R_CW1_NormApp_vs_Exact.csv")

################################################################################################
# Simulation, t=2, One-stage Normal Appr. 

indicator<-function(condition) 
{ifelse(condition,1,0)}

#### K=3, P(CS0) ####
p.cs0.sim <- function(nsim,n,c,p0,delta1){
  k=3; t=2
  set.seed(752)
  count <- 0
  for (l in 1:nsim){
    
    z0 <- rnorm(1,mean=asin(sqrt(p0)),sd=1/(2*sqrt(n)))
    z1 <- rnorm(1,mean=asin(sqrt(p0+delta1)),sd=1/(2*sqrt(n)))
    z2 <- rnorm(1,mean=asin(sqrt(p0+delta1)),sd=1/(2*sqrt(n)))
    z3 <- rnorm(1,mean=asin(1),sd=1/(2*sqrt(n)))
    z.vec <- c(z1,z2,z3)
    test <- sort(z.vec, TRUE)[t]  #find the 2nd largest order statistic
    if (test<=z0+c*sqrt(1/(2*n))){
      count=count+1
    }
  }
  p.cs0 <- count/nsim
  p.cs0
}

p.cs0.sim(nsim=1000000,n=58,c=1.4,p0=0.2,delta1=0.05)

#### K=3, P(CS1) ####
p.cs1.sim <- function(nsim,n,c,p0,delta1,delta2){
  k=3; t=2
  set.seed(752)
  count <- 0
  for (l in 1:nsim){
    
    z0 <- rnorm(1,mean=asin(sqrt(p0)),sd=1/(2*sqrt(n)))
    z1 <- rnorm(1,mean=asin(sqrt(p0+delta1)),sd=1/(2*sqrt(n)))
    z2 <- rnorm(1,mean=asin(sqrt(p0+delta2)),sd=1/(2*sqrt(n)))
    z3 <- rnorm(1,mean=asin(sqrt(p0+delta2)),sd=1/(2*sqrt(n)))
    
    if (min(z2,z3)>z0+c*sqrt(1/(2*n)) && z1 < min(z2,z3)){
      count=count+1
    }
  }
  p.cs0 <- count/nsim
  p.cs0
}

p.cs1.sim(nsim=100000,p0=0.2,n=58,c=1.400,delta1=0.05,delta2=0.2)




#################################################################################################
######## Subset selection one-stage design allowing less than t, optimal design using exact calc
######## Goal: select a subset containing the t best given that at least t are acceptable; 
########        reject marginal populations if g<t are acceptable.
######## Procedure: (d>=0)
######## if X_(k-t+1)-X0>=c, select all populations for which  Xi>=X_(k-t+1)-d
######## else, reject all populations for which X_i<X0+c
#################################################################################################
## Original version: 02NOV2022
library(dplyr)
source(file="~/Dropbox/AAAA Paper/Codes/Ch4_Subset Selection Approach/419_R_CW4_mod_PCS_ES_exact.R")
source(file="~/Dropbox/AAAA Paper/Codes/Ch1_One stage FSS/107_OneStage_FSS_search value.R")
source(file="~/Dropbox/AAAA Paper/Codes/Ch4_Subset Selection Approach/418_R_CW4_mod_simulation.R")

setwd("~/Dropbox/AAAA Paper/Codes/Ch4_Subset Selection Approach/Outputs/420_R_CW4_mod_optimal_design_exact")

RCW4.fvalue <- function(k,t,p0star,p1star,p0, delta1=0, delta2=0.2,
                        c_L,c_U,c_step=1,d_L,d_U,d_step=1,n_L,n_U,n_step=1){
  out<-NULL
  done <- FALSE
  for (c in seq(c_L,c_U,by=c_step)){
    for (n in seq(n_L,n_U,by=n_step)){
      pcs0 <- RCW1.pcs.seq(k,t,n,c,p0,delta1,delta2)[1]
      if (pcs0>=p0star){
        for (d in seq(d_L,d_U,by=d_step)){
          pv.mod <- PV.mod.exact(k,t,n,c,d,p0,delta=delta2)
          if(pv.mod>=p1star){
            N_val <- (k+1)*n
            pv.con <- PV.exact(k,t,n,c,d,p0,delta=delta2)
            pcs1 <- PCS1.exact(k,t,n,c,d,p0,delta=delta2)
            ES1<-ES.1.exact(k,t,n,c,d,p0,delta1=delta2,delta2)$ES1
            ES<-ES.exact(k,t,n,c,d,p0,delta=delta2)
            tmp<-data.frame(k,t,p0,n,c,d,ES1,ES,N_val,pcs0,pcs1,pv.con,pv.mod,delta1,delta2,p0star,p1star)
            print(tmp)
            out<-rbind(out,tmp)
            
            done<-TRUE
            break ## no need to continue in current loop as d will increase and thus E(s) will increase
          }
        }
      }
      else 
        break  ## P(CS0) will only be smaller when increasing n while c is fixed
    }
    #if (done){break} ## after a set of (n,c,d) found, if continue to increase c, the appropriate n will only be larger
  }
  out
}

###################################
#### k=2, t=1, p0=0.2 or 0.5 or 0.7
###################################
#p0=0.2
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=2,t=1,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,n_L=70,n_U=70,c_L=10,c_U=12,d_L=8,d_U=14)  #obtain (70,10,13)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

#p0=0.5
RCW1.fvalue(k=2,t=1,p0star=0.95,p1star=0.8,p0=0.5,delta1=0,delta2=0.2,c_L=11,c_U=20,n_L=75,n_U=100,n_step=1) #RCW1: n=96,c=14
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=2,t=1,p0star=0.95,p1star=0.8,p0=0.5,n_L=96,n_U=96,c_L=14,c_U=15,d_L=12,d_U=18)  #obtain (96,14,14)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

#p0=0.7
RCW1.fvalue(k=2,t=1,p0star=0.95,p1star=0.8,p0=0.7,delta1=0,delta2=0.2,c_L=8,c_U=15,n_L=60,n_U=90,n_step=1) #RCW1: n=78,c=12
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=2,t=1,p0star=0.95,p1star=0.8,p0=0.7,n_L=78,n_U=78,c_L=12,c_U=13,d_L=8,d_U=15)  #obtain (78,12,9)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

###################################
#### k=3, t=1, p0=0.2 or 0.5 or 0.7
###################################
#p0=0.2
RCW1.fvalue(k=3,t=1,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,c_L=10,c_U=13,n_L=70,n_U=90,n_step=1) #RCW1: n=76,c=11
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=3,t=1,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,n_L=76,n_U=76,c_L=11,c_U=12,d_L=10,d_U=15)  #obtain (76,11,15)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

#p0=0.5
RCW1.fvalue(k=3,t=1,p0star=0.95,p1star=0.8,p0=0.5,delta1=0,delta2=0.2,c_L=12,c_U=17,n_L=80,n_U=120,n_step=1) #RCW1: n=108,c=16
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=3,t=1,p0star=0.95,p1star=0.8,p0=0.5,delta1=0,delta2=0.2,n_L=108,n_U=108,c_L=16,c_U=18,d_L=14,d_U=22)  #obtain (108,16,15)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

#p0=0.7
RCW1.fvalue(k=3,t=1,p0star=0.95,p1star=0.8,p0=0.7,delta1=0,delta2=0.2,c_L=10,c_U=13,n_L=70,n_U=90,n_step=1) #RCW1: n=84,c=13
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=3,t=1,p0star=0.95,p1star=0.8,p0=0.7,delta1=0,delta2=0.2,n_L=84,n_U=84,c_L=13,c_U=14,d_L=7,d_U=17)  #obtain (84,13,9)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

###################################
#### k=4, t=1, p0=0.2 or 0.5 or 0.7
###################################
#p0=0.2
RCW1.fvalue(k=4,t=1,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,c_L=10,c_U=13,n_L=70,n_U=90,n_step=1) #RCW1: n=82,c=12
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=4,t=1,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,n_L=82,n_U=82,c_L=12,c_U=12,d_L=10,d_U=17)  #obtain (82,12,16)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

#p0=0.5
RCW1.fvalue(k=4,t=1,p0star=0.95,p1star=0.8,p0=0.5,delta1=0,delta2=0.2,c_L=13,c_U=17,n_L=85,n_U=120,n_step=1) #RCW1: n=113,c=17
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=4,t=1,p0star=0.95,p1star=0.8,p0=0.5,delta1=0,delta2=0.2,n_L=113,n_U=113,c_L=17,c_U=18,d_L=17,d_U=22)  #obtain (113,17,19)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

#p0=0.7
RCW1.fvalue(k=4,t=1,p0star=0.95,p1star=0.8,p0=0.7,delta1=0,delta2=0.2,c_L=10,c_U=15,n_L=80,n_U=100,n_step=1) #RCW1: n=90,c=14
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=4,t=1,p0star=0.95,p1star=0.8,p0=0.7,delta1=0,delta2=0.2,n_L=90,n_U=90,c_L=14,c_U=15,d_L=9,d_U=16)  #obtain (90,14,10)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

#####################
#### k=3, t=2, p0=0.1
#####################
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=3,t=2,p0star=0.95,p1star=0.8,p0=0.1,n_L=64,n_U=64,c_L=8,c_U=10,d_L=3,d_U=10)  #obtain (64,8,9)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=3,t=2,p0star=0.95,p1star=0.9,p0=0.1,n_L=83,n_U=83,c_L=9,c_U=12,d_L=7,d_U=12)  #obtain (83,9,11)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

#####################
#### k=3, t=2, p0=0.2
#####################
RCW1.fvalue(k=3,t=2,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,c_L=9,c_U=12,n_L=65,n_U=100,n_step=1) #RCW1: n=94,c=12
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=3,t=2,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,n_L=94,n_U=94,c_L=12,c_U=13,d_L=10,d_U=15)  #obtain (94,12,14)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

#####################
#### k=3, t=2, p0=0.3
#####################
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=3,t=2,p0star=0.95,p1star=0.8,p0=0.3,n_L=116,n_U=116,c_L=14,c_U=17,d_L=13,d_U=17)  #obtain (116,15,14)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=3,t=2,p0star=0.95,p1star=0.9,p0=0.3,n_L=148,n_U=148,c_L=17,c_U=20,d_L=20,d_U=25)  #obtain (148,17,23)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

#####################
#### k=3, t=2, p0=0.5
#####################
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=3,t=2,p0star=0.95,p1star=0.8,p0=0.5,n_L=127,n_U=127,c_L=17,c_U=20,d_L=15,d_U=25)  #obtain (127,17,16)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=3,t=2,p0star=0.95,p1star=0.9,p0=0.5,n_L=161,n_U=161,c_L=19,c_U=22,d_L=16,d_U=25)  #obtain (161,19,17)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

#####################
#### k=3, t=2, p0=0.7
#####################
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=3,t=2,p0star=0.95,p1star=0.8,p0=0.7,n_L=98,n_U=98,c_L=14,c_U=15,d_L=7,d_U=15)  #obtain (98,14,9)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=3,t=2,p0star=0.95,p1star=0.9,p0=0.7,n_L=119,n_U=119,c_L=15,c_U=16,d_L=8,d_U=17)  #obtain (119,15,10)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

#####################
#### k=4, t=2, p0=0.1
#####################
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=4,t=2,p0star=0.95,p1star=0.8,p0=0.1,n_L=64,n_U=64,c_L=8,c_U=10,d_L=5,d_U=12)  #obtain (64,8,10)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=4,t=2,p0star=0.95,p1star=0.9,p0=0.1,n_L=83,n_U=83,c_L=9,c_U=10,d_L=9,d_U=15)  #obtain (83,9,12)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

#####################
#### k=4, t=2, p0=0.2
#####################
RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,c_L=12,c_U=17,n_L=83,n_U=128,n_step=1) #RCW1: n=101,c=13
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=4,t=2,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,n_L=101,n_U=101,c_L=13,c_U=13,d_L=10,d_U=15)  #obtain (101,13,13)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

#####################
#### k=4, t=2, p0=0.3
#####################
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=4,t=2,p0star=0.95,p1star=0.8,p0=0.3,n_L=128,n_U=128,c_L=17,c_U=18,d_L=15,d_U=23)  #obtain (128,17,17)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=4,t=2,p0star=0.95,p1star=0.9,p0=0.3,n_L=155,n_U=155,c_L=18,c_U=19,d_L=17,d_U=23)  #obtain (155,18,21)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

#####################
#### k=4, t=2, p0=0.5
#####################
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=4,t=2,p0star=0.95,p1star=0.8,p0=0.5,n_L=139,n_U=139,c_L=19,c_U=20,d_L=15,d_U=22)  #obtain (139,19,19)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=4,t=2,p0star=0.95,p1star=0.9,p0=0.5,n_L=171,n_U=171,c_L=21,c_U=22,d_L=23,d_U=28)  #could not find solution 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

#####################
#### k=4, t=2, p0=0.7
#####################
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=4,t=2,p0star=0.95,p1star=0.8,p0=0.7,n_L=104,n_U=104,c_L=15,c_U=16,d_L=8,d_U=16)  #obtain (104,15,10)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

#####################
#### k=5, t=2
#####################
RCW1.fvalue(k=5,t=2,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,c_L=12,c_U=17,n_L=100,n_U=145,n_step=1) #RCW1: n=107,c=14
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=5,t=2,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,n_L=107,n_U=107,c_L=14,c_U=15,d_L=10,d_U=16)  #obtain (107,14,15)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

RCW1.fvalue(k=5,t=2,p0star=0.95,p1star=0.8,p0=0.5,delta1=0,delta2=0.2,c_L=18,c_U=21,n_L=140,n_U=160,n_step=1) #RCW1: n=145,c=20
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=5,t=2,p0star=0.95,p1star=0.8,p0=0.5,delta1=0,delta2=0.2,n_L=145,n_U=145,c_L=20,c_U=20,d_L=15,d_U=25)  #obtain (145,20,20)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

RCW1.fvalue(k=5,t=2,p0star=0.95,p1star=0.8,p0=0.7,delta1=0,delta2=0.2,c_L=14,c_U=17,n_L=90,n_U=120,n_step=1) #RCW1: n=110,c=16
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=5,t=2,p0star=0.95,p1star=0.8,p0=0.7,delta1=0,delta2=0.2,n_L=110,n_U=110,c_L=16,c_U=16,d_L=10,d_U=16)  #obtain (110,16,11)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")


#####################
#### k=5, t=3
#####################
RCW1.fvalue(k=5,t=3,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,c_L=12,c_U=17,n_L=100,n_U=145,n_step=1) #RCW1: n=120,c=15
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=5,t=3,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.2,n_L=120,n_U=120,c_L=15,c_U=16,d_L=13,d_U=20)  #obtain (120,15,16)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

RCW1.fvalue(k=5,t=3,p0star=0.95,p1star=0.8,p0=0.5,delta1=0,delta2=0.2,c_L=20,c_U=21,n_L=155,n_U=160,n_step=1) #RCW1: n=160,c=21
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=5,t=3,p0star=0.95,p1star=0.8,p0=0.5,delta1=0,delta2=0.2,n_L=160,n_U=160,c_L=21,c_U=22,d_L=15,d_U=25)  #obtain (160,21,16)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

RCW1.fvalue(k=5,t=3,p0star=0.95,p1star=0.8,p0=0.7,delta1=0,delta2=0.2,c_L=14,c_U=17,n_L=120,n_U=135,n_step=1) #RCW1: n=121,c=17
tmp<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
out<-RCW4.fvalue(k=5,t=3,p0star=0.95,p1star=0.8,p0=0.7,delta1=0,delta2=0.2,n_L=121,n_U=121,c_L=17,c_U=18,d_L=8,d_U=20)  #obtain (120,15,16)* 
out<-out %>% arrange(ES1)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW4_mod_optimal_designs.csv")

#######################################################################################
#### For optimal designs found above, calculate E_1(s) at slippage configuration
#### Slippage conf: p1=...=p_{k-t}=p0+delta1
####                p_{k-t+1}=...=p_k=p0+delta2
#######################################################################################
opt.dsgn<-read.csv(file="R_CW4_mod_optimal_designs.csv")[,-1]
tmp<-NULL
tmp<-read.csv(file="E1S_SPC.csv")[,-1]
for (i in 29:nrow(opt.dsgn)){ 
  param<-opt.dsgn[i,]
  if (param$k<5){ES1<-ES.1.exact(k=param$k,t=param$t,n=param$n,c=param$c,d=param$d,p0=param$p0,delta1=0,delta2=0.2)$ES1}
  else {# k>=5 cannot be calculated in reasonable time, so it is calcualted via simulation
    ES1<-subset.sim.onestage.calc(nsim=100000,k=param$k,t=param$t,g=param$t,n=param$n,c=param$c,d=param$d,p0=param$p0,
                                  p_lower=param$p0,p_higher=param$p0+0.2)$ES.1
    }
  out<-data.frame(k=param$k,t=param$t,n=param$n,c=param$c,d=param$d,p0=param$p0,delta1=0,delta2=0.2,ES1)
  print(out)
  tmp<-rbind(tmp,out)
}
write.csv(tmp,file="E1S_SPC.csv")

###################################################################
#### Example Sec 5.6
#### For k6t3, n=20, delta1=-0.3, delta2=-0.1
#### i) find smallest c to satisfy P0*=0.75
#### ii) for the c found, find smallest d to satisfy P1*=0.75
#### iii) calculate ES1 at EPC (p1=...=pk=p0+delta2) 
####      and at SPC (p1=...=p_{k-t}=p0+delta1, other pi=p0+delta2)
###################################################################
out<-RCW4.fvalue(k=6,t=3,p0star=0.7,p1star=0.7,p0=0.9,delta1=-0.3,delta2=0,
                 n_L=25,n_U=25,c_L=-3,c_U=3,d_L=1,d_U=10)  #obtain (25,-2,4)* 
ES1<-subset.sim.onestage.calc(nsim=100000,k=6,t=3,g=3,n=25,c=-2,d=3,p0=param$p0,
                              p_lower=0.9,p_higher=0.6)$ES.1

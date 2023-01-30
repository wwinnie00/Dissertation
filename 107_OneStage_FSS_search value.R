###################################################################
### One-Stage FSS design allowing less than t
### If X_{k-t+1} >= X_0+c, select the t largest X_i
### Otherwise, select only X_i >= X_0+c
##########################################################
### Original version: 2022/07/29

source(file="~/Dropbox/AAAA Paper/Codes/Ch1_One stage FSS/105_OneStage_FSS_PCS_exact.R")
setwd("~/Dropbox/AAAA Paper/Codes/Ch1_One stage FSS/Outputs/107_one stage FSS optimal designs_new")

# Search values, P(CS0) and P(CS1) requirements

RCW1.fvalue <- function(k,t,p0star,p1star,p0, delta1=0.05 ,delta2=0.2,
                          c_L,c_U,c_step=1,n_L,n_U,n_step=1){
  out<-NULL
  done<-FALSE
  for (c in seq(c_L,c_U,by=c_step)){
    for (n in seq(n_L,n_U,by=n_step)){
      if (t==1){pcs1 <- dun.pcs1.t1(k,n,c,p0,delta1,delta2)}
        else if (t==2) {pcs1 <- dun.pcs1.t2(k,n,c,p0,delta1,delta2)}
        else if (t==3) {pcs1 <- dun.pcs1.t3(k,n,c,p0,delta1,delta2)}
        else {stop('Undefined value of t!')}
      
        if (pcs1>=p1star){
          pcs0 <- RCW1.pcs.seq(k,t,n,c,p0,delta1,delta2)[1]
          if(pcs0>=p0star){
            N_val <- (k+1)*n
            out<-data.frame(k,t,p0,delta1,delta2,c,n,N_val,pcs0,pcs1)
            print(out)

            done<-TRUE
            break ## no need to continue in current loop as n will be larger
          }
          else 
            break  ## P(CS0) will only be smaller when increasing n while c is fixed
        }
    }
    if (done){break} ## after a set of (n,c) found, if continue to increase c, the appropriate n will only be larger
  }
  out
}

# #####################
# #### k=3, t=2, p0=0.1
# #####################
# # P^* = 0.80
# RCW1.fvalue(k=3,t=2,p0star = 0.8,p1star = 0.8,p0=0.1,c_L=10,c_U=15,n_L=80,n_U=100,n_step=1)  #obtain (11,83)* N=332
# out<-onestage.design.out(k=3,t=2,n=83,c=11,p0=0.1,delta1=0.05,delta2=0.2)
# tmp<-NULL
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.85
# RCW1.fvalue(k=3,t=2,p0star = 0.85,p1star = 0.85,p0=0.1,c_L=12,c_U=17,n_L=95,n_U=110,n_step=1)  #obtain (13,101)* N=404
# out<-onestage.design.out(k=3,t=2,n=101,c=13,p0=0.1,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.90
# RCW1.fvalue(k=3,t=2,p0star = 0.9,p1star = 0.9,p0=0.1,c_L=15,c_U=20,n_L=125,n_U=150,n_step=1)  #obtain (17,134)* N=536
# out<-onestage.design.out(k=3,t=2,n=134,c=17,p0=0.1,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.80, delta1=0
# RCW1.fvalue(k=3,t=2,p0star = 0.95,p1star = 0.8,p0=0.1,delta1=0,c_L=7,c_U=8,n_L=55,n_U=70,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=64,c=8,p0=0.1,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.85, delta1=0
# RCW1.fvalue(k=3,t=2,p0star = 0.95,p1star = 0.85,p0=0.1,delta1=0,c_L=7,c_U=8,n_L=55,n_U=70,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=69,c=8,p0=0.1,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.9, delta1=0
# RCW1.fvalue(k=3,t=2,p0star = 0.95,p1star = 0.9,p0=0.1,delta1=0,c_L=7,c_U=9,n_L=70,n_U=85,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=83,c=9,p0=0.1,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# ######################
# ##### k=3, t=2, p0=0.2
# ######################
# # P^* = 0.80, delta1=0.05
# RCW1.fvalue(k=3,t=2,p0star = 0.8,p1star = 0.8,p0=0.2,c_L=13,c_U=17,n_L=100,n_U=120,n_step=1) #obtain (15,113)* N=452
# 
# out<-onestage.design.out(k=3,t=2,n=113,c=15,p0=0.2,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.85, delta1=0.05
# RCW1.fvalue(k=3,t=2,p0star = 0.85,p1star = 0.85,p0=0.2,c_L=16,c_U=19,n_L=130,n_U=150,n_step=1) #obtain (19,145)* N=580
# out<-onestage.design.out(k=3,t=2,n=145,c=19,p0=0.2,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.90, delta1=0.05
# RCW1.fvalue(k=3,t=2,p0star = 0.9,p1star = 0.9,p0=0.2,c_L=22,c_U=25,n_L=175,n_U=195,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=187,c=24,p0=0.2,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# ######################
# ##### k=3, t=2, p0=0.3
# ######################
# # P^* = 0.80
# RCW1.fvalue(k=3,t=2,p0star = 0.8,p1star = 0.8,p0=0.3,c_L=15,c_U=18,n_L=120,n_U=140,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=135,c=18,p0=0.3,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.85
# RCW1.fvalue(k=3,t=2,p0star=0.85,p1star=0.85,p0=0.3,c_L=19,c_U=22,n_L=150,n_U=170,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=168,c=22,p0=0.3,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.90
# RCW1.fvalue(k=3,t=2,p0star = 0.9,p1star = 0.9,p0=0.3,c_L=26,c_U=29,n_L=205,n_U=225,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=223,c=29,p0=0.3,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.8, delta1=0
# RCW1.fvalue(k=3,t=2,p0star = 0.95,p1star = 0.8,p0=0.3,delta1=0,c_L=12,c_U=16,n_L=100,n_U=125,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=116,c=15,p0=0.3,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.85, delta1=0
# RCW1.fvalue(k=3,t=2,p0star = 0.95,p1star = 0.85,p0=0.3,delta1=0,c_L=15,c_U=17,n_L=120,n_U=140,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=131,c=16,p0=0.3,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.9, delta1=0
# RCW1.fvalue(k=3,t=2,p0star = 0.95,p1star = 0.9,p0=0.3,delta1=0,c_L=15,c_U=17,n_L=130,n_U=150,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=148,c=17,p0=0.3,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# ######################
# ##### k=3, t=2, p0=0.4
# ######################
# # P^* = 0.80
# RCW1.fvalue(k=3,t=2,p0star = 0.8,p1star = 0.8,p0=0.4,c_L=18,c_U=20,n_L=135,n_U=150,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=142,c=19,p0=0.4,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.85
# RCW1.fvalue(k=3,t=2,p0star = 0.85,p1star = 0.85,p0=0.4,c_L=23,c_U=25,n_L=175,n_U=190,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=181,c=24,p0=0.4,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.90
# RCW1.fvalue(k=3,t=2,p0star = 0.9,p1star = 0.9,p0=0.4,c_L=30,c_U=32,n_L=230,n_U=245,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=243,c=32,p0=0.4,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# ######################
# ##### k=3, t=2, p0=0.5
# ######################
# # P^* = 0.80
# RCW1.fvalue(k=3,t=2,p0star = 0.8,p1star = 0.8,p0=0.5,c_L=17,c_U=20,n_L=135,n_U=150,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=140,c=19,p0=0.5,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.85
# RCW1.fvalue(k=3,t=2,p0star = 0.85,p1star = 0.85,p0=0.5,c_L=22,c_U=25,n_L=175,n_U=190,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=179,c=24,p0=0.5,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.90
# RCW1.fvalue(k=3,t=2,p0star = 0.9,p1star = 0.9,p0=0.5,c_L=29,c_U=32,n_L=230,n_U=250,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=235,c=31,p0=0.5,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.8, delta1=0
# RCW1.fvalue(k=3,t=2,p0star = 0.95,p1star = 0.8,p0=0.5,delta1=0,c_L=15,c_U=18,n_L=120,n_U=135,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=127,c=17,p0=0.5,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.85, delta1=0
# RCW1.fvalue(k=3,t=2,p0star = 0.95,p1star = 0.85,p0=0.5,delta1=0,c_L=15,c_U=19,n_L=130,n_U=150,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=142,c=18,p0=0.5,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.9, delta1=0
# RCW1.fvalue(k=3,t=2,p0star = 0.95,p1star = 0.9,p0=0.5,delta1=0,c_L=17,c_U=19,n_L=145,n_U=165,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=161,c=19,p0=0.5,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# ######################
# ##### k=3, t=2, p0=0.6
# ######################
# # P^* = 0.80
# RCW1.fvalue(k=3,t=2,p0star = 0.8,p1star = 0.8,p0=0.6,c_L=16,c_U=19,n_L=120,n_U=130,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=124,c=17,p0=0.6,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.85
# RCW1.fvalue(k=3,t=2,p0star = 0.85,p1star = 0.85,p0=0.6,c_L=20,c_U=22,n_L=160,n_U=175,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=162,c=22,p0=0.6,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.90
# RCW1.fvalue(k=3,t=2,p0star = 0.9,p1star = 0.9,p0=0.6,c_L=27,c_U=29,n_L=215,n_U=230,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=216,c=29,p0=0.6,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# ######################
# ##### k=3, t=2, p0=0.7
# ######################
# # P^* = 0.80
# RCW1.fvalue(k=3,t=2,p0star = 0.8,p1star = 0.8,p0=0.7,c_L=13,c_U=15,n_L=95,n_U=115,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=98,c=14,p0=0.7,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.85
# RCW1.fvalue(k=3,t=2,p0star = 0.85,p1star = 0.85,p0=0.7,c_L=16,c_U=18,n_L=125,n_U=140,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=128,c=18,p0=0.7,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.90
# RCW1.fvalue(k=3,t=2,p0star = 0.9,p1star = 0.9,p0=0.7,c_L=21,c_U=23,n_L=165,n_U=180,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=167,c=23,p0=0.7,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.8, delta1=0
# RCW1.fvalue(k=3,t=2,p0star = 0.95,p1star = 0.8,p0=0.7,delta1=0,c_L=12,c_U=14,n_L=90,n_U=100,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=98,c=14,p0=0.7,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.85, delta1=0
# RCW1.fvalue(k=3,t=2,p0star = 0.95,p1star = 0.85,p0=0.7,delta1=0,c_L=12,c_U=15,n_L=100,n_U=130,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=110,c=15,p0=0.7,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.9, delta1=0
# RCW1.fvalue(k=3,t=2,p0star = 0.95,p1star = 0.9,p0=0.7,delta1=0,c_L=14,c_U=16,n_L=110,n_U=130,n_step=1)
# out<-onestage.design.out(k=3,t=2,n=119,c=15,p0=0.7,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# ######################
# ##### k=4, t=2, p0=0.1
# ######################
# # P^* = 0.80
# RCW1.fvalue(k=4,t=2,p0star = 0.8,p1star = 0.8,p0=0.1,c_L=11,c_U=14,n_L=85,n_U=105,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=90,c=12,p0=0.1,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.85
# RCW1.fvalue(k=4,t=2,p0star = 0.85,p1star = 0.85,p0=0.1,c_L=14,c_U=17,n_L=110,n_U=130,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=113,c=15,p0=0.1,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.90
RCW1.fvalue(k=4,t=2,p0star = 0.9,p1star = 0.9,p0=0.1,c_L=18,c_U=20,n_L=145,n_U=165,n_step=1)
out<-onestage.design.out(k=4,t=2,n=146,c=19,p0=0.1,delta1=0.05,delta2=0.2)
tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
tmp<-rbind(tmp,out)
write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.8, delta1=0
# RCW1.fvalue(k=4,t=2,p0star = 0.95,p1star = 0.8,p0=0.1,delta1=0,c_L=6,c_U=13,n_L=50,n_U=80,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=64,c=8,p0=0.1,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.85, delta1=0
# RCW1.fvalue(k=4,t=2,p0star = 0.95,p1star = 0.85,p0=0.1,delta1=0,c_L=7,c_U=9,n_L=50,n_U=80,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=76,c=9,p0=0.1,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.9, delta1=0
# RCW1.fvalue(k=4,t=2,p0star = 0.95,p1star = 0.9,p0=0.1,delta1=0,c_L=8,c_U=10,n_L=80,n_U=100,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=83,c=9,p0=0.1,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# ######################
# ##### k=4, t=2, p0=0.2
# ######################
# # P^* = 0.80
# RCW1.fvalue(k=4,t=2,p0star = 0.8,p1star = 0.8,p0=0.2,c_L=16,c_U=18,n_L=125,n_U=140,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=126,c=17,p0=0.2,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.85
# RCW1.fvalue(k=4,t=2,p0star = 0.85,p1star = 0.85,p0=0.2,c_L=20,c_U=22,n_L=155,n_U=170,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=158,c=21,p0=0.2,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.90
# RCW1.fvalue(k=4,t=2,p0star = 0.9,p1star = 0.9,p0=0.2,c_L=25,c_U=27,n_L=200,n_U=215,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=205,c=27,p0=0.2,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# ######################
# ##### k=4, t=2, p0=0.3
# ######################
# # P^* = 0.80
# RCW1.fvalue(k=4,t=2,p0star = 0.8,p1star = 0.8,p0=0.3,c_L=19,c_U=21,n_L=145,n_U=160,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=147,c=20,p0=0.3,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.85
# RCW1.fvalue(k=4,t=2,p0star = 0.85,p1star = 0.85,p0=0.3,c_L=24,c_U=26,n_L=185,n_U=200,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=186,c=25,p0=0.3,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.90
# RCW1.fvalue(k=4,t=2,p0star = 0.9,p1star = 0.9,p0=0.3,c_L=30,c_U=32,n_L=235,n_U=250,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=242,c=32,p0=0.3,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.8, delta1=0
# RCW1.fvalue(k=4,t=2,p0star = 0.95,p1star = 0.8,p0=0.3,delta1=0,c_L=15,c_U=17,n_L=110,n_U=130,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=128,c=17,p0=0.3,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.85, delta1=0
# RCW1.fvalue(k=4,t=2,p0star = 0.95,p1star = 0.85,p0=0.3,delta1=0,c_L=16,c_U=18,n_L=130,n_U=150,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=137,c=17,p0=0.3,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.8, delta1=0
# RCW1.fvalue(k=4,t=2,p0star = 0.95,p1star = 0.9,p0=0.3,delta1=0,c_L=16,c_U=18,n_L=140,n_U=160,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=155,c=18,p0=0.3,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# ######################
# ##### k=4, t=2, p0=0.4
# ######################
# # P^* = 0.80
# RCW1.fvalue(k=4,t=2,p0star = 0.8,p1star = 0.8,p0=0.4,c_L=20,c_U=22,n_L=150,n_U=170,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=154,c=21,p0=0.4,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.85
# RCW1.fvalue(k=4,t=2,p0star = 0.85,p1star = 0.85,p0=0.4,c_L=25,c_U=27,n_L=190,n_U=210,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=193,c=26,p0=0.4,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.90
# RCW1.fvalue(k=4,t=2,p0star = 0.9,p1star = 0.9,p0=0.4,c_L=33,c_U=35,n_L=255,n_U=270,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=256,c=34,p0=0.4,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# ######################
# ##### k=4, t=2, p0=0.5
# ######################
# # P^* = 0.80
# RCW1.fvalue(k=4,t=2,p0star = 0.8,p1star = 0.8,p0=0.5,c_L=20,c_U=22,n_L=150,n_U=165,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=152,c=21,p0=0.5,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.85
# RCW1.fvalue(k=4,t=2,p0star = 0.85,p1star = 0.85,p0=0.5,c_L=25,c_U=27,n_L=190,n_U=210,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=191,c=26,p0=0.5,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.90
# RCW1.fvalue(k=4,t=2,p0star = 0.9,p1star = 0.9,p0=0.5,c_L=32,c_U=34,n_L=245,n_U=260,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=253,c=34,p0=0.5,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.8, delta1=0
# RCW1.fvalue(k=4,t=2,p0star = 0.95,p1star = 0.8,p0=0.5,delta1=0,c_L=17,c_U=19,n_L=120,n_U=140,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=139,c=19,p0=0.5,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.85, delta1=0
# RCW1.fvalue(k=4,t=2,p0star = 0.95,p1star = 0.85,p0=0.5,delta1=0,c_L=18,c_U=20,n_L=140,n_U=160,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=155,c=20,p0=0.5,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.9, delta1=0
# RCW1.fvalue(k=4,t=2,p0star = 0.95,p1star = 0.9,p0=0.5,delta1=0,c_L=19,c_U=21,n_L=155,n_U=180,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=173,c=21,p0=0.5,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# ######################
# ##### k=4, t=2, p0=0.6
# ######################
# # P^* = 0.80
# RCW1.fvalue(k=4,t=2,p0star = 0.8,p1star = 0.8,p0=0.6,c_L=18,c_U=20,n_L=135,n_U=150,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=136,c=19,p0=0.6,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.85
# RCW1.fvalue(k=4,t=2,p0star = 0.85,p1star = 0.85,p0=0.6,c_L=23,c_U=25,n_L=170,n_U=190,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=174,c=24,p0=0.6,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.90
# RCW1.fvalue(k=4,t=2,p0star = 0.9,p1star = 0.9,p0=0.6,c_L=30,c_U=32,n_L=225,n_U=240,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=228,c=31,p0=0.6,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# ######################
# ##### k=4, t=2, p0=0.7
# ######################
# # P^* = 0.80
# RCW1.fvalue(k=4,t=2,p0star = 0.8,p1star = 0.8,p0=0.7,c_L=14,c_U=16,n_L=100,n_U=125,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=104,c=15,p0=0.7,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.85
# RCW1.fvalue(k=4,t=2,p0star = 0.85,p1star = 0.85,p0=0.7,c_L=18,c_U=21,n_L=125,n_U=140,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=134,c=19,p0=0.7,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P^* = 0.90
# RCW1.fvalue(k=4,t=2,p0star = 0.9,p1star = 0.9,p0=0.7,c_L=24,c_U=26,n_L=170,n_U=190,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=179,c=25,p0=0.7,delta1=0.05,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.8, delta1=0
# RCW1.fvalue(k=4,t=2,p0star = 0.95,p1star = 0.8,p0=0.7,delta1=0,c_L=13,c_U=15,n_L=90,n_U=105,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=104,c=15,p0=0.7,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.85, delta1=0
# RCW1.fvalue(k=4,t=2,p0star = 0.95,p1star = 0.85,p0=0.7,delta1=0,c_L=14,c_U=16,n_L=105,n_U=125,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=116,c=16,p0=0.7,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
# # P0^* = 0.95, P1^*=0.9, delta1=0
# RCW1.fvalue(k=4,t=2,p0star = 0.95,p1star = 0.9,p0=0.7,delta1=0,c_L=15,c_U=17,n_L=120,n_U=140,n_step=1)
# out<-onestage.design.out(k=4,t=2,n=131,c=17,p0=0.7,delta1=0,delta2=0.2)
# tmp<-read.csv(file="R_CW1_mod_optimal_designs.csv")[,-1]
# tmp<-rbind(tmp,out)
# write.csv(tmp,file="R_CW1_mod_optimal_designs.csv")
# 
###############################################################################
### Example: Section 2.5, Table 2.3
### k=4, t=2
### p0=0.2, 0.25, 0.3
### P0*=0.95, P1*=0.8 or 0.85 or 0.9
###############################################################################
# RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.5,c_L=3,c_U=9,n_L=15,n_U=20,n_step=1) #n=17,c=6
# RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.85,p0=0.2,delta1=0,delta2=0.5,c_L=3,c_U=9,n_L=15,n_U=20,n_step=1) #n=19,c=6
# RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.9,p0=0.2,delta1=0,delta2=0.5,c_L=3,c_U=9,n_L=15,n_U=25,n_step=1) #n=23,c=7
# 
# RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.8,p0=0.25,delta1=0,delta2=0.5,c_L=3,c_U=9,n_L=15,n_U=20,n_step=1) #n=20,c=7
# RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.85,p0=0.25,delta1=0,delta2=0.5,c_L=3,c_U=9,n_L=20,n_U=30,n_step=1) #n=21,c=7
# RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.9,p0=0.25,delta1=0,delta2=0.5,c_L=3,c_U=9,n_L=20,n_U=25,n_step=1) #n=23,c=7
# 
# RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.8,p0=0.3,delta1=0,delta2=0.5,c_L=3,c_U=9,n_L=15,n_U=20,n_step=1) #n=20,c=7
# RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.85,p0=0.3,delta1=0,delta2=0.5,c_L=3,c_U=9,n_L=20,n_U=30,n_step=1) #n=21,c=7
# RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.9,p0=0.3,delta1=0,delta2=0.5,c_L=3,c_U=9,n_L=20,n_U=30,n_step=1) #n=25,c=8

# # delta2=0.3
# RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.8,p0=0.15,delta1=0,delta2=0.3,c_L=7,c_U=11,n_L=30,n_U=50,n_step=1) #n=41,c=8
# RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.85,p0=0.15,delta1=0,delta2=0.3,c_L=7,c_U=10,n_L=42,n_U=60,n_step=1) #n=44,c=8
# RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.9,p0=0.15,delta1=0,delta2=0.3,c_L=8,c_U=11,n_L=45,n_U=60,n_step=1) #n=52,c=9
# 
# RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.8,p0=0.2,delta1=0,delta2=0.3,c_L=7,c_U=10,n_L=40,n_U=60,n_step=1) #n=46,c=9
# RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.85,p0=0.2,delta1=0,delta2=0.3,c_L=10,c_U=15,n_L=48,n_U=60,n_step=1) #n=53,c=10
# RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.9,p0=0.2,delta1=0,delta2=0.3,c_L=11,c_U=15,n_L=55,n_U=70,n_step=1) #n=62,c=11
# 
# RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.8,p0=0.3,delta1=0,delta2=0.3,c_L=9,c_U=15,n_L=50,n_U=70,n_step=1) #n=55,c=11
# RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.85,p0=0.3,delta1=0,delta2=0.3,c_L=11,c_U=15,n_L=56,n_U=70,n_step=1) #n=63,c=12
# RCW1.fvalue(k=4,t=2,p0star=0.95,p1star=0.9,p0=0.3,delta1=0,delta2=0.3,c_L=12,c_U=15,n_L=65,n_U=80,n_step=1) #n=72,c=13

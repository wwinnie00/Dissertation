###### Exact calculation of P(CS) & E(N) for fixed-sample-size TSE two-stage allowing less than t #######
###### Procedure:
######    Stage1. If X_{[k-t+1],1} >= X_{01}+c1, select the t best populations and use randomization to break a tie, then continue to S2.
######            Otherwise, select only populations with X_{i1}>=X_{01}+c1 and continue to S2.
######            If no population has X_{i1}>=X_{01}+c1, select none and terminate the experiment.
######    Stage2. Select only populations with X_{i1}+X_{i2}>=X_{01}+X_{02}+c2.

source("/Users/Dropbox/310_TwoStage_PCS_EN_exact.R")

library(dplyr)

optimal_twostage <- function(k,pstar,alpha,p0, delta1 ,delta2,
                             c1_L,c1_U,c1_step=1,
                             c2_L,c2_U,c2_step=1,
                             n1_L,n1_U,n1_step=1,
                             n2_L,n2_U,n2_step=1){
  
  m <- 1
  plist <- list()
  outall<-read.csv(file="R_CW3_mod_designs_gridsearch.csv")[,-1] #comment this out if not saving outputs
  out<-NULL
  
  for (c1 in seq(c1_L,c1_U,by=c1_step)){
    for (c2 in seq(c2_L,c2_U,by=c2_step)){
      for (n1 in seq(n1_L,n1_U,by=n1_step)){
        for (n2 in seq(n2_L,n2_U,by=n2_step)){
          myindex <- ifelse(m==1,1,myfunction(n1,n2,c1))
          if (myindex==1)
          {
            typ1err <- 1-pcs0(k,n1,n2,c1,c2,p0)
            if (typ1err>=alpha){break} # type I error will only be larger when increasing n2 while the other paramaters are fixed
            else{
              power <- pcs1(k,n1,n2,c1,c2,p0,delta1,delta2)
              if(power>=pstar){
                plist[[m]] <- c(n1,n2,c1,c2)
                N_vec <- EN(k,n1,n2,c1,p0,delta1,delta2)
                EN.avr <- N_vec$EN.avr
                EN.H0<-N_vec$EN.H0
                EN.LFC<-N_vec$EN.LFC
                Nmax <- N_vec$Nmax
                tau <- N_vec$tau.H0
                tmp<-data.frame(k,t=2,p0,n1,n2,c1,c2,delta1,delta2,power,typ1err,EN.avr,EN.H0,EN.LFC,Nmax,tau)
                print(tmp)
                out<-rbind(out,tmp)
                outall<-rbind(outall,tmp)
                write.csv(outall,file="R_CW3_mod_designs_gridsearch.csv")
                
                m<-m+1
                
                myfunction <- function(n1,n2,c1){ 
                  # return value 0 if (c1,n1,n2) necessarily yields EN no less than any exsisting (c1',c2',n1',n2'), which happens when c1<=c1' and n1>=n1' and n2>=n2' 
                  
                  ss<-1
                  for (l in 1:(m-1)){
                    count1 <- (c1>plist[[l]][3] | n1<plist[[l]][1] | n2<plist[[l]][2])
                    ss <- ss*count1
                  }
                  ss
                }
                break
              } 
            }
          }
        }
      }
    }
  }
  out
}

#################################################################
### k=3, t=2, alpha=0.05, p0=0.1
#################################################################

# P*=0.8
out<-optimal_twostage(k=3,pstar=0.8,alpha=0.05,p0=0.1,delta1=0,delta2=0.2,
                 c1_L=1,c1_U=3,c2_L=8,c2_U=8,n1_L=20,n1_U=36,n1_step=2,n2_L=30,n2_U=60,n2_step=5)
out<-optimal_twostage(k=3,pstar=0.8,alpha=0.05,p0=0.1,delta1=0,delta2=0.2,
                 c1_L=2,c1_U=2,c2_L=8,c2_U=8,n1_L=29,n1_U=36,n1_step=1,n2_L=30,n2_U=49,n2_step=1)
out<-out%>%arrange(EN.avr)
tmp<-NULL
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW3_mod_optimal_designs.csv")

# P*=0.85 
optimal_twostage(k=3,pstar=0.85,alpha=0.05,p0=0.1,delta1=0,delta2=0.2,
                 c1_L=0,c1_U=3,c2_L=8,c2_U=10,n1_L=30,n1_U=50,n1_step=1,n2_L=25,n2_U=45,n2_step=5)
optimal_twostage(k=3,pstar=0.85,alpha=0.05,p0=0.1,delta1=0,delta2=0.2,
                 c1_L=1,c1_U=2,c2_L=8,c2_U=8,n1_L=35,n1_U=45,n1_step=1,n2_L=25,n2_U=35,n2_step=1)
out<-optimal_twostage(k=3,pstar=0.85,alpha=0.05,p0=0.1,delta1=0,delta2=0.2,
                 c1_L=2,c1_U=3,c2_L=8,c2_U=8,n1_L=40,n1_U=50,n1_step=2,n2_L=20,n2_U=60,n2_step=5)
out<-out%>%arrange(EN.avr)
tmp<-read.csv(file="R_CW3_mod_optimal_designs.csv")[,-1]
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW3_mod_optimal_designs.csv")
  
# P*=0.9
out<-optimal_twostage(k=3,pstar=0.9,alpha=0.05,p0=0.1,delta1=0,delta2=0.2,
                 c1_L=0,c1_U=4,c2_L=8,c2_U=10,n1_L=40,n1_U=50,n1_step=1,n2_L=25,n2_U=45,n2_step=5)
out<-optimal_twostage(k=3,pstar=0.9,alpha=0.05,p0=0.1,delta1=0,delta2=0.2,
                 c1_L=1,c1_U=2,c2_L=9,c2_U=9,n1_L=40,n1_U=50,n1_step=1,n2_L=37,n2_U=47,n2_step=1)
out<-optimal_twostage(k=3,pstar=0.9,alpha=0.05,p0=0.1,delta1=0,delta2=0.2,
                      c1_L=2,c1_U=3,c2_L=9,c2_U=9,n1_L=40,n1_U=50,n1_step=1,n2_L=30,n2_U=60,n2_step=5)
out<-optimal_twostage(k=3,pstar=0.9,alpha=0.05,p0=0.1,delta1=0,delta2=0.2,
                      c1_L=3,c1_U=3,c2_L=9,c2_U=9,n1_L=47,n1_U=50,n1_step=1,n2_L=36,n2_U=50,n2_step=1)
out<-optimal_twostage(k=3,pstar=0.9,alpha=0.05,p0=0.1,delta1=0,delta2=0.2,
                      c1_L=2,c1_U=2,c2_L=9,c2_U=9,n1_L=40,n1_U=50,n1_step=1,n2_L=36,n2_U=55,n2_step=1)
out<-out%>%arrange(EN.avr)
tmp<-read.csv(file="R_CW3_mod_optimal_designs.csv")[,-1]
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW3_mod_optimal_designs.csv")  

#################################################################
### k=3, t=2, alpha=0.05, p0=0.3
#################################################################

# P*=0.8
out<-optimal_twostage(k=3,pstar=0.8,alpha=0.05,p0=0.3,delta1=0,delta2=0.2, #did not store outputs
                 c1_L=2,c1_U=4,c2_L=14,c2_U=16,n1_L=60,n1_U=75,n1_step=1,n2_L=40,n2_U=70,n2_step=5)
out<-optimal_twostage(k=3,pstar=0.8,alpha=0.05,p0=0.3,delta1=0,delta2=0.2, 
                 c1_L=2,c1_U=3,c2_L=15,c2_U=15,n1_L=60,n1_U=70,n1_step=1,n2_L=45,n2_U=55,n2_step=1)
out<-out%>%arrange(EN.avr)
tmp<-read.csv(file="R_CW3_mod_optimal_designs.csv")[,-1]
# out<-twostage.out(k=3,n1=64,n2=55,c1=3,c2=15,p0=0.3,delta1=0,delta2=0.2)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW3_mod_optimal_designs.csv")

# P*=0.85
out<-optimal_twostage(k=3,pstar=0.85,alpha=0.05,p0=0.3,delta1=0,delta2=0.2, #did not store outputs
                 c1_L=2,c1_U=4,c2_L=15,c2_U=16,n1_L=68,n1_U=80,n1_step=1,n2_L=65,n2_U=75,n2_step=5)
out<-optimal_twostage(k=3,pstar=0.85,alpha=0.05,p0=0.3,delta1=0,delta2=0.2, 
                 c1_L=3,c1_U=4,c2_L=16,c2_U=16,n1_L=61,n1_U=71,n1_step=1,n2_L=62,n2_U=90,n2_step=1)
out<-out%>%arrange(EN.avr)
tmp<-read.csv(file="R_CW3_mod_optimal_designs.csv")[,-1]
# out<-twostage.out(k=3,n1=62,n2=77,c1=3,c2=16,p0=0.3,delta1=0,delta2=0.2)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW3_mod_optimal_designs.csv")

# P*=0.9
out<-optimal_twostage(k=3,pstar=0.9,alpha=0.05,p0=0.3,delta1=0,delta2=0.2, #did not store outputs
                      c1_L=2,c1_U=4,c2_L=16,c2_U=18,n1_L=74,n1_U=78,n1_step=1,n2_L=80,n2_U=90,n2_step=5)
out<-optimal_twostage(k=3,pstar=0.9,alpha=0.05,p0=0.3,delta1=0,delta2=0.2, 
                      c1_L=3,c1_U=3,c2_L=17,c2_U=17,n1_L=73,n1_U=75,n1_step=1,n2_L=75,n2_U=100,n2_step=1)
out<-out%>%arrange(EN.avr)
tmp<-read.csv(file="R_CW3_mod_optimal_designs.csv")[,-1]
# out<-twostage.out(k=3,n1=74,n2=82,c1=3,c2=17,p0=0.3,delta1=0,delta2=0.2)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW3_mod_optimal_designs.csv")

#################################################################
### k=3, t=2, alpha=0.05, p0=0.5
#################################################################

# P*=0.8
out<-optimal_twostage(k=3,pstar=0.8,alpha=0.05,p0=0.5,delta1=0,delta2=0.2, #did not store outputs
                      c1_L=3,c1_U=4,c2_L=16,c2_U=18,n1_L=63,n1_U=75,n1_step=1,n2_L=55,n2_U=80,n2_step=5)
out<-optimal_twostage(k=3,pstar=0.8,alpha=0.05,p0=0.5,delta1=0,delta2=0.2, 
                      c1_L=3,c1_U=3,c2_L=17,c2_U=17,n1_L=56,n1_U=67,n1_step=1,n2_L=61,n2_U=85,n2_step=1)
out<-out%>%arrange(EN.avr)
tmp<-read.csv(file="R_CW3_mod_optimal_designs.csv")[,-1]
# out<-twostage.out(k=3,n1=60,n2=73,c1=3,c2=17,p0=0.5,delta1=0,delta2=0.2)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW3_mod_optimal_designs.csv")

# P*=0.85
out<-optimal_twostage(k=3,pstar=0.85,alpha=0.05,p0=0.5,delta1=0,delta2=0.2, #did not store outputs
                      c1_L=3,c1_U=5,c2_L=16,c2_U=18,n1_L=64,n1_U=80,n1_step=2,n2_L=60,n2_U=100,n2_step=5)
out<-optimal_twostage(k=3,pstar=0.85,alpha=0.05,p0=0.5,delta1=0,delta2=0.2, 
                      c1_L=3,c1_U=4,c2_L=18,c2_U=18,n1_L=66,n1_U=71,n1_step=1,n2_L=75,n2_U=95,n2_step=1)
out<-out%>%arrange(EN.avr)
tmp<-read.csv(file="R_CW3_mod_optimal_designs.csv")[,-1]
# out<-twostage.out(k=3,n1=67,n2=85,c1=4,c2=18,p0=0.5,delta1=0,delta2=0.2)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW3_mod_optimal_designs.csv")

# P*=0.9
out<-optimal_twostage(k=3,pstar=0.9,alpha=0.05,p0=0.5,delta1=0,delta2=0.2, #did not store outputs
                      c1_L=3,c1_U=4,c2_L=19,c2_U=20,n1_L=80,n1_U=100,n1_step=2,n2_L=80,n2_U=120,n2_step=5)
out<-optimal_twostage(k=3,pstar=0.9,alpha=0.05,p0=0.5,delta1=0,delta2=0.2, 
                      c1_L=4,c1_U=4,c2_L=20,c2_U=20,n1_L=75,n1_U=83,n1_step=1,n2_L=90,n2_U=110,n2_step=1)
out<-optimal_twostage(k=3,pstar=0.9,alpha=0.05,p0=0.5,delta1=0,delta2=0.2, 
                      c1_L=3,c1_U=3,c2_L=19,c2_U=19,n1_L=75,n1_U=86,n1_step=1,n2_L=60,n2_U=95,n2_step=1)
out<-out%>%arrange(EN.avr)
tmp<-read.csv(file="R_CW3_mod_optimal_designs.csv")[,-1]
# out<-twostage.out(k=3,n1=78,n2=88,c1=3,c2=19,p0=0.5,delta1=0,delta2=0.2) #EN=508.21
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW3_mod_optimal_designs.csv")

#################################################################
### k=3, t=2, alpha=0.05, p0=0.7
#################################################################

# P*=0.8
out<-optimal_twostage(k=3,pstar=0.8,alpha=0.05,p0=0.7,delta1=0,delta2=0.2, #did not store outputs
                      c1_L=2,c1_U=4,c2_L=13,c2_U=15,n1_L=40,n1_U=60,n1_step=2,n2_L=35,n2_U=55,n2_step=5)
out<-optimal_twostage(k=3,pstar=0.8,alpha=0.05,p0=0.7,delta1=0,delta2=0.2, 
                      c1_L=3,c1_U=4,c2_L=14,c2_U=14,n1_L=35,n1_U=50,n1_step=1,n2_L=55,n2_U=75,n2_step=1)
out<-out%>%arrange(EN.avr)
tmp<-read.csv(file="R_CW3_mod_optimal_designs.csv")[,-1]
# out<-twostage.out(k=3,n1=38,n2=70,c1=3,c2=14,p0=0.7,delta1=0,delta2=0.2)
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW3_mod_optimal_designs.csv")

# P*=0.85
out<-optimal_twostage(k=3,pstar=0.85,alpha=0.05,p0=0.7,delta1=0,delta2=0.2, #did not store outputs
                      c1_L=3,c1_U=5,c2_L=14,c2_U=16,n1_L=40,n1_U=50,n1_step=2,n2_L=70,n2_U=100,n2_step=5)
out<-optimal_twostage(k=3,pstar=0.85,alpha=0.05,p0=0.7,delta1=0,delta2=0.2, 
                      c1_L=3,c1_U=3,c2_L=15,c2_U=15,n1_L=35,n1_U=45,n1_step=1,n2_L=70,n2_U=100,n2_step=1)
out<-out%>%arrange(EN.avr)
tmp<-read.csv(file="R_CW3_mod_optimal_designs.csv")[,-1]
tmp<-rbind(tmp,out[1,])
# out<-twostage.out(k=3,n1=44,n2=74,c1=3,c2=15,p0=0.7,delta1=0,delta2=0.2)
write.csv(tmp,file="R_CW3_mod_optimal_designs.csv")

# P*=0.9
out<-optimal_twostage(k=3,pstar=0.9,alpha=0.05,p0=0.7,delta1=0,delta2=0.2, #did not store outputs
                      c1_L=3,c1_U=4,c2_L=15,c2_U=16,n1_L=45,n1_U=55,n1_step=2,n2_L=80,n2_U=120,n2_step=5)
out<-optimal_twostage(k=3,pstar=0.9,alpha=0.05,p0=0.7,delta1=0,delta2=0.2, # did not store outputs
                      c1_L=3,c1_U=3,c2_L=15,c2_U=15,n1_L=46,n1_U=60,n1_step=1,n2_L=60,n2_U=90,n2_step=1)
out%>%arrange(EN.avr)
tmp<-read.csv(file="R_CW3_mod_optimal_designs.csv")[,-1]
out<-twostage.out(k=3,n1=51,n2=75,c1=3,c2=15,p0=0.7,delta1=0,delta2=0.2)
tmp<-rbind(tmp,out)
write.csv(tmp,file="R_CW3_mod_optimal_designs.csv")

#################################################################
### k=4, t=2, alpha=0.05, p0=0.3
#################################################################

# P*=0.8
out<-optimal_twostage(k=4,pstar=0.8,alpha=0.05,p0=0.3,delta1=0,delta2=0.2, #did not store outputs
                      c1_L=2,c1_U=4,c2_L=17,c2_U=18,n1_L=55,n1_U=75,n1_step=1,n2_L=60,n2_U=90,n2_step=5)
out<-optimal_twostage(k=4,pstar=0.8,alpha=0.05,p0=0.3,delta1=0,delta2=0.2, #did not store outputs
                      c1_L=3,c1_U=4,c2_L=17,c2_U=17,n1_L=55,n1_U=65,n1_step=1,n2_L=70,n2_U=100,n2_step=1)
out<-out%>%arrange(EN.avr)
tmp<-read.csv(file="R_CW3_mod_optimal_designs.csv")[,-1]
out<-twostage.out(k=4,n1=55,n2=86,c1=3,c2=17,p0=0.3,delta1=0,delta2=0.2) #EN=466.37
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW3_mod_optimal_designs.csv")

# P*=0.85
out<-optimal_twostage(k=4,pstar=0.85,alpha=0.05,p0=0.3,delta1=0,delta2=0.2, #did not store outputs
                      c1_L=3,c1_U=4,c2_L=17,c2_U=18,n1_L=56,n1_U=70,n1_step=2,n2_L=70,n2_U=100,n2_step=5)
out<-optimal_twostage(k=4,pstar=0.85,alpha=0.05,p0=0.3,delta1=0,delta2=0.2, #did not store outputs
                      c1_L=3,c1_U=4,c2_L=18,c2_U=18,n1_L=58,n1_U=65,n1_step=1,n2_L=80,n2_U=105,n2_step=1)
out<-out%>%arrange(EN.avr)
tmp<-read.csv(file="R_CW3_mod_optimal_designs.csv")[,-1]
out<-twostage.out(k=4,n1=63,n2=92,c1=3,c2=18,p0=0.3,delta1=0,delta2=0.2) #EN=523.27
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW3_mod_optimal_designs.csv")


#################################################################
### k=4, t=2, alpha=0.05, p0=0.5
#################################################################

# P*=0.8
out<-optimal_twostage(k=4,pstar=0.8,alpha=0.05,p0=0.5,delta1=0,delta2=0.2, #did not store outputs
                      c1_L=2,c1_U=4,c2_L=19,c2_U=20,n1_L=55,n1_U=75,n1_step=1,n2_L=55,n2_U=90,n2_step=5)
out<-optimal_twostage(k=4,pstar=0.8,alpha=0.05,p0=0.5,delta1=0,delta2=0.2, #did not store outputs
                      c1_L=3,c1_U=4,c2_L=19,c2_U=19,n1_L=58,n1_U=68,n1_step=1,n2_L=76,n2_U=100,n2_step=1)
out<-out%>%arrange(EN.avr)
tmp<-read.csv(file="R_CW3_mod_optimal_designs.csv")[,-1]
out<-twostage.out(k=4,n1=60,n2=93,c1=4,c2=19,p0=0.5,delta1=0,delta2=0.2) #EN=498.0951
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW3_mod_optimal_designs.csv")

# P*=0.85
out<-optimal_twostage(k=4,pstar=0.85,alpha=0.05,p0=0.5,delta1=0,delta2=0.2, #did not store outputs
                      c1_L=4,c1_U=4,c2_L=20,c2_U=21,n1_L=60,n1_U=74,n1_step=2,n2_L=90,n2_U=120,n2_step=5)
out<-optimal_twostage(k=4,pstar=0.85,alpha=0.05,p0=0.5,delta1=0,delta2=0.2, #did not store outputs
                      c1_L=4,c1_U=4,c2_L=20,c2_U=20,n1_L=65,n1_U=70,n1_step=1,n2_L=90,n2_U=105,n2_step=1)
out<-out%>%arrange(EN.avr)
tmp<-read.csv(file="R_CW3_mod_optimal_designs.csv")[,-1]
out<-twostage.out(k=4,n1=,n2=,c1=4,c2=,p0=0.5,delta1=0,delta2=0.2) #EN=498.0951
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW3_mod_optimal_designs.csv")

# P*=0.9
out<-optimal_twostage(k=4,pstar=0.9,alpha=0.05,p0=0.5,delta1=0,delta2=0.2, 
                      c1_L=4,c1_U=5,c2_L=21,c2_U=22,n1_L=76,n1_U=100,n1_step=2,n2_L=100,n2_U=130,n2_step=5)
out<-optimal_twostage(k=4,pstar=0.9,alpha=0.05,p0=0.5,delta1=0,delta2=0.2, 
                      c1_L=4,c1_U=4,c2_L=21,c2_U=21,n1_L=75,n1_U=81,n1_step=1,n2_L=98,n2_U=120,n2_step=1)
out<-out%>%arrange(EN.avr)
tmp<-read.csv(file="R_CW3_mod_optimal_designs.csv")[,-1]
#out<-twostage.out(k=4,n1=75,n2=115,c1=4,c2=21,p0=0.5,delta1=0,delta2=0.2) 
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW3_mod_optimal_designs.csv")

#################################################################
### k=4, t=2, alpha=0.05, p0=0.7
#################################################################

# P*=0.8
out<-optimal_twostage(k=4,pstar=0.8,alpha=0.05,p0=0.7,delta1=0,delta2=0.2, 
                      c1_L=2,c1_U=4,c2_L=15,c2_U=16,n1_L=40,n1_U=60,n1_step=2,n2_L=60,n2_U=90,n2_step=5)
out<-optimal_twostage(k=4,pstar=0.8,alpha=0.05,p0=0.7,delta1=0,delta2=0.2, 
                      c1_L=3,c1_U=3,c2_L=15,c2_U=15,n1_L=32,n1_U=42,n1_step=1,n2_L=65,n2_U=85,n2_step=1)
out<-out%>%arrange(EN.avr)
tmp<-read.csv(file="R_CW3_mod_optimal_designs.csv")[,-1]
#out<-twostage.out(k=4,n1=38,n2=79,c1=3,c2=15,p0=0.7,delta1=0,delta2=0.2) 
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW3_mod_optimal_designs.csv")

# P*=0.85
out<-optimal_twostage(k=4,pstar=0.85,alpha=0.05,p0=0.7,delta1=0,delta2=0.2, 
                      c1_L=3,c1_U=4,c2_L=16,c2_U=17,n1_L=39,n1_U=49,n1_step=2,n2_L=75,n2_U=95,n2_step=5)
out<-optimal_twostage(k=4,pstar=0.85,alpha=0.05,p0=0.7,delta1=0,delta2=0.2, 
                      c1_L=3,c1_U=4,c2_L=16,c2_U=16,n1_L=39,n1_U=49,n1_step=1,n2_L=77,n2_U=100,n2_step=1)
out<-out%>%arrange(EN.avr)
tmp<-read.csv(file="R_CW3_mod_optimal_designs.csv")[,-1]
#out<-twostage.out(k=4,n1=44,n2=83,c1=3,c2=16,p0=0.7,delta1=0,delta2=0.2) 
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW3_mod_optimal_designs.csv")

# P*=0.9
out<-optimal_twostage(k=4,pstar=0.9,alpha=0.05,p0=0.7,delta1=0,delta2=0.2, 
                      c1_L=3,c1_U=4,c2_L=17,c2_U=18,n1_L=40,n1_U=56,n1_step=2,n2_L=75,n2_U=105,n2_step=5)
out<-optimal_twostage(k=4,pstar=0.9,alpha=0.05,p0=0.7,delta1=0,delta2=0.2, 
                      c1_L=3,c1_U=3,c2_L=17,c2_U=17,n1_L=47,n1_U=55,n1_step=1,n2_L=80,n2_U=105,n2_step=1)
out<-out%>%arrange(EN.avr)
tmp<-read.csv(file="R_CW3_mod_optimal_designs.csv")[,-1]
#out<-twostage.out(k=4,n1=50,n2=92,c1=3,c2=17,p0=0.7,delta1=0,delta2=0.2) 
tmp<-rbind(tmp,out[1,])
write.csv(tmp,file="R_CW3_mod_optimal_designs.csv")

################################################
#### Example: Section 4.5
################################################

# p0=0.15
out<-optimal_twostage(k=4,pstar=0.8,alpha=0.05,p0=0.15,delta1=0,delta2=0.3,
                      c1_L=-2,c1_U=2,c2_L=7,c2_U=9,n1_L=10,n1_U=30,n1_step=2,n2_L=10,n2_U=35,n2_step=5)
out<-optimal_twostage(k=4,pstar=0.8,alpha=0.05,p0=0.15,delta1=0,delta2=0.3,
                      c1_L=0,c1_U=2,c2_L=8,c2_U=9,n1_L=14,n1_U=26,n1_step=1,n2_L=15,n2_U=32,n2_step=1)
out%>%arrange(EN.avr)
twostage.out(k=4,n1=19,n2=28,c1=2,c2=8,p0=0.15,delta1=0,delta2=0.3) #EN=154.2475, Nmax=179

# p0=0.2
out<-optimal_twostage(k=4,pstar=0.8,alpha=0.05,p0=0.2,delta1=0,delta2=0.3,
                      c1_L=1,c1_U=3,c2_L=9,c2_U=10,n1_L=20,n1_U=30,n1_step=1,n2_L=25,n2_U=40,n2_step=1)
out<-optimal_twostage(k=4,pstar=0.8,alpha=0.05,p0=0.2,delta1=0,delta2=0.3,
                      c1_L=1,c1_U=2,c2_L=9,c2_U=9,n1_L=22,n1_U=27,n1_step=1,n2_L=20,n2_U=35,n2_step=1)
out%>%arrange(EN.avr)
twostage.out(k=4,n1=24,n2=25,c1=2,c2=9,p0=0.2,delta1=0,delta2=0.3) #EN=175.7028, Nmax=195

# p0=0.3
out<-optimal_twostage(k=4,pstar=0.8,alpha=0.05,p0=0.3,delta1=0,delta2=0.3,
                      c1_L=1,c1_U=3,c2_L=11,c2_U=12,n1_L=26,n1_U=36,n1_step=1,n2_L=20,n2_U=40,n2_step=5)
out<-optimal_twostage(k=4,pstar=0.8,alpha=0.05,p0=0.3,delta1=0,delta2=0.3,
                      c1_L=1,c1_U=3,c2_L=11,c2_U=11,n1_L=25,n1_U=33,n1_step=1,n2_L=20,n2_U=40,n2_step=1)
out%>%arrange(EN.avr)
twostage.out(k=4,n1=29,n2=29,c1=3,c2=11,p0=0.3,delta1=0,delta2=0.3) #EN=205.7909, Nmax=232

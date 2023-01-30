#########################################################################################
## Curtailed version of modified one-stage design, allows <t populations to be selected
## Vector-at-a-time sampling


source(file="~/Dropbox/107_OneStage_FSS_search value.R")

# Calculate power, size and EN
RCW2_pcs_EN_sim <- function(k,t,nsim, n, c, p0, delta1, delta2,g){
### g: true # of populations better than control 
  
  if (t!=2 && !g%in%c(0,t) ){
    stop ("Specified number of best populations to be selected is not legitimate for selected value of t. Please choose g=0 or t")
  }
  
  pb = txtProgressBar(min = 0, max = nsim, style = 3) 
    
  count<-matrix(0,nsim,3) #count[,1]=1 if control is selected, count[,2]=1 if g<t good are selected, count[,3]=1 if t best are selected

  #trytest1 <- trytest2 <- trytest3 <-trytest4<-trytest5<-trytest6<-trytest7<-count
  set.seed(752)
  Nval <- vector()
  for (l in 1:nsim){

    grp <- c(1:(k+1)) # index for each group, control=1, trt from 2 to k+1
    n.obs <- 0
    # define numbers of successes & trials
    y.vec <- rep(0,k+1)
    
    p.vec <- c(p0,rep(p0+delta1,k-g),rep(p0+delta2,g))
    k.tmp <- k
    
    for (M in 1:n){ #M is the number of vectors observed
      if (1 %in% grp){ #control is not eliminated from sampling
        M0<-M #M0 is the # of obs. from control
      }
      b <- vector()
      n.obs <- n.obs+length(grp)
      
      # generate binary data
      for (i in grp){
        b[i] <- rbinom(1,1,p.vec[i])
        y.vec[i] <- y.vec[i]+b[i]
      }
      
      y0 <- y.vec[1]
      
      y.trt <- y.vec[-1]
      
      # rank y values from remaining treatments, same y values at a tie are first come first assign
      sort.data <- data.frame(trt=grp[grp!=1]-1,y=y.vec[grp[grp!=1]],rank=rank(y.vec[grp[grp!=1]],ties.method = "first"))
      
      # tth largest number of successes
      y.t <- sort.data$y[sort.data$rank==k.tmp-t+1] 
      #print(c(y0,y.t))
      
      #(t-1)th largest number of successes
      y.g <- sort.data$y[sort.data$rank==k.tmp-t+2] 
      
      #### Termination Rule (stopping the experiment) 
      # Accept H0: no pop. can possibly exceed y0+c in the end
      if (y0+c>max(y.trt)+n-M){ #y0+c>y_[k]+n-M
        trt.s <- NULL
        count[l,1]<-1
        break
      }
      
      # at least t pop. are better than control
      else if (y.t>=y0+c+n-M0){ 
        if (k.tmp==t){ # only t treatments left
          trt.s<-sort.data$trt[sort.data$rank>=k.tmp-t+1] #select all t treatments
          break
        }
        else{ # more than t remaining treatments
          
          if (M!=n){
            if (y.t>=sort.data$y[sort.data$rank==k.tmp-t]+n-M){# bad trt cannot exceed good trt, strong curtailment applied
              #trytest3=trytest3+1
              trt.s <- sort.data$trt[sort.data$rank>=k.tmp-t+1] # select the t largest treatments
              break
            }
            else{ grp<-grp[-1] } # drop the control and continue sampling from the remaining trts
          }
          else{ # M = n
            if (y.t > sort.data$y[sort.data$rank==k.tmp-t]){ # no tie
              #trytest5=trytest5+1
              trt.s <- sort.data$trt[sort.data$rank>=k.tmp-t+1] # select the t largest treatments
              break
            }
            else{# there is a tie
              #trytest6=trytest6+1
              trt.tmp <- sort.data$trt[sort.data$y>y.t]
              grp.tie <- sort.data$trt[sort.data$y==y.t]
              random.ind <- sample(grp.tie,t-length(trt.tmp))
              trt.s <- c(trt.tmp,random.ind)
              break
            }
          }
        }
      }
      
      # (t-1) pop. are better than the control and the other (k-t+1) pop. are worse than the control
      else if (y.g>=y0+c+n-M0 && y.t+n-M<y0+c){
        trt.s <- c(sort.data$trt[sort.data$y>=y.g])
        break
      }
      ### end of Termination Rule
      
      ### Stop sampling from bad treatments
      if (k.tmp>t){
        # Drop the inferior treatments
        for (j in 1:(k.tmp-t)){
          # compare the t largest with any other treatment
          if (y.t>sort.data$y[sort.data$rank==j]+n-M){ # t largest better than the jth smallest treatment
            #trytest2<-trytest2+1
            trt.omit <- sort.data$trt[sort.data$rank==j] # drop the jth smallest treatment
            grp<-grp[grp!=(trt.omit+1)] #grp exludes eliminated trts
          }
        }
        k.tmp <- length(grp[grp!=1])
      }
    }
  
    truebest <- (k-g+1):k
    
    if (g==t) {count[l,3]<-ifelse(all(trt.s%in%truebest) & all(truebest%in%trt.s),1,0)}
      
    if (g<t && g>0) {count[l,2]<-ifelse(all(trt.s%in%c(0,truebest)),1,0)}
    
    Nval[l] <- n.obs
    setTxtProgressBar(pb,l)
  }
  
  
  if (t==2){# for t=2, print P(CS00) for control only, P(CS01) for selecting less than t, and P(CS1)
    pcs<-sum(count[,g+1])/nsim
  }
  else{# for t=1 or >2, print P(CS00) and P(CS1), skip P(CS01)
    pcs<-ifelse(g==0,sum(count[,1])/nsim,sum(count[,3])/nsim)
  }
  
  output <- data.frame(pcs=pcs, EN=mean(Nval),
                 #trytest1=trytest1,trytest2=trytest2,trytest3=trytest3,trytest4=trytest4,trytest5=trytest5,trytest6=trytest6,trytest7=trytest7,k=k,M=M
                 se=sd(Nval)/sqrt(nsim))
  close(pb)
  output
}

 RCW2_pcs_EN_sim(k=3,t=2,nsim=10000, n=83, c=11, p0=0.1, delta1=0.05, delta2=0.2,g=0)
# RCW2_pcs_EN_sim(k=3,t=2,nsim=10000, n=83, c=11, p0=0.1, delta1=0.05, delta2=0.2,g=1)
# RCW2_pcs_EN_sim(k=3,t=2,nsim=10000, n=83, c=11, p0=0.1, delta1=0.05, delta2=0.2,g=2)

# Calculate average EN, power and size
RCW2_sim <- function(nsim, k,t,n, c, p0, delta1, delta2){
  result.tbest <- RCW2_pcs_EN_sim(nsim=nsim,k=k,t=t, n=n, c=c, p0=p0, delta1=delta1, delta2=delta2,g=t)
  result.nobest <- RCW2_pcs_EN_sim(nsim=nsim,k=k,t=t, n=n, c=c, p0=p0,delta1=delta1, delta2=delta2,g=0)
  
  # calculate E(N)=1/2*(E(N|g=t)+E(N|g=0))
  EN.avr<-1/2*(result.tbest$EN + result.nobest$EN)
  EN.avr.se<-sqrt(1/4*((result.tbest$se)^2+(result.nobest$se)^2))
  
  ## For t=2, calculate E(N)=1/3*(E(N|g=t)+E(N|g=1)+E(N|g=0))
  # result.gbest <- RCW2_pcs_EN_sim(nsim=nsim,k=k,t=t, n=n, c=c, p0=p0, delta1=delta1, delta2=delta2,g=t-1)
  # EN.avr<-1/3*(result.tbest$EN + result.gbest$EN + result.nobest$EN)
  # EN.avr.se<-sqrt(1/9*((result.tbest$se)^2+(result.gbest$se)^2+(result.nobest$se)^2))
   
  out <- data.frame(nsim,k,t,p0,n,c,N.max=(k+1)*n,EN.avr,EN.avr.se,EN.LFC1=result.tbest$EN,EN.LFC1.se=result.tbest$se,
                       EN.LFC0=result.nobest$EN,EN.LFC0.se=result.nobest$se,PCS0=result.nobest$pcs,PCS1=result.tbest$pcs,delta1,delta2)
  out
}

############################
### k3t2, p0=0.1
############################
## detla1=0.05, P0*=P1*
out<-NULL
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=83,c=11,p0=0.1,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=101,c=13,p0=0.1,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=134,c=17,p0=0.1,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

## delta1=0, P0*=0.95
out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=64,c=8,p0=0.1,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=69,c=8,p0=0.1,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=83,c=9,p0=0.1,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

############################
### k3t2, p0=0.3
############################
## detla1=0.05, P0*=P1*
out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=135,c=18,p0=0.3,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=168,c=22,p0=0.3,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=223,c=29,p0=0.3,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

## delta1=0, P0*=0.95
out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=116,c=15,p0=0.3,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=131,c=16,p0=0.3,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=148,c=17,p0=0.3,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

############################
### k3t2, p0=0.4
############################
## detla1=0.05, P0*=P1*
out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=142,c=19,p0=0.4,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=181,c=24,p0=0.4,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=243,c=32,p0=0.4,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

############################
### k3t2, p0=0.5
############################
## detla1=0.05, P0*=P1*
out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=140,c=19,p0=0.5,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=179,c=24,p0=0.5,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=235,c=31,p0=0.5,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

## delta1=0, P0*=0.95
out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=127,c=17,p0=0.5,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=142,c=18,p0=0.5,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=161,c=19,p0=0.5,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

############################
### k3t2, p0=0.7
############################
## detla1=0.05, P0*=P1*
out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=98,c=14,p0=0.7,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=128,c=18,p0=0.7,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=167,c=23,p0=0.7,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

## delta1=0, P0*=0.95
out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=98,c=14,p0=0.7,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=110,c=15,p0=0.7,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=119,c=15,p0=0.7,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

############################
### k4t2, p0=0.1
############################
## detla1=0.05, P0*=P1*
out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=90,c=12,p0=0.1,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=113,c=15,p0=0.1,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=134,c=17,p0=0.1,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

## delta1=0, P0*=0.95
out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=64,c=8,p0=0.1,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=76,c=9,p0=0.1,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=83,c=9,p0=0.1,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

############################
### k4t2, p0=0.3
############################
## detla1=0.05, P0*=P1*
out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=147,c=20,p0=0.3,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=186,c=25,p0=0.3,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=242,c=32,p0=0.3,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

## delta1=0, P0*=0.95
out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=128,c=17,p0=0.3,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=137,c=17,p0=0.3,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=155,c=18,p0=0.3,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

############################
### k4t2, p0=0.4
############################
## detla1=0.05, P0*=P1*
out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=154,c=21,p0=0.4,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=193,c=26,p0=0.4,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=256,c=34,p0=0.4,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

############################
### k4t2, p0=0.5
############################
## detla1=0.05, P0*=P1*
out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=152,c=21,p0=0.5,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=191,c=26,p0=0.5,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=253,c=34,p0=0.5,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

## delta1=0, P0*=0.95
out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=139,c=19,p0=0.5,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=155,c=20,p0=0.5,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=171,c=21,p0=0.5,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

############################
### k4t2, p0=0.7
############################
## detla1=0.05, P0*=P1*
out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=104,c=15,p0=0.7,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=134,c=19,p0=0.7,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=179,c=25,p0=0.7,delta1=0.05,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

## delta1=0, P0*=0.95
out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=104,c=15,p0=0.7,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=116,c=16,p0=0.7,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

out<-read.csv(file="R_CW2_optimal_designs_sim.csv")[,-1]
tmp<-RCW2_sim(nsim=10000,k=4,t=2,n=131,c=17,p0=0.7,delta1=0,delta2=0.2) 
out<-rbind(out,tmp)
write.csv(out,file="R_CW2_optimal_designs_sim.csv")

###########################################################
### Example: Section 3.5
###########################################################
# optimal parameters using R_CW1
RCW1.fvalue(k=3,t=2,p0star=0.8,p1star=0.8,p0=0.7,delta1=-0.05,delta2=0.15,
            c_L=0,c_U=5,n_L=60,n_U=70,n_step=1) #n=67,c=5
# calculate EN using R_CW2
tmp<-RCW2_sim(nsim=10000,k=3,t=2,n=67,c=5,p0=0.7,delta1=-0.05,delta2=0.15) 




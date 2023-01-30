###### Exact calculation of P(CS) & E(N) for fixed-sample-size TSE two-stage allowing less than t Version A #######
###### Procedure:
######    Stage1. If X_{[k-t+1],1} >= X_{01}+c1, select the t best populations and use randomization to break a tie, then continue to S2.
######            Otherwise, select only populations with X_{i1}>=X_{01}+c1 and continue to S2.
######            If no population has X_{i1}>=X_{01}+c1, select none and terminate the experiment.
######    Stage2. Select only populations with X_{i1}+X_{i2}>=X_{01}+X_{02}+c2.

####################################### Version History ##########################
######  Original - 03OCT2022

## P(CS1) from the original two stage design combined for different k values
pcs1 <- function(k,n1,n2,c1,c2,p0,delta1,delta2){
  if (k<3 | k>5){stop('choose k value between 3 and 5')}
  s <- 0
  
  for (x01 in 0:n1){ # x01
    for (x11 in 0:n1){  # x11 the largest
      for (x21 in 0:n1){ # x21 the 2nd largest
        if (min(x11,x21)-x01>=c1)
        {
          a <- dbinom(x01,n1,p0)*dbinom(x11,n1,p0+delta2)*dbinom(x21,n1,p0+delta2)
          if (k==3){
            b <- pbinom(min(x11-1,x21-1),n1,p0+delta1)+ 1/2*dbinom(min(x11,x21),n1,p0+delta1)*(x11!=x21)+
                1/3* dbinom(x11,n1,p0+delta1)*(x11==x21)
          }
          else if (k==4){
            b<-(pbinom(min(x11-1,x21-1),n1,p0+delta1))^2+ 
              2*(1/2*(x11!=x21)+1/3*(x11==x21))*pbinom(min(x11-1,x21-1),n1,p0+delta1)*dbinom(min(x11,x21),n1,p0+delta1)+
              (1/3*(x11!=x21)+1/6*(x11==x21))*(dbinom(min(x11,x21),n1,p0+delta1))^2
          }
          else if (k==5){
            b<-(pbinom(min(x11-1,x21-1),n1,p0+delta1))^3+ 
              3*(1/2*(x11!=x21)+1/3*(x11==x21))*pbinom(min(x11-1,x21-1),n1,p0+delta1)^2*dbinom(min(x11,x21),n1,p0+delta1)+
              3*(1/3*(x11!=x21)+1/6*(x11==x21))*pbinom(min(x11-1,x21-1),n1,p0+delta1)*dbinom(min(x11,x21),n1,p0+delta1)^2+
              (1/4*(x11!=x21)+1/10*(x11==x21))*dbinom(min(x11,x21),n1,p0+delta1)^3
          }
          
          
          dsum <- 0
          for (x02 in 0:n2){   # x02
            dd <- (1-pbinom(x01+x02+c2-x11-1,n2,p0+delta2)) * (1-pbinom(x01+x02+c2-x21-1,n2,p0+delta2)) * dbinom(x02,n2,p0)
            dsum <- dsum + dd
          }
          s <- s + a*b*dsum
        }
      }
    }
  }
  s
}

#pcs1(k=3,n1=40,n2=20,c1=5,c2=8,p0=0.1,delta1=0,delta2=0.2)


## P(CS0 | H0) for t=2, k from 3 to 5
pcs0 <- function(k,n1,n2,c1,c2,p0){
  if (k<3 | k>5){stop('choose k value between 3 and 5')}
  
  s1<-0
  for (x01 in 0:n1){
    m1<-dbinom(x01,n1,p0)
    s2<-(pbinom(x01+c1-1,n1,p0))^k
    
    for (x11 in 0:n1){
      m2<-dbinom(x11,n1,p0)
      m3<-k*(x11>=x01+c1)*(pbinom(x01+c1-1,n1,p0))^(k-1)

      m4<-0
      for (x02 in 0:n2){
        m4 <- m4+dbinom(x02,n2,p0)*pbinom(x01+x02+c2-x11-1,n2,p0)
      }
      
      s3<-0
      for (x21 in 0:n1){
        m5<-k*(k-1)/2 *dbinom(x21,n1,p0)*(min(x11,x21)>=x01+c1)
        if (k==3){
          m6<-pbinom(min(x11,x21)-1,n1,p0) + 
              (1/2*(x11!=x21)+1/3*(x11==x21))*dbinom(min(x11,x21),n1,p0) # P(E_V)
        }
        else if (k==4){
          m6<-pbinom(min(x11,x21)-1,n1,p0)^2+ 
            ((x11!=x21)+2/3*(x11==x21))*pbinom(min(x11,x21)-1,n1,p0)*dbinom(min(x11,x21),n1,p0)+
            (1/3*(x11!=x21)+1/6*(x11==x21))*dbinom(min(x11,x21),n1,p0)^2  # P(E_V)
        }
        else if (k==5){
          m6<-pbinom(min(x11,x21)-1,n1,p0)^3+ 
            (3/2*(x11!=x21)+(x11==x21))*pbinom(min(x11,x21)-1,n1,p0)^2*dbinom(min(x11,x21),n1,p0)+
            ((x11!=x21)+1/2*(x11==x21))*pbinom(min(x11,x21)-1,n1,p0)*dbinom(min(x11,x21),n1,p0)^2+
            (1/4*(x11!=x21)+1/10*(x11==x21))*dbinom(min(x11,x21),n1,p0)^3
        }
        m7<-0
        for (x02 in 0:n2){
          m7 <- m7+dbinom(x02,n2,p0)*pbinom(x01+x02+c2-x11-1,n2,p0)*pbinom(x01+x02+c2-x21-1,n2,p0)
        }
        s3<-s3+m5*m6*m7
      }
      
      s2<-s2+m2*(m3*m4+s3)
    }
    s1<-s1+m1*s2
  }
  s1
}

# alternative function (more time-consuming)
# pcs0_alt <- function(k,n1,n2,c1,c2,p0){
#   if (k!=3){stop('k value must be 3')}
#   s1<-0 #prob of selecting no pop at stage 1
#   for (x01 in 0:n1){
#     m1<-dbinom(x01,n1,p0)*(pbinom(x01+c1-1,n1,p0))^k
#     s1<-s1+m1
#   }
#   
#   s2<-0 #prob of selecting 1 pop at stage 1, then fail at stage 2
#   for (x01 in 0:n1){  
#     for (x11 in 0:n1){
#       m2<-k*dbinom(x01,n1,p0)*dbinom(x11,n1,p0)*(x11>=x01+c1)*(pbinom(x01+c1-1,n1,p0))^(k-1)
#       
#       m3<-0
#       for (x02 in 0:n2){
#         m3 <- m3+dbinom(x02,n2,p0)*pbinom(x01+x02+c2-x11-1,n2,p0)
#       }
#       s2<-s2+m2*m3
#     }
#   }
#       
#   s3<-0 #prob of selecting 2 pop at stage 1, then fail at stage 2
#   for (x01 in 0:n1){  
#     for (x11 in 0:n1){
#       for (x21 in 0:n1){
#         m4<-k*(k-1)/2 *dbinom(x01,n1,p0)*dbinom(x11,n1,p0)*dbinom(x21,n1,p0)*(min(x11,x21)>=x01+c1)
#         m5<-pbinom(min(x11,x21)-1,n1,p0) + 1/2*dbinom(min(x11,x21),n1,p0)*(x11!=x21) 
#         + 1/3*dbinom(x11,n1,p0)*(x11==x21)  # P(E_V)
#         m6<-0
#         for (x02 in 0:n2){
#           m6 <- m6+dbinom(x02,n2,p0)*pbinom(x01+x02+c2-x11-1,n2,p0)*pbinom(x01+x02+c2-x21-1,n2,p0)
#         }
#         s3<-s3+m4*m5*m6
#       }
#     }
#   }
#   s1+s2+s3
# }

#pcs0(k=3,n1=40,n2=20,c1=5,c2=8,p0=0.1)

#########################################################
####### Probability of early termination
#########################################################
tau <-function(k,n1,n2,c1,p0,delta1,delta2){
  s1<-0 #prob of selecting no pop at stage 1
  for (x01 in 0:n1){
    m1<-dbinom(x01,n1,p0)*pbinom(x01+c1-1,n1,p0+delta1)^(k-2)*pbinom(x01+c1-1,n1,p0+delta2)^2
    s1<-s1+m1
  }
  s1
}
  
  
#########################################################################
####### E(S): # of exp. populations selected at stage 1  & E(N) when t=2
#########################################################################
ES.s1 <-function(k,n1,n2,c1,p0,delta1,delta2){
  # t=2 only, g=2 or 0 (when setting delta2=0)
  
  s2<-0 #prob of selecting 1 pop at stage 1
  for (x01 in 0:n1){
    for (x11 in 0:n1){
      m2<-(k-2)*dbinom(x01,n1,p0)*dbinom(x11,n1,p0+delta1)*(x11>=x01+c1)*pbinom(x01+c1-1,n1,p0+delta1)^(k-3)*pbinom(x01+c1-1,n1,p0+delta2)^2+
        2*dbinom(x01,n1,p0)*dbinom(x11,n1,p0+delta2)*(x11>=x01+c1)*pbinom(x01+c1-1,n1,p0+delta1)^(k-2)*pbinom(x01+c1-1,n1,p0+delta2)
      
      s2<-s2+m2
    }
  }
  
  s3<-0 #prob of selecting 2 pop at stage 1
  for (x01 in 0:n1){
    for (x11 in 0:n1){
      for (x21 in 0:n1){
        if (k==3){
          m4<-
            # one good, one bad
            2*dbinom(x01,n1,p0)*dbinom(x11,n1,p0+delta1)*dbinom(x21,n1,p0+delta2)*(min(x11,x21)>=x01+c1)*
            (pbinom(min(x11,x21)-1,n1,p0+delta2) + 1/2*dbinom(min(x11,x21),n1,p0+delta2)*(x11!=x21)+ 1/3*dbinom(x11,n1,p0+delta2)*(x11==x21))+
            # both are good
            dbinom(x01,n1,p0)*dbinom(x11,n1,p0+delta2)*dbinom(x21,n1,p0+delta2)*(min(x11,x21)>=x01+c1)*
            (pbinom(min(x11,x21)-1,n1,p0+delta1) + 1/2*dbinom(min(x11,x21),n1,p0+delta1)*(x11!=x21)+ 1/3*dbinom(x11,n1,p0+delta1)*(x11==x21))
        }
        else if (k==4){
          m4<-
            # both are bad
            dbinom(x01,n1,p0)*dbinom(x11,n1,p0+delta1)*dbinom(x21,n1,p0+delta1)*(min(x11,x21)>=x01+c1)*
                (pbinom(min(x11,x21)-1,n1,p0+delta2)^2+ 
                ((x11!=x21)+2/3*(x11==x21))*pbinom(min(x11,x21)-1,n1,p0+delta2)*dbinom(min(x11,x21),n1,p0+delta2)+
                (1/3*(x11!=x21)+1/6*(x11==x21))*dbinom(min(x11,x21),n1,p0+delta2)^2) +
            # one good, one bad
            4*dbinom(x01,n1,p0)*dbinom(x11,n1,p0+delta1)*dbinom(x21,n1,p0+delta2)*(min(x11,x21)>=x01+c1)*
                (pbinom(min(x11,x21)-1,n1,p0+delta1)*pbinom(min(x11,x21)-1,n1,p0+delta2)+ 
                (1/2*(x11!=x21)+1/3*(x11==x21))*pbinom(min(x11,x21)-1,n1,p0+delta1)*dbinom(min(x11,x21),n1,p0+delta2)+
                (1/2*(x11!=x21)+1/3*(x11==x21))*pbinom(min(x11,x21)-1,n1,p0+delta2)*dbinom(min(x11,x21),n1,p0+delta1)+
                (1/3*(x11!=x21)+1/6*(x11==x21))*dbinom(min(x11,x21),n1,p0+delta1)*dbinom(min(x11,x21),n1,p0+delta2)) +
            # both are good
            dbinom(x01,n1,p0)*dbinom(x11,n1,p0+delta2)*dbinom(x21,n1,p0+delta2)*(min(x11,x21)>=x01+c1)*
                (pbinom(min(x11,x21)-1,n1,p0+delta1)^2+ 
                ((x11!=x21)+2/3*(x11==x21))*pbinom(min(x11,x21)-1,n1,p0+delta1)*dbinom(min(x11,x21),n1,p0+delta1)+
                (1/3*(x11!=x21)+1/6*(x11==x21))*dbinom(min(x11,x21),n1,p0+delta1)^2)
        }
        else if (k==5){
          m4<-
            # both are bad
            3*dbinom(x01,n1,p0)*dbinom(x11,n1,p0+delta1)*dbinom(x21,n1,p0+delta1)*(min(x11,x21)>=x01+c1)* 
                (pbinom(min(x11,x21)-1,n1,p0+delta1)*pbinom(min(x11,x21)-1,n1,p0+delta2)^2+ 
                (1/2*(x11!=x21)+1/3*(x11==x21))*pbinom(min(x11,x21)-1,n1,p0+delta2)^2*dbinom(min(x11,x21),n1,p0+delta1)+
                2*(1/2*(x11!=x21)+1/3*(x11==x21))*pbinom(min(x11,x21)-1,n1,p0+delta1)*pbinom(min(x11,x21)-1,n1,p0+delta2)*dbinom(min(x11,x21),n1,p0+delta2)+
                (1/3*(x11!=x21)+1/6*(x11==x21))*pbinom(min(x11,x21)-1,n1,p0+delta1)*dbinom(min(x11,x21),n1,p0+delta2)^2+
                2*(1/3*(x11!=x21)+1/6*(x11==x21))*pbinom(min(x11,x21)-1,n1,p0+delta2)*dbinom(min(x11,x21),n1,p0+delta1)*dbinom(min(x11,x21),n1,p0+delta2)+
                (1/4*(x11!=x21)+1/10*(x11==x21))*dbinom(min(x11,x21),n1,p0+delta1)*dbinom(min(x11,x21),n1,p0+delta2)^2) +
            # one good, one bad
            6*dbinom(x01,n1,p0)*dbinom(x11,n1,p0+delta1)*dbinom(x21,n1,p0+delta2)*(min(x11,x21)>=x01+c1)*
                (pbinom(min(x11,x21)-1,n1,p0+delta1)^2*pbinom(min(x11,x21)-1,n1,p0+delta2)+ 
                (1/2*(x11!=x21)+1/3*(x11==x21))*pbinom(min(x11,x21)-1,n1,p0+delta1)^2*dbinom(min(x11,x21),n1,p0+delta2)+
                2*(1/2*(x11!=x21)+1/3*(x11==x21))*pbinom(min(x11,x21)-1,n1,p0+delta1)*pbinom(min(x11,x21)-1,n1,p0+delta2)*dbinom(min(x11,x21),n1,p0+delta1)+
                2*(1/3*(x11!=x21)+1/6*(x11==x21))*pbinom(min(x11,x21)-1,n1,p0+delta1)*dbinom(min(x11,x21),n1,p0+delta1)*dbinom(min(x11,x21),n1,p0+delta2) +
                (1/3*(x11!=x21)+1/6*(x11==x21))*pbinom(min(x11,x21)-1,n1,p0+delta2)*dbinom(min(x11,x21),n1,p0+delta1)*dbinom(min(x11,x21),n1,p0+delta1) +
                (1/4*(x11!=x21)+1/10*(x11==x21))*dbinom(min(x11,x21),n1,p0+delta2)*dbinom(min(x11,x21),n1,p0+delta1)^2) +
            # both are good
            dbinom(x01,n1,p0)*dbinom(x11,n1,p0+delta2)*dbinom(x21,n1,p0+delta2)*(min(x11,x21)>=x01+c1)*
            (pbinom(min(x11,x21)-1,n1,p0+delta1)^3+ 
               3*(1/2*(x11!=x21)+1/3*(x11==x21))*pbinom(min(x11,x21)-1,n1,p0+delta1)^2*dbinom(min(x11,x21),n1,p0+delta1)+
               3*(1/3*(x11!=x21)+1/6*(x11==x21))*pbinom(min(x11,x21)-1,n1,p0+delta1)*dbinom(min(x11,x21),n1,p0+delta1)^2+
               (1/4*(x11!=x21)+1/10*(x11==x21))*dbinom(min(x11,x21),n1,p0+delta1)^3)
        }
        s3<-s3+m4
      }
    }
  }
  ES.s1<-s2+2*s3 # E(S)= P(S=1)+2*P(S=2)+...+t*P(S=t)
  ES.s1
}

# ES.s1(k=4,n1=34,n2=44,c1=2,p0=0.1,delta1=0,delta2=0.2)
# ES.s1(k=4,n1=34,n2=44,c1=2,p0=0.1,delta1=0,delta2=0)


EN <- function(k,n1,n2,c1,p0,delta1,delta2){
  s2=0
  # EN=(k+1)*n1 + n2*E(S) + n2*P(max_{x_i1}>=X_{01}+c_1), S:# of exp. pop. seleted at stage 1 
  ES.s1.H0 <- ES.s1(k,n1,n2,c1,p0,delta1=0,delta2=0)
  tau.H0<-tau(k,n1,n2,c1,p0,delta1=0,delta2=0)
  ES.s1.LFC<- ES.s1(k,n1,n2,c1,p0,delta1,delta2)
  tau.LFC<-tau(k,n1,n2,c1,p0,delta1,delta2)
  
  EN.H0<-(k+1)*n1+n2*ES.s1.H0+n2*(1-tau.H0)
  EN.LFC <- (k+1)*n1+n2*ES.s1.LFC+n2*(1-tau.LFC)
  EN.avr <- 1/2*(EN.H0+EN.LFC)  # EN = 1/2[(k+1)*n1+(t+1)*n2*P(T1>=c1|H0)]+ 1/2*[(k+1)*n1+(t+1)*n2*P(T1>=c1|p*)]
  Nmax <- (k+1)*n1+3*n2
  data.frame(EN.avr,EN.H0,EN.LFC,Nmax,tau.H0)
}

#EN(k=3,n1=40,n2=20,c1=5,p0=0.1,delta1=0,delta2=0.2)
#EN(k=5,n1=34,n2=44,c1=2,p0=0.1,delta1=0,delta2=0.2)

############################################################
###### Combine outputs 
############################################################
twostage.out <- function(k,n1,n2,c1,c2,p0,delta1,delta2){
  power <- pcs1(k,n1,n2,c1,c2,p0,delta1,delta2)
  typ1err<-1-pcs0(k,n1,n2,c1,c2,p0)
  EN_vec<-EN(k,n1,n2,c1,p0,delta1,delta2)
  
  data.frame(k,t=2,p0,n1,n2,c1,c2,delta1,delta2,power,typ1err,
             EN.avr=EN_vec$EN.avr,EN.H0=EN_vec$EN.H0,EN.LFC=EN_vec$EN.LFC,Nmax=EN_vec$Nmax,tau=EN_vec$tau.H0)
}

#twostage.out(k=3,n1=40,n2=20,c1=5,c2=8,p0=0.1,delta1=0,delta2=0.2)
#twostage.out(k=4,n1=34,n2=44,c1=2,c2=9,p0=0.1,delta1=0,delta2=0.2)

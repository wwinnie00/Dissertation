##############################################################################################
######## Subset selection one-stage design (Version 6), exact P(CS1), P(V), P(V') and E(S) 
######## Goal: select a subset containing the t best given that at least t are acceptable; 
########        reject marginal populations if g<t are acceptable.
######## Procedure: (d>=0)
######## if X_(k-t+1)-X0>=c, select all populations for which  Xi>=X_(k-t+1)-d
######## else, reject all populations for which X_i<X0+c
##############################################################################################
## Original version: 31OCT2022

###################################################
####     P(CS1) at p*: p=...=pk=p0+delta2      ####
###################################################
PCS1.exact <- function(k,t,n,c,d,p0,delta){
  # p1=...pk=p0+delta
  # d>=0
  s1<-0
  for (z in 0:n){ # X_{[k-t+1]}=z
    for (m in 0:(t-1)){ # number of pop. with X_i>z regardless of good or bad
      for (l in 0:(k-t)){ # number of pop with X_i<z regardless of good or bad
        for (u in 0:m){ # number of pop. among the t best with X_i>z
          for (v in 0:min(t-u,l)){ # number of pop. among the t best with X_i<z
            a1<-choose(t,u)*choose(k-t,m-u)*choose(t-u,v)*choose(k-t-m+u,l-v)*
              (1-pbinom(z,n,p0+delta))^m *pbinom(z-1,n,p0+delta)^(l-v) *dbinom(z,n,p0+delta)^(k-m-l)*
              (pbinom(z-1,n,p0+delta)-pbinom(z-d-1,n,p0+delta))^v *pbinom(z-c,n,p0)
            s1<-s1+a1
          }
        }
      } 
    }
  }
  s1
}
# PCS1.exact(k=3,t=2,n=30,c=5,d=1,p0=0.3,delta=0.2)

###################################################
#### Conservative P(V) at p1=...=pk=p0+delta2 #####
###################################################

# refer to Thm 5.2
PV.exact <- function(k,t,n,c,d,p0,delta){
  # p1=...pk=p0+delta
  # d>=0
  # any k, t=1 or 2
  
  sum<-0
  if (t==1){
    for (xk in 0:n){
      tmp<-(pbinom(xk+d,n,p0+delta))^(k-t) *pbinom(xk-c,n,p0) *
        dbinom(xk,n,p0+delta)
      sum<-sum+tmp
    }
  }
  
  if (t==2){
    for (xj in 0:n){ #x_{k-1}
      for (xk in 0:n){
        xmin <- min(xj,xk)
        tmp<-(pbinom(xmin+d,n,p0+delta))^(k-t) *pbinom(xmin-c,n,p0) *
          dbinom(xj,n,p0+delta) *dbinom(xk,n,p0+delta)
        sum<-sum+tmp
      }
    }
  }
  
  if (t==3){
    for (xi in 0:n){ #x_{k-2}
      for (xj in 0:n){ #x_{k-1}
        for (xk in 0:n){
          xmin <- min(xi,xj,xk)
          tmp<-(pbinom(xmin+d,n,p0+delta))^(k-t) *pbinom(xmin-c,n,p0) *
            dbinom(xi,n,p0+delta) *dbinom(xj,n,p0+delta) *dbinom(xk,n,p0+delta)
          sum<-sum+tmp
        }
      }
    }
  }
  sum
}
#PV.exact(k=3,t=2,n=30,c=5,d=1,p0=0.3,delta=0.2)


###################################################
#### Modified LB P(V') at p=...=pk=p0+delta2  #####
###################################################
PV.mod.exact <- function(k,t,n,c,d,p0,delta){
  # p1=...pk=p0+delta
  # d>=0

  sum<-0
  for (x0 in 0:n){
    for (z in 0:n){ # X_{[k-t+1]}=z
      for (m in 0:(t-1)){ # number of pop. with X_i>z regardless of good or bad
        for (l in 0:(k-t)){ # number of pop with X_i<z regardless of good or bad
          for (u in 0:m){ # number of pop. among the t best with X_i>z
            for (v in 0:min(t-u,l)){ # number of pop. among the t best with X_i<z
              tmp<-choose(t,u)*choose(k-t,m-u)*choose(t-u,v)*choose(k-t-m+u,l-v)*
                (1-pbinom(z,n,p0+delta))^(m-u) *(1-pbinom(max(z,x0+c-1),n,p0+delta))^u *
                pbinom(z-1,n,p0+delta)^(l-v) *dbinom(z,n,p0+delta)^(k-m-l)*
                (pbinom(z-1,n,p0+delta)-pbinom(max(x0+c-1,z-d-1),n,p0+delta))^v *(x0+c<=z) *dbinom(x0,n,p0)
              sum<-sum+tmp
              #print(c(a1,s1))
            }
          }
        } 
      }
    }
  }
  sum
}

# PV.mod.exact(k=2,t=1,n=30,c=5,d=1,p0=0.3,delta=0.2)

## altervative calculation of P(V') using the derived formula in thesis
# PV.mod.exact <- function(k,t,n,c,d,p0,delta){
#   # p1=...pk=p0+delta
#   # d>=0
#   # k=3,4
#   # t=2
#   sum<-0
#   if (k<3|k>4|t!=2){stop('Choose valid values of k and t!')}
#   else{
#     if (k==3){
#       for (x1 in (0:n)){
#         for (x2 in (0:n)){
#           for (x3 in (0:n)){
#             xt <- sort(c(x1,x2,x3))[k-t+1] #tth largest Xi
#             minx.A <- min(x2,x3) # minX_i for i in A={k-t+1,...,k}
#             tmp<-dbinom(x1,n,p0+delta)*dbinom(x2,n,p0+delta)*dbinom(x3,n,p0+delta)*pbinom(minx.A-c,n,p0)*(minx.A>=xt-d)
#             sum<-sum+tmp
#           }
#         }
#       }
#     }
# 
#     if (k==4){
#       for (x1 in (0:n)){
#         for (x2 in (0:n)){
#           for (x3 in (0:n)){
#             for (x4 in (0:n)){
#               xt <- sort(c(x1,x2,x3,x4))[k-t+1] #tth largest Xi
#               minx.A <- min(x3,x4) # minX_i for i in A={k-t+1,...,k}
#               tmp<-dbinom(x1,n,p0+delta)*dbinom(x2,n,p0+delta)*dbinom(x3,n,p0+delta)*dbinom(x4,n,p0+delta)*pbinom(minx.A-c,n,p0)*(minx.A>=xt-d)
#               sum<-sum+tmp
#             }
#           }
#         }
#       }
#     }
#   }
#   sum
# }
# PV.mod.exact(k=4,t=2,n=30,c=5,d=1,p0=0.3,delta=0.2)

###################################################
####     E(s) conditional on s>t: E(S|S>=t)    ####
###################################################
ES.1.exact <- function(k,t,n,c,d,p0,delta1=delta2,delta2){
  # d>=0
  # k>=1, t>=1
  num<-denom<-0
  
  if (delta1==delta2){# at p*: p1=...pk=p0+delta2

      
    if (t==1){ 
      for (z in 0:n){ # X_{[k]}=z
        s1<-s2<-s4<-0
        
        for (l in 0:(k-2)){ # number of pop with X_i<z except X1 
          a1<-choose(k-1,l) *pbinom(z-1,n,p0+delta2)^l *(pbinom(z-1,n,p0+delta2)-pbinom(z-d-1,n,p0+delta2)) *
            dbinom(z,n,p0+delta2)^(k-l-1)
          s1<-s1+a1
        } 
        for (l in 0:(k-1)){ # number of pop with X_i<z  
          a2<-choose(k-1,l) *pbinom(z-1,n,p0+delta2)^l *dbinom(z,n,p0+delta2)^(k-l)
          s2<-s2+a2
        }

        for (l in 0:(k-1)){ # number of pop with X_i<z  
          a4<-choose(k,l) *pbinom(z-1,n,p0+delta2)^l *dbinom(z,n,p0+delta2)^(k-l)
          s4<-s4+a4
          }
        
        num<-num+k*(s1+s2)*pbinom(z-c,n,p0)
        denom<-denom+s4*pbinom(z-c,n,p0)
      }
    }
    
    else{ #t>=2
      for (z in 0:n){ # X_{[k-t+1]}=z
        s1<-s2<-s3<-s4<-0
        for (m in 0:(t-2)){ # number of pop. with X_i>z except X1
          for (l in 0:(k-t)){ # number of pop with X_i<z  
          a1<-choose(k-1,m)*choose(k-m-1,l)*(1-pbinom(z,n,p0+delta2))^(m+1)*
            pbinom(z-1,n,p0+delta2)^l * dbinom(z,n,p0+delta2)^(k-m-l-1)
          s1<-s1+a1
          } 
        }
        for (m in 0:(t-1)){ # number of pop. with X_i>z 
          for (l in 0:(k-t-1)){ # number of pop with X_i<z except X1 
            a2<-choose(k-1,m)*choose(k-m-1,l)*(1-pbinom(z,n,p0+delta2))^m*
              pbinom(z-1,n,p0+delta2)^l *(pbinom(z-1,n,p0+delta2)-pbinom(z-d-1,n,p0+delta2))*
              dbinom(z,n,p0+delta2)^(k-m-l-1)
            s2<-s2+a2
          }
        }
        for (m in 0:(t-1)){ # number of pop. with X_i>z 
          for (l in 0:(k-t)){ # number of pop with X_i<z  
            a3<-choose(k-1,m)*choose(k-m-1,l)*(1-pbinom(z,n,p0+delta2))^m*
              pbinom(z-1,n,p0+delta2)^l * dbinom(z,n,p0+delta2)^(k-m-l)
            s3<-s3+a3
            
            a4<-choose(k,m)*choose(k-m,l)*(1-pbinom(z,n,p0+delta2))^m*
              pbinom(z-1,n,p0+delta2)^l * dbinom(z,n,p0+delta2)^(k-m-l)
            s4<-s4+a4
          }
        }
        num<-num+k*(s1+s2+s3)*pbinom(z-c,n,p0)
        denom<-denom+s4*pbinom(z-c,n,p0)
      }
    }
  }
  
  # alternative calculation for E(s) cost more time
  else{ # p1=...=p_{k-t}=p0+delta1, p_{k-t+1}=...=p_k=delta2
      if (k<2|k>5){stop('Choose a valid value of k between 2 and 5!')}
      else{
        p.vec<-c(rep(p0+delta1,k-t),rep(p0+delta2,t))
        
        if (k==2){
          for (x1 in (0:n)){
            for (x2 in (0:n)){
              xt <- max(x1,x2) #tth largest Xi
              aa<-dbinom(x1,n,p.vec[1])*dbinom(x2,n,p.vec[2])*pbinom(xt-c,n,p0)
              num<-num+aa*((x1>=xt-d)+(x2>=xt-d))
              denom<-denom+aa
            }
          }
        }
        
        if (k==3){
          for (x1 in (0:n)){
            for (x2 in (0:n)){
              for (x3 in (0:n)){
                xt <- sort(c(x1,x2,x3))[k-t+1] #tth largest Xi
                aa<-dbinom(x1,n,p.vec[1])*dbinom(x2,n,p.vec[2])*dbinom(x3,n,p.vec[3])*pbinom(xt-c,n,p0)
                num<-num+aa*((x1>=xt-d)+(x2>=xt-d)+(x3>=xt-d))
                denom<-denom+aa
              }
            }
          }
        }

        if (k==4){
          for (x1 in (0:n)){
            for (x2 in (0:n)){
              for (x3 in (0:n)){
                for (x4 in (0:n)){
                  xt <- sort(c(x1,x2,x3,x4))[k-t+1] #tth largest Xi
                  aa<-dbinom(x1,n,p.vec[1])*dbinom(x2,n,p.vec[2])*dbinom(x3,n,p.vec[3])*dbinom(x4,n,p.vec[4])*pbinom(xt-c,n,p0)
                  num<-num+aa*((x1>=xt-d)+(x2>=xt-d)+(x3>=xt-d)+(x4>=xt-d))
                  denom<-denom+aa
                }
              }
            }
          }
        }

        if (k==5){
          for (x1 in (0:n)){
            for (x2 in (0:n)){
              for (x3 in (0:n)){
                for (x4 in (0:n)){
                  for (x5 in (0:n)){
                  xt <- sort(c(x1,x2,x3,x4,x5))[k-t+1] #tth largest Xi
                  aa<-dbinom(x1,n,p.vec[1])*dbinom(x2,n,p.vec[2])*dbinom(x3,n,p.vec[3])*dbinom(x4,n,p.vec[4])*dbinom(x5,n,p.vec[5])*pbinom(xt-c,n,p0)
                  num<-num+aa*((x1>=xt-d)+(x2>=xt-d)+(x3>=xt-d)+(x4>=xt-d)+(x5>=xt-d))
                  denom<-denom+aa
                  }
                }
              }
            }
          }
        }
      }
  }
      
  data.frame(num,denom,ES1=num/denom)
}

#ES.1.exact(k=3,t=2,n=30,c=5,d=1,p0=0.3,delta2=0.2)

###################################################
####     E(s) conditional on s<t: E(s|s<t)     ####
###################################################

ES.0.exact <- function(k,t,n,c,d,p0,delta){
  # p1=...pk=p0+delta
  # d>=0
  # k>=3, t>=2
  
  if (t<2){stop('Choose a valid value of t>=2 !')}
  
  num<-denom<-0

  for (z in 0:n){ # X_{[k-t+1]}=z
    s1<-s2<-0
    
    for (x0 in (0:n)){
      for (m in 0:(t-2)){ # number of pop. with X_i>z except X1
        for (l in 0:(k-t)){ # number of pop with X_i<z  
          a1<-choose(k-1,m)*choose(k-m-1,l)*(1-pbinom(z,n,p0+delta))^m*
            pbinom(z-1,n,p0+delta)^l * dbinom(z,n,p0+delta)^(k-m-l-1)*
            dbinom(x0,n,p0) * (1-pbinom(x0+c-1,n,p0+delta)) * (z<x0+c)
          s1<-s1+a1
        } 
      }
    }

    for (m in 0:(t-1)){ # number of pop. with X_i>z 
      for (l in 0:(k-t)){ # number of pop with X_i<z  
        a2<-choose(k,m)*choose(k-m,l)*(1-pbinom(z,n,p0+delta))^m*
          pbinom(z-1,n,p0+delta)^l * dbinom(z,n,p0+delta)^(k-m-l)*
          (1-pbinom(z-c,n,p0))
        s2<-s2+a2
      }
    }
    
    num<-num+k*s1
    denom<-denom+s2
  }
  
  data.frame(num,denom,ES0=num/denom)
}
# ES.0.exact(k=3,t=2,n=30,c=5,d=0,p0=0.3,delta=0.2)


###################################################
####            Unconditional E(S)             ####
###################################################
ES.exact <- function(k,t,n,c,d,p0,delta){
  # p1=...pk=p0+delta
  # d>=0
  # k>=3, t>=2
  ES1.num<-ES.1.exact(k,t,n,c,d,p0,delta2=delta)$num
  ES0.num<-ifelse(t==1,0,ES.0.exact(k,t,n,c,d,p0,delta)$num)
  ES <- ES1.num+ES0.num
  ES
}

# ES.exact(k=3,t=2,n=30,c=5,d=0,p0=0.3,delta=0.2)


# ES.0.exact <- function(k,t,n,c,d,p0,delta){
#   # d>=0
#   # k=3,4,5
#   # t<=k
#   num<-denom<-0
#   if (k!=3){stop('Choose a valid value of k!')}
#   else{
#     if (k==3){
#       for (x1 in (0:n)){
#         for (x2 in (0:n)){
#           for (x3 in (0:n)){
#             xt <- sort(c(x1,x2,x3))[k-t+1] #tth largest Xi
#             aa<-dbinom(x1,n,p0+delta)*dbinom(x2,n,p0+delta)*dbinom(x3,n,p0+delta)*
#               ( (x1>xt)*(pbinom(x1-c,n,p0)-pbinom(xt-c,n,p0))+
#                   (x2>xt)*(pbinom(x2-c,n,p0)-pbinom(xt-c,n,p0))+
#                   (x3>xt)*(pbinom(x3-c,n,p0)-pbinom(xt-c,n,p0)) )
#             num<-num+aa
#             
#             bb<-dbinom(x1,n,p0+delta)*dbinom(x2,n,p0+delta)*dbinom(x3,n,p0+delta)*(1-pbinom(xt-c,n,p0))
#             denom<-denom+bb
#           }
#         }
#       }
#     }
#   }
#   data.frame(num,denom,ES0=num/denom)
# }
# 
# ES.0.exact(k=3,t=2,n=30,c=5,d=0,p0=0.3,delta=0.2)

###################################################################
######## R_CW1: modified for selecting less than t
######## g: number of populations p_i>=p_0+delta_2
###################################################################
################## Original version: 2022/07/29 ###################

# P(CS_{0,g}) when g<t, any k,t
RCW1.pcs.seq <- function(k,t,n,c,p0,delta1,delta2){
  
  pcs.seq<-numeric()
  
  for (g in 0:(t-1)){
    prob <- 0
    for (x0 in 0:n){# x0 value
      #pp <- (pbinom(x0+c-1,n,p0+delta1))^(k-g)*(1-pbinom(x0+c-1,n,p0+delta2))^g*dbinom(x0,n,p0) #select all good ones and exclude all bad ones
      pp <- (pbinom(x0+c-1,n,p0+delta1))^(k-g)*dbinom(x0,n,p0) #exclude bad ones
      prob <- prob + pp
    }
    pcs.seq[g+1]<-prob
  }
  pcs.seq
}

#RCW1.pcs.seq(k=3,t=2,n=88,c=11,p0=0.1,delta1=0.05,delta2=0.2)
#RCW1.pcs.seq(k=2,t=1,n=56,c=9,p0=0.2,delta1=0,delta2=0.2)

# Exact P(CS1) when t=1
dun.pcs1.t1 <- function(k,n,c,p0,delta1,delta2){
  sum <-0
  for (xk in 0:n){ 
    
    inner.sum <- 0
    for (m in 0:(k-1)){ #m is number of bad populations having exactly x_k successes (tied with pi_k)
      inner.sum <- inner.sum + 1/(m+1)*(dbinom(xk,n,p0+delta1))^m*(pbinom(xk-1,n,p0+delta1))^(k-1-m)
    }
    sum <- sum + dbinom(xk,n,p0+delta2)*pbinom(xk-c,n,p0)*inner.sum
  }
  sum
}
#dun.pcs1.t1(k=3,n=50, c=6, p0=0.2, delta1=0.05, delta2=0.2)

# Exact P(CS) when t=2
dun.pcs1.t2 <- function(k,n,c,p0,delta1,delta2){
  
  sum <- 0
  for (y in 0:n){ #x_{k-1}
    for (xk in 0:n){ 
      z <- min(y,xk)
      s<-ifelse(y==xk,2,1) #s is the number of good populations having exactly z successes
      
      inner.sum <- 0
      for (m in 0:(k-2)){ #m is number of bad populations having exactly z successes (tied with good populations)
        inner.sum <- inner.sum + 1/choose(m+s,s)*(dbinom(z,n,p0+delta1))^m*(pbinom(z-1,n,p0+delta1))^(k-2-m)
      }
      sum <- sum + dbinom(y,n,p0+delta2)*dbinom(xk,n,p0+delta2)*inner.sum*pbinom(z-c,n,p0)
    }
  }
  sum
}

# Exact P(CS) when t=3
dun.pcs1.t3 <- function(k,n,c,p0,delta1,delta2){
  
  sum <- 0
  for (y1 in 0:n){ #x_{k-2}
    for (y2 in 0:n){ #x_{k-1}
      for (xk in 0:n){ 
        y.vec<-c(y1,y2,xk)
        z <- min(y.vec)
        s<-length(which(y.vec==z)) #s is the number of good populations having exactly z successes
        
        inner.sum <- 0
        for (m in 0:(k-3)){ #m is number of bad populations having exactly z successes (tied with good populations)
          inner.sum <- inner.sum + 1/choose(m+s,s)*(dbinom(z,n,p0+delta1))^m*(pbinom(z-1,n,p0+delta1))^(k-3-m)
        }
        sum <- sum + dbinom(y1,n,p0+delta2)*dbinom(y2,n,p0+delta2)*dbinom(xk,n,p0+delta2)*inner.sum*pbinom(z-c,n,p0)
      }
    }
  }
  sum
}

#dun.pcs1.t3(k=4,n=50, c=6, p0=0.2, delta1=0.05, delta2=0.2)

### Print results based on exact probability calculations
onestage.design.out <- function(k,t,n,c, p0, delta1=0.05, delta2=0.2){
  # delta1: minimal improvement
  # delta2: acceptable improvement
  pcs.seq <- RCW1.pcs.seq(k,t,n,c,p0,delta1,delta2)
  
  if (t>3){stop("Undefined value of t!")}
  else{
    if (t==1){
      pcs1<-dun.pcs1.t1(k,n, c, p0, delta1, delta2)
    }
    if (t==2){
      pcs1<-dun.pcs1.t2(k,n, c, p0, delta1, delta2)
    }
    if (t==3){
      pcs1<-dun.pcs1.t3(k,n, c, p0, delta1, delta2)
    }
  }
  final.out <- data.frame(k,t,p0,delta1,delta2,n,c,N=(k+1)*n,pcs00=pcs.seq[1],pcs02=pcs.seq[2],pcs1=pcs1)
  final.out
}

#onestage.design.out(k=3,t=2,n=83,c=11,p0=0.1,delta1=0.05,delta2=0.2)
#onestage.design.out(k=4,t=2,n=90,c=12,p0=0.1,delta1=0.05,delta2=0.2)

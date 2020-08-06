# Created by Allan Strand
#
# broken stick regression model
#
#
bs <- function(X,Y) 
  {

#    initial values
    p <- rep(0,4)
    p[1:2] <- c(quantile(X,0.25),mean(Y)) #point 1
    p[3:4] <- c(quantile(X,0.75),mean(Y)) #point 2
    optim(p,llike.bs,X=X,Y=Y)
  }

bs.new <- function(X,Y)
  {

#    initial values
    p <- rep(0,5)
    p[1:2] <- c(quantile(X,0.25),mean(Y)-0.1) #point 1
    p[3:4] <- c(quantile(X,0.75),mean(Y)+0.1) #point 2
    p[5] <- 1 #sd
    optim(p,llike.bs,X=X,Y=Y)
  }

flt.old <- function(X,Y)
  {
    optim(p=c(0.5,1),llike.flt,X=X,Y=Y)
  }

flt <- function(X,Y)
  {
    optim(p=c(0.5),llike.flt,X=X,Y=Y)
  }


llike.bs <- function(p,X,Y)# current objective function but sets sd to
                              #the sd of allele freqs across cline
  {
    p1 <- p[1:2]
    p2 <- p[3:4]
    sdY <- sd(Y)
#    sd=0.2
    if ((p1[1]<p2[1])&(sdY>0)) {ok <- T} else {ok <- F}# points not swapped, sd ok
    if (ok & ((p1[1]>min(X))&(p2[1]<max(X)))) {ok <- T} else {ok <- F}# points not at end of line
    if (ok & (min(c(p1[2],p2[2])>=0))) {ok <- T} else {ok <- F}
    if (ok & (max(c(p1[2],p2[2])<=1))) {ok <- T} else {ok <- F}
    
    if (ok)
      {
        leftX <- X[which(X<p1[1])]
        leftY <- Y[which(X<p1[1])]
        rightX <- X[which(X>p2[1])]
        rightY <- Y[which(X>p2[1])]
        middleX <- X[which((X>=p1[1])&(X<=p2[1]))]
        middleY <- Y[which((X>=p1[1])&(X<=p2[1]))]
        
        resleft <- p1[2]-leftY
        resright <- p2[2]-rightY
        
        slp <- (p2-p1)[2]/(p2-p1)[1] #rise over run!
        mnp <- (p2+p1)/2
        int = mnp[2]-mnp[1]*slp
        resmid = middleY-(middleX*slp + int)
        -sum(c(dnorm(c(resleft),mean=0,sd=sdY,log=T),dnorm(c(resmid),mean=0,sd=sdY,log=T),dnorm(c(resright),mean=0,sd=sdY,log=T)))
      } else {Inf}
  }

llike.bs.new <- function(p,X,Y)
  {
    p1 <- p[1:2]
    p2 <- p[3:4]
    sd <- p[5]
#    sd=0.2
    if ((p1[1]<p2[1])&(sd>0)) {ok <- T} else {ok <- F}# points not swapped, sd ok
    if (ok & ((p1[1]>min(X))&(p2[1]<max(X)))) {ok <- T} else {ok <- F}# points not at end of line
#    if (ok & (min(c(p1[2],p2[2])>=0))) {ok <- T} else {ok <- F}
#    if (ok & (max(c(p1[2],p2[2])<=1))) {ok <- T} else {ok <- F}

    slp <- (p2-p1)[2]/(p2-p1)[1] #rise over run!
    mnp <- (p2+p1)/2
    int = mnp[2]-mnp[1]*slp
    if (slp==0) {ok <- FALSE} else {ok <- TRUE}

    if (ok)
      {
        leftX <- X[which(X<p1[1])]
        leftY <- Y[which(X<p1[1])]
        rightX <- X[which(X>p2[1])]
        rightY <- Y[which(X>p2[1])]
        middleX <- X[which((X>=p1[1])&(X<=p2[1]))]
        middleY <- Y[which((X>=p1[1])&(X<=p2[1]))]
        
        resleft <- p1[2]-leftY
        resright <- p2[2]-rightY
        
        resmidY = middleY-(middleX*slp + int)
        resmidX = ((middleY - int)/slp) - middleX
        -sum(c(dnorm(c(resleft),mean=0,sd=sd,log=T),
               dnorm(c(resmidY),mean=0,sd=sd,log=T),dnorm(resmidX,mean=0,sd=sd,log=T),
               dnorm(c(resright),mean=0,sd=sd,log=T)))
      } else {Inf}
  }

llike.flt.old <- function(p,X,Y)
  {
    if (p[2]>0)
      {
        -sum(dnorm(Y-p[1],mean=0,sd=p[2],log=T))
      } else Inf
  }


llike.flt <- function(p,X,Y)
  {
    sdY <- sd(Y)
    -sum(dnorm(Y-p[1],mean=0,sd=sdY,log=T))
  }

extract.bs <- function(bs.par)
  {
    slope <- (bs.par[2]-bs.par[4])/(bs.par[1]-bs.par[3])
    mid <- mean(bs.par[c(1,3)])
    c(slope=slope,mid=mid)
  }



################non broken stick models

sig <-  function(X,Y) 
  {

#    initial values
    p <- rep(0,3)
    p[1] <- 17 #c
    p[2] <- 10 #w
    p[3] <- 1 #sd
    optim(p,llike.sig,X=X,Y=Y)
  }

llike.sig <- function(p, X, Y)
  {
    c <- p[1]
    w <- p[2]
    sd <- p[3]
    if (sd>0) {
      resid <- unlist(Y - (1+tanh(2*(X-c)/w))/2)
      -sum(dnorm(resid,mean=0,sd=sd,log=T))
    } else
    {
      -Inf
    }
  }

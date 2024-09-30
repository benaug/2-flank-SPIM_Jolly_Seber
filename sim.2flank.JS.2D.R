e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.2flank.JS.2D <- function(lambda.y1=NA,gamma=NA,n.year=NA,
                       phi=phi,p0.B=NA,p0.L=NA,p0.R=NA,sigma=NA,X=NA,buff=buff,K=NA,
                       K2D=NA,J.cams=NA,n.fixed=NA,sigma.move=NULL,seed=NULL){
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  #Population dynamics
  N <- rep(NA,n.year)
  N.recruit <- N.survive <- ER <- rep(NA,n.year-1)
  N[1] <- rpois(1,lambda.y1)
  
  #Easiest to increase dimension of z as we simulate bc size not known in advance.
  z <- matrix(0,N[1],n.year)
  z[1:N[1],1] <- 1
  for(g in 2:n.year){
    #Simulate recruits
    ER[g-1] <- N[g-1]*gamma[g-1]
    N.recruit[g-1] <- rpois(1,ER[g-1])
    #add recruits to z
    z.dim.old <- nrow(z)
    z <- rbind(z,matrix(0,nrow=N.recruit[g-1],ncol=n.year))
    z[(z.dim.old+1):(z.dim.old+N.recruit[g-1]),g] <- 1
    
    #Simulate survival
    idx <- which(z[,g-1]==1)
    z[idx,g] <- rbinom(length(idx),1,phi[g-1])
    N.survive[g-1] <- sum(z[,g-1]==1&z[,g]==1)
    N[g] <- N.recruit[g-1]+N.survive[g-1]
  }
  
  if(any(N.recruit+N.survive!=N[2:n.year]))stop("Simulation bug")
  if(any(colSums(z)!=N))stop("Simulation bug")
  
  #detection
  #get maximal x and y extent across yearly grids plus buffer
  xlim  <-  c(max(unlist(lapply(X,function(x){min(x[,1])}))),max(unlist(lapply(X,function(x){max(x[,1])})))) + c(-buff,buff)
  ylim  <-  c(max(unlist(lapply(X,function(x){min(x[,2])}))),max(unlist(lapply(X,function(x){max(x[,2])})))) + c(-buff,buff)
  J  <-  unlist(lapply(X,nrow)) #extract number of traps per year
  J.max <- max(J)
  K.max <- max(K)
  
  #simulate activity centers - fixed through time
  N.super <- nrow(z)
  library(truncnorm)
  if(!is.null(sigma.move)){
    print("simulating mobile ACs")
    s <- array(NA,dim=c(N.super,n.year,2))
    s[,1,] <- cbind(runif(N.super, xlim[1],xlim[2]), runif(N.super,ylim[1],ylim[2]))
    for(g in 2:n.year){
      s[,g,1] <- rtruncnorm(N.super,s[,g-1,1],sd=sigma.move,a=xlim[1],b=xlim[2])
      s[,g,2] <- rtruncnorm(N.super,s[,g-1,2],sd=sigma.move,a=ylim[1],b=ylim[2])
    }
  }else{
    print("simulating fixed ACs (provide sigma.move for mobile)")
    s <- cbind(runif(N.super, xlim[1],xlim[2]), runif(N.super,ylim[1],ylim[2]))
  }
  
  kern <- pd.B <- pd.L <- pd.R <- array(0,dim=c(N.super,n.year,J.max))
  y.B <- y.L <- y.R <- array(0,dim=c(N.super,n.year,J.max))
  for(g in 1:n.year){
    if(!is.null(sigma.move)){
      D <- e2dist(s[,g,],X[[g]])
    }else{
      D <- e2dist(s,X[[g]])
    }
    kern[,g,1:J[g]] <- exp(-D*D/(2*sigma[g]*sigma[g]))
    pd.B[,g,1:J[g]] <- p0.B[g]*kern[,g,1:J[g]]
    pd.L[,g,1:J[g]] <- p0.L[g]*kern[,g,1:J[g]]
    pd.R[,g,1:J[g]] <- p0.R[g]*kern[,g,1:J[g]]
    for(i in 1:N.super){
      if(z[i,g]==1){
        for(j in 1:J[g]){
          if(J.cams[g,j]==1){
            y.L[i,g,j] <- rbinom(1,K2D[g,j],pd.L[i,g,j])
            y.R[i,g,j] <- rbinom(1,K2D[g,j],pd.R[i,g,j])
          }else{ #2 cams
            #P(A or B)=P(A)+P(B)-P(A and B)
            y.L[i,g,j] <- rbinom(1,K2D[g,j],2*pd.L[i,g,j]-pd.L[i,g,j]^2) #single side p. two chances for capture with 2 cameras
            y.R[i,g,j] <- rbinom(1,K2D[g,j],2*pd.R[i,g,j]-pd.R[i,g,j]^2)
            y.B[i,g,j] <- rbinom(1,K2D[g,j],pd.B[i,g,j])
          }
        }
      }
    }
  }
  
  #store true data for model buildling/debugging
  truth <- list(y.B=y.B,y.L=y.L,y.R=y.R,N=N,N.recruit=N.recruit,N.survive=N.survive,z=z,s=s)
  
  #Process flanks
  ID.B <- which(apply(y.B,1,sum)>0)
  n.B <- length(ID.B)
  ID.L <- which(rowSums(y.L)>0)
  ID.R <- which(rowSums(y.R)>0)
  ID.L <- setdiff(ID.L,ID.B)
  ID.R <- setdiff(ID.R,ID.B)
  ID.L <- c(ID.B,ID.L)
  ID.R <- c(ID.B,ID.R)
  if(!is.na(n.fixed)){
    if(n.fixed<n.B){
      print("More both-side captured individuals than n.fixed. Changing n.fixed to n.B")
      n.fixed <- n.B
    }
    #add uncaptured but known left and right flanks
    ID.L <- sort(unique(c(ID.L,1:n.fixed)))
    ID.R <- sort(unique(c(ID.R,1:n.fixed)))
  }
  
  Nknown <- length(ID.B)
  n <- sum(apply(y.B+y.L+y.R,1,sum)>0)
  
  #remove uncaptured individuals, put ID.B at top of y.L and y.R
  y.B.obs <- y.B[ID.B,,]
  y.L.obs <- y.L[ID.L,,]
  y.R.obs <- y.R[ID.R,,]
  n.B <- length(ID.B)
  n.L <- length(ID.L)
  n.R <- length(ID.R)
  if(n.B==1){
    y.B.obs <- array(y.B.obs,dim=c(1,g,J,K))
  }
  if(n.L==1){
    y.L.obs <- array(y.L.obs,dim=c(1,g,J,K))
  }
  if(n.R==1){
    y.R.obs <- array(y.R.obs,dim=c(1,g,J,K))
  }
  if(is.na(n.fixed)){
    n.fixed <- n.B
  }
  
  out <-list(N=N,N.recruit=N.recruit,N.survive=N.survive,X=X,J=J,K=K,K2D=K2D,J.cams=J.cams,n.year=n.year,
             xlim=xlim,ylim=ylim,y.B.obs=y.B.obs,y.L.obs=y.L.obs,y.R.obs=y.R.obs,
             y.B=y.B,y.L=y.L,y.R=y.R,s=s,n=n,ID.B=ID.B,ID.L=ID.L,ID.R=ID.R,n.fixed=n.fixed,
             truth=truth,seed=seed)
  return(out)
}

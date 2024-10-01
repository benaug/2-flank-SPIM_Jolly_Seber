e2dist<-function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.2flank.JS.mobileAC.2D <- function(data=data,M=M,n.fixed=NA,initTrue=FALSE){
  X <- data$X
  J <- unlist(lapply(X,nrow)) #traps per year
  J.max <- max(J)
  K.max <- data$K
  n.year <- dim(data$y.L.obs)[2]
  
  n.B <- nrow(data$y.B.obs)
  if(is.na(n.fixed)|(n.B==n.fixed)){
    if(n.B>0){
      print("Assuming the only known flank matches are for the n.B both side captured individuals.")
      n.fixed <- n.B
    }else{
      print("Assuming no known flank linkages because no both side captures and n.fixed not supplied or is 0.")
      n.fixed <- 0
    }
  }else{
    if(n.fixed<n.B)stop("n.fixed should be greater than n.B. We should know flank matches for n.B individuals.")
    n.diff <- n.fixed-n.B
    print("Assuming 1:n.fixed individual flank matches are known instead of 1:n.B because n.fixed supplied and is larger than n.B.")
  }
  
  if(initTrue){#initialize flanks to true matches?
    if(!all(c("y.B","y.L","y.R")%in%names(data)))stop("To init to truth, true data y.B, y.L, and y.R must be in data object.")
    y.B <- data$y.B.obs
    y.L <- data$y.L.obs
    y.R <- data$y.R.obs
    ID.B.true <- data$ID.B
    ID.L.true <- data$ID.L
    ID.R.true <- data$ID.R
    n.B <- nrow(y.B)
    n.L <- nrow(y.L)
    n.R <- nrow(y.R)
    
    y.true <- array(0,dim=c(M,n.year,J.max,3))
    ID.L <- rep(NA,n.L)
    ID.R <- rep(NA,n.R)
    if(n.fixed>0){
      if(n.B>0){
        y.true[1:n.B,,,1] <- y.B[1:n.B,,]
      }
      y.true[1:n.fixed,,,2] <- y.L[1:n.fixed,,]
      y.true[1:n.fixed,,,3] <- y.R[1:n.fixed,,]
      ID.L[1:n.fixed] <- ID.R[1:n.fixed] <- 1:n.fixed
    }
    if(n.fixed>0){
      ID.remain <- unique(c(ID.L.true[-c(1:n.fixed)],ID.R.true[-c(1:n.fixed)]))
    }else{
      ID.remain <- unique(c(ID.L.true,ID.R.true))
    }
    for(i in 1:length(ID.remain)){
      idx <- which(ID.L.true==ID.remain[i])
      if(length(idx)>0){
        y.true[i+n.fixed,,,2] <- y.L[idx,,]
        ID.L[idx] <- i+n.fixed
        
      }
      idx <- which(ID.R.true==ID.remain[i])
      if(length(idx)>0){
        y.true[i+n.fixed,,,3] <- y.R[idx,,]
        ID.R[idx] <- i+n.fixed
      }
    }
  }else{ #initialize without knowing truth
    #observed data, not using known data to initialize
    y.B <- data$y.B.obs
    y.L <- data$y.L.obs
    y.R <- data$y.R.obs
    #caps per type
    n.B <- nrow(y.B)
    n.L <- nrow(y.L)
    n.R <- nrow(y.R)
    y.true <- array(0,dim=c(M,n.year,J.max,3))
    #initialize lefts and rights to "both side" order
    if(n.B>0){
      y.true[1:n.B,,,1] <- y.B
    }
    #initialized unknown flank matches using spatial locations of captures
    #reformatting data so I do not have to modify LRmatch()
    #paste along year dimension
    y.L2D <- array(0,dim=c(n.L,sum(J)))
    y.R2D <- array(0,dim=c(n.R,sum(J)))
    X.all <- matrix(NA,sum(J),2)
    idx1 <- 1
    for(g in 1:n.year){
      idx2 <- idx1 + J[g] - 1
      y.L2D[,idx1:idx2] <- y.L[,g,]
      y.R2D[,idx1:idx2] <- y.R[,g,]
      X.all[idx1:idx2,] <- X[[g]]
      idx1 <- idx2 + 1
    }
    flank.init <- LRmatch(M=M, left=y.L2D, nleft=n.L-n.fixed,
                         right=y.R2D,nright=n.R-n.fixed, X=X.all, Nfixed=n.fixed)
    ID.L <- flank.init$ID.L
    ID.R <- flank.init$ID.R
    for(i in 1:length(ID.L)){
      y.true[ID.L[i],,,2] <- y.L[i,,]
    }
    for(i in 1:length(ID.R)){
      y.true[ID.R[i],,,3] <- y.R[i,,]
    }
  }
  #for simulated data, test to make sure initialization algorithm is correct
  if(all(c("y.B","y.L","y.R")%in%names(data))){
    if(n.B>0){
      for(i in 1:n.B){
        if(!all(y.true[i,,,1]==y.B[i,,]))stop("Error initializing y.B")
      }
    }
    for(i in 1:n.L){
      if(!all(y.true[ID.L[i],,,2]==y.L[i,,]))stop("Error initializing y.L")
    }
    for(i in 1:n.R){
      if(!all(y.true[ID.R[i],,,3]==y.R[i,,]))stop("Error initializing y.R")
    }
  }
  
  #JS stuff now
  y.nim <- apply(y.true,c(1,2,3),sum)
  y.nim2D <- apply(y.nim,c(1,2),sum)
  z.super.init <- 1*(rowSums(y.nim2D)>0)
  N.super.init <- sum(z.super.init)
  
  #initialize z, start with observed guys
  z.init <- 1*(y.nim>0)
  z.init <- matrix(0,M,n.year)
  z.start.init <- z.stop.init <- rep(NA,M)
  #initialize detected guys
  detected.inds <- which(rowSums(y.nim2D)>0)
  for(i in detected.inds){
    det.idx <- which(y.nim2D[i,]>0)
    if(length(det.idx)>0){ #some all 0 histories if init from truth
      z.start.init[i] <- min(det.idx)
      z.stop.init[i] <- max(det.idx)
      z.init[i,z.start.init[i]:z.stop.init[i]] <- 1
    }
  }
  #initialize undetected guys
  undetected.inds <- which(rowSums(y.nim2D)==0)
  for(i in undetected.inds){
    start <- sample(1:n.year,1) #random recruit year
    if(start<(n.year-1)){ #random death year
      stop <- sample(start:n.year,1)
    }else{ #unless recruited one year before end
      stop <- n.year
    }
    z.init[i,start:stop] <- 1
    z.start.init[i] <- start
    z.stop.init[i] <- stop
  }
  z.obs <- 1*(rowSums(y.nim)>0) #indicator for "ever observed"
  
  #initialize N structures from z.init
  N.init <- colSums(z.init[z.super.init==1,])
  N.survive.init <- N.recruit.init <- rep(NA,n.year-1)
  for(g in 2:n.year){
    N.survive.init[g-1] <- sum(z.init[,g-1]==1&z.init[,g]==1&z.super.init==1)
    N.recruit.init[g-1] <- N.init[g]-N.survive.init[g-1]
  }

  #remaining SCR stuff to initialize
  #put X in ragged array
  X.nim <- array(0,dim=c(n.year,J.max,2))
  for(g in 1:n.year){
    X.nim[g,1:J[g],1:2] <- X[[g]]
  }
  
  #pull out state space with buffer around maximal trap dimensions
  xlim <- data$xlim
  ylim <- data$ylim
  s.init <- array(NA,dim=c(M,n.year,2))
  for(g in 1:n.year){
    s.init[,g,]=cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
    idx <- which(rowSums(y.nim[,g,])>0) #switch for those actually caught
    for(i in idx){
      trps <- matrix(X.nim[g,which(y.nim[i,g,]>0),],ncol=2,byrow=FALSE)
      if(nrow(trps)>1){
        s.init[i,g,] <- c(mean(trps[,1]),mean(trps[,2]))
      }else{
        s.init[i,g,] <- trps
      }
    }
  }
  
  return(list(y.true=y.true,ID.L=ID.L,ID.R=ID.R,n.B=n.B,n.L=n.L,n.R=n.R,
              xlim=xlim,ylim=ylim,X.nim=X.nim,
              s=s.init,N=N.init,N.survive=N.survive.init,N.recruit=N.recruit.init,
              N.super=N.super.init,z.start=z.start.init,z.stop=z.stop.init,
              z=z.init,z.super=z.super.init,M=M))
}
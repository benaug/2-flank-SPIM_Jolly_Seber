LRmatch <-  function(M, left, nleft, right,nright, X, Nfixed){
  # This function takes a left and right data set and ad hocly associates unknown right guys
  # with unknown left guys to minimize the total separation distance.
  #Unknown left and right only guys can't be captured both guys because then they would be known.
  #Function modified from previous version by Andy Royle
  
  #assumes more left flanks than right flanks. If not, swaps them to initialize and swaps back at the end
  swapback <- FALSE
  if(nright>nleft){
    nright.tmp <- nright
    nright <- nleft
    nleft <- nright.tmp
    right.tmp <- right
    right <- left
    left <- right.tmp
    swapback <- TRUE
  }

  # function needs to spit out initial ID and activity centers
  #put in flag for cases where left or right are all known
  idkpl <- (Nfixed+1):dim(left)[1]
  idkpr <- (Nfixed+1):dim(right)[1]
  #Extract unknown individuals
  ld <- left[idkpl,]
  rd <-right[idkpr,]
  X <- as.matrix(X[,1:2])

   #matrices to store initial activity centers
  sbar.left <- matrix(NA,nrow=nleft,ncol=2)
  sbar.right <- matrix(NA,nrow=nright,ncol=2)

  #record average location of capture
  for(i in 1:nleft){
    if(sum(ld[i,])>0){  # should always be satisfied... unless you didn't remove uncaptured individuals from simulated data set.
      pos.traps <- (1:nrow(X))[ld[i,]>0] #Where was this guy left captured?
      pos.traps <- rep(pos.traps, ld[i,][ld[i,]>0]) #rep each trap by number of captures
      sbar.left[i,] <- apply(matrix(X[pos.traps,],ncol=2,byrow=FALSE),2,mean) #record average location of capture
    }
    else{
      sbar.left[i,] <- X[sample(1:nrow(X),1),] #if not captured, pick random location.  But this shouldn't happen.
    }
  }
  for(i in 1:nright){
    if(sum(rd[i,])>0){
      pos.traps <- (1:nrow(X))[rd[i,]>0]
      pos.traps <- rep(pos.traps, rd[i,][rd[i,]>0])
      sbar.right[i,] <- apply(matrix(X[pos.traps,],ncol=2,byrow=FALSE),2,mean)
    }
    else{
      sbar.right[i,] <- X[sample(1:nrow(X),1),]
    }
  }

  D <- e2dist(sbar.right,sbar.left)
  # optimization problem here is to put a 1 in each row such that sum of all distance is small
  ID.R2L <- sample(1:nleft, nright)
  Q <- sum( D[cbind(1:nright,ID.R2L)]  )

  if(nleft > nright){ #Should always happen
    for(loop in 1:20){
      for(i in 1:nrow(D)){
        # if there are unused left guys then try to make a swap there first
        notused <- (1:nleft)[is.na(match(1:nleft,ID.R2L))]
        curr.spot <- ID.R2L[i]
        Qtmp <- rep(NA,length(notused))
        for(k in 1:length(notused)){
          ID.R2L[i] <- notused[k]
          Qtmp[k] <- sum(  D[cbind(1:nright,ID.R2L)] )
        }
        if(min(Qtmp) < Q ){
          # Make the swap
          swap.in <- Qtmp==min(Qtmp)
          ID.R2L[i] <- notused[Qtmp==min(Qtmp)][1]  # Just use the first one
          Q <- min(Qtmp)
        }
        else{
          ID.R2L[i] <- curr.spot
        }
        #cat("new Q: ", Q, fill=TRUE)
      }
    }

  }
  ## for the last loop no other point could change with the available points.
  for(loop in 1:20){
    for(i in 1:nrow(D)){
      curr.spot <- ID.R2L[i]
      Qtmp <- rep(NA,length(ID.R2L))  # loop over EACH other
      for(k in 1:length(ID.R2L)){
        ID.R2L[i] <-  ID.R2L[k]
        ID.R2L[k] <- curr.spot
        Qtmp[k] <- sum(  D[cbind(1:nright,ID.R2L)] )
        ID.R2L[k] <- ID.R2L[i] # set it back to where it was after computing criterion
        ID.R2L[i] <- curr.spot
      }
      if(min(Qtmp) < Q ){
        # Make the swap
        which <- (1:length(ID.R2L))[Qtmp==min(Qtmp)][1]
        swap.in <- ID.R2L[which]
        ID.R2L[i] <- swap.in
        ID.R2L[which]<- curr.spot
        Q <- min(Qtmp)
      }
      else{  # Make sure to put it back if no swap was made
        ID.R2L[i] <- curr.spot
      }
      #cat("new Q: ", Q, fill=TRUE)
    }
  }
  #Randomly match lefts to boths not captured
  ID.L <- sample((Nfixed+1):M,nleft)
  #Translate new left IDs to right IDs
  ID.R <- ID.L[ID.R2L]
  #Add back to known individuals
  if(Nfixed>0){
    ID.L <- c((1:Nfixed) , ID.L )
    ID.R <- c((1:Nfixed) , ID.R )
  }
  if(swapback){
    ID.L.tmp <- ID.L
    ID.L <- ID.R
    ID.R <- ID.L.tmp
  }
  
  return(list(ID.L=ID.L,ID.R=ID.R))
}

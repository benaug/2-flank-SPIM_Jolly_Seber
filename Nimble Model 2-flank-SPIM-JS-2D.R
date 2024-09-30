NimModel <- nimbleCode({
  ##Abundance##
  lambda.y1 ~ dunif(0,1000) #Expected starting population size
  N[1] ~ dpois(lambda.y1) #Realized starting population size
  for(g in 2:n.year){
    N[g] <- N.survive[g-1] + N.recruit[g-1] #yearly abundance
    #N.recruit and N.survive information also contained in z/z.start + z.stop
    #N.recruit has distributions assigned below, but survival distributions defined on z
  }
  N.super <- N[1] + sum(N.recruit[1:(n.year-1)]) #size of superpopulation
  
  #Recruitment
  gamma.fixed ~ dunif(0,2) #share per capita recruitment rate across years
  for(g in 1:(n.year-1)){
    gamma[g] <- gamma.fixed
    # gamma[g] ~ dunif(0,2)
    ER[g] <- N[g]*gamma[g] #yearly expected recruits
    N.recruit[g] ~ dpois(ER[g]) #yearly realized recruits
  }
  
  #Individual covariates
  for(i in 1:M){
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
  }
  
  #Survival (phi must have M x n.year - 1 dimension for custom updates to work)
  #without individual or year effects, use for loop to plug into phi[i,g]
  phi.fixed ~ dunif(0,1)
  for(i in 1:M){
    for(g in 1:(n.year-1)){ #plugging same individual phi's into each year for custom update
      phi[i,g] <- phi.fixed #individual by year survival
    }
    #survival likelihood (bernoulli) that only sums from z.start to z.stop
    z[i,1:n.year] ~ dSurvival(phi=phi[i,1:(n.year-1)],z.start=z.start[i],z.stop=z.stop[i])
  }
  
  ##Detection##
  sigma ~ dunif(0,10) #fixing sigma across years
  for(g in 1:n.year){
    #p0 varies by year
    p0.B[g] ~ dunif(0,1)
    p0.S[g] ~ dunif(0,1)
    p0.L[g] <- p0.S[g]
    p0.R[g] <- p0.S[g]
    for(i in 1:M){ #only compute d2, pd.flank, and y when z.super[i]=1&z[i,g]=1
      d2[i,g,1:J[g]] <- GetD2(s=s[i,1:2],X=X[g,1:J[g],1:2],z=z[i,g],z.super=z.super[i])
      pd.B[i,g,1:J[g]] <- GetDetectionProb(p0=p0.B[g],sigma=sigma,d2=d2[i,g,1:J[g]],z=z[i,g],z.super=z.super[i])
      pd.L[i,g,1:J[g]] <- GetDetectionProb(p0=p0.L[g],sigma=sigma,d2=d2[i,g,1:J[g]],z=z[i,g],z.super=z.super[i])
      pd.R[i,g,1:J[g]] <- GetDetectionProb(p0=p0.R[g],sigma=sigma,d2=d2[i,g,1:J[g]],z=z[i,g],z.super=z.super[i])
      y.B.true[i,g,1:J[g]] ~ dBernoulliVectorBoth(pd.B[i,g,1:J[g]],K2D=K2D[g,1:J[g]],J.cams=J.cams[g,1:J[g]],
                                                  z=z[i,g],z.super=z.super[i])
      y.L.true[i,g,1:J[g]] ~ dBernoulliVectorSingle(pd.L[i,g,1:J[g]],K2D=K2D[g,1:J[g]],J.cams=J.cams[g,1:J[g]],
                                                    z=z[i,g],z.super=z.super[i])
      y.R.true[i,g,1:J[g]] ~ dBernoulliVectorSingle(pd.R[i,g,1:J[g]],K2D=K2D[g,1:J[g]],J.cams=J.cams[g,1:J[g]],
                                                    z=z[i,g],z.super=z.super[i])
    }
  }
  IDdummy <- IDdummyfun(ID.L=ID.L[1:n.L],ID.R=ID.R[1:n.R])
})
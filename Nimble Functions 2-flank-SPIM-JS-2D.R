IDdummyfun <- nimbleFunction(
  run = function(ID.L=double(1),ID.R=double(1)){ 
    returnType(double(0))
    return(0)
  }
)
GetD2 <- nimbleFunction(
  run = function(s=double(1),X=double(2), z=double(0), z.super=double(0)){ 
    returnType(double(1))
    J <- nimDim(X)[1]
    if(z.super==0 | z.super==1&z==0){
      return(rep(0,J)) #skip calculation if not is superpop, or in superpop, but not alive in this primary occasion
    }
    if(z==1){ #otherwise calculate
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      return(d2)
    }
  }
)
#make pd 3D here and in PMF functions if adding occasion effects to model. 
#Must modify custom ID updates to reflect this change.
GetDetectionProb <- nimbleFunction(
  run = function(d2=double(1), p0=double(0),sigma=double(0), z=double(0), z.super=double(0)){ 
    returnType(double(1))
    J <- nimDim(d2)[1]
    if(z.super==0 | z.super==1&z==0){
      return(rep(0,J)) #skip calculation if not is superpop, or in superpop, but not alive in this primary occasion
    }
    if(z==1){ #otherwise calculate
      pd <-p0*exp(-d2/(2*sigma^2))
      return(pd)
    }
  }
)
dBernoulliVectorBoth <- nimbleFunction(
  run = function(x = double(1), pd=double(1), K2D=double(1), J.cams = double(1), z = double(0), z.super=double(0),log = integer(0)) {
    returnType(double(0))
    if(z.super==0 | z.super==1&z==0){#skip calculation if not is superpop, or in superpop, but not alive in this primary occasion
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      J <- nimDim(K2D)[1]
      logProb <- 0
      for(j in 1:J){
        if(J.cams[j]==2){
          logProb <- logProb + dbinom(x[j],prob=pd[j],size=K2D[j],log=TRUE)
        }
      }
      return(logProb)
    }
  }
)
dBernoulliVectorSingle <- nimbleFunction(
  run = function(x = double(1), pd=double(1), K2D=double(1), J.cams=double(1), z = double(0), z.super=double(0),log = integer(0)) {
    returnType(double(0))
    if(z.super==0 | z.super==1&z==0){#skip calculation if not is superpop, or in superpop, but not alive in this primary occasion
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      J <- nimDim(K2D)[1]
      logProb <- 0
      for(j in 1:J){
        if(J.cams[j]==1){
          logProb <- logProb + dbinom(x[j],prob=pd[j],size=K2D[j],log=TRUE)
        }else{
          logProb <- logProb + dbinom(x[j],prob=2*pd[j]-pd[j]^2,size=K2D[j],log=TRUE)
        }
      }
      return(logProb)
    }
  }
)
#dummy function to make this work in parallel
rBernoulliVectorBoth <- nimbleFunction(
  run = function(n=integer(0),pd=double(1), K2D=double(1), J.cams=double(1), z = double(0), z.super=double(0)) {
    returnType(double(1))
    J <- nimDim(K2D)[1]
    y.true <- rep(0,J)
    return(y.true)
  }
)
#dummy function to make this work in parallel
rBernoulliVectorSingle <- nimbleFunction(
  run = function(n=integer(0),pd=double(1), K2D=double(1), J.cams=double(1), z = double(0), z.super=double(0)) {
    returnType(double(1))
    J <- nimDim(K2D)[1]
    y.true <- rep(0,J)
    return(y.true)
  }
)

#this is used to restrict likelihood evaluation to only the primary occasions relevant for survival for each individual
dSurvival <- nimbleFunction(
  run = function(x = double(1), phi = double(1), z.start = double(0), z.stop = double(0), z.super = double(0),
                 log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    if(z.super==1){
      n.primary <- length(phi)+1
      #extract first and last survival event primary occasions
      surv.start <- z.start+1
      surv.stop <- z.stop+1 #count death events, first z[i,]=0
      if(surv.start <= n.primary){ #if surv.start beyond last primary occasion, no survival events, logProb=0
        if(surv.stop > n.primary){ #but can't survive past n.primary
          surv.stop <- n.primary 
        }
        for(g in surv.start:surv.stop){ #sum logprob over survival event primary occasions
          logProb <- logProb + dbinom(x[g], size = 1, p = phi[g-1], log = TRUE)
        }
      }
    }
    return(logProb)
  }
)

#make dummy random vector generator to make nimble happy
rSurvival <- nimbleFunction(
  run = function(n = integer(0),phi = double(1), z.start = double(0), z.stop = double(0), z.super = double(0)) {
    returnType(double(1))
    n.primary <- length(phi)+1
    return(rep(0,n.primary))
  }
)

#all z updates live here
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    M <- control$M
    J <- control$J
    z.super.ups <- control$z.super.ups
    n.primary <- control$n.primary
    z.nodes <- control$z.nodes
    d2.nodes <- control$d2.nodes
    y.B.nodes <- control$y.B.nodes
    y.L.nodes <- control$y.L.nodes
    y.R.nodes <- control$y.R.nodes
    pd.B.nodes <- control$pd.B.nodes
    pd.L.nodes <- control$pd.L.nodes
    pd.R.nodes <- control$pd.R.nodes
    N.nodes <- control$N.nodes
    ER.nodes <- control$ER.nodes
    N.survive.nodes <- control$N.survive.nodes
    N.recruit.nodes <- control$N.recruit.nodes
    calcNodes <- control$calcNodes
  },
  run = function(){
    #compute y2D and z.obs for current data
    y2D <- matrix(0,M,n.primary)
    z.obs <- rep(0,M)
    for(i in 1:M){
      for(g in 1:n.primary){
        y2D[i,g] <- sum(model$y.B.true[i,g,]) + sum(model$y.L.true[i,g,]) + sum(model$y.R.true[i,g,])
      }
      if(sum(y2D[i,])>0){
        z.obs[i] <- 1
      }
    }
    #precompute entry counts
    entry.counts.curr <- rep(0,n.primary+1)
    for(i in 1:M){
      if(model$z.super[i]==1){
        entry.counts.curr[model$z.start[i]] <- entry.counts.curr[model$z.start[i]]+1
      }else{
        entry.counts.curr[n.primary+1] <- entry.counts.curr[n.primary+1]+1
      }
    }
    #1) Detected guy updates: z.start, z.stop
    # 1a) z start update (z.stop update below): Gibbs, compute full conditional
    for(i in 1:M){
      if(z.obs[i]==1&y2D[i,1]==0){ #for detected guys, skip if observed 1st primary occasion
        z.curr <- model$z[i,]
        z.start.curr <- model$z.start[i]
        N.curr <- model$N
        N.recruit.curr <- model$N.recruit
        dets <- which(y2D[i,]>0)
        first.det <- min(dets)
        lp.start <- rep(-Inf,n.primary)
        i.idx <- seq(i,M*n.primary,M) #used to reference correct y and pd nodes
        #only y, d2, and pd nodes before first detection can change across z.start candidates
        y.idx <- i.idx[1:(first.det-1)]
        #remove focal individual from entry counts once. The candidate-specific part of the
        #multinomial coefficient is then just log(entry.counts.minus[g]+1)
        entry.counts.minus <- entry.counts.curr
        entry.counts.minus[z.start.curr] <- entry.counts.minus[z.start.curr] - 1
        #all z.start > 1 candidates have the same N[1], so only calculate that logProb once
        lp.N1.not1 <- 0
        # Here, we are looping over all valid recruit dates and storing the logProb for each
        for(g in 1:first.det){ #must be recruited in primary occasion with first detection or before
          z.start.prop <- g
          model$z.start[i] <<- z.start.prop
          z.prop <- rep(0,n.primary)
          z.prop[g:first.det] <- 1 #must be alive until first detection
          if(first.det < n.primary){
            z.prop[(first.det+1):n.primary] <- z.curr[(first.det+1):n.primary] #fill in remaining current z values, keeping death event the same
          }
          model$z[i,] <<- z.prop
          #update N, N.recruit, N.survive. These individuals always in superpopulation
          #1) Update N
          model$N <<- N.curr - z.curr + z.prop
          #2) Update N.recruit
          model$N.recruit <<- N.recruit.curr #set back to original first
          if(z.start.curr > 1){ #if wasn't in pop in primary occasion 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in primary occasion 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          #only ER nodes before first detection can change across z.start candidates
          #model$calculate(ER.nodes) #update ER when N updated
          model$calculate(ER.nodes[1:(first.det-1)])
          #only d2 and pd nodes before first detection can change across z.start candidates
          #model$calculate(d2.nodes[i.idx]); model$calculate(pd.B.nodes[i.idx]); model$calculate(pd.L.nodes[i.idx]); model$calculate(pd.R.nodes[i.idx]) #update d2 and pd nodes when a z changes
          model$calculate(d2.nodes[y.idx])
          model$calculate(pd.B.nodes[y.idx])
          model$calculate(pd.L.nodes[y.idx])
          model$calculate(pd.R.nodes[y.idx])
          #get these logProbs
          #lp.N1 <- model$calculate(N.nodes[1])
          #there are only two possible N[1] values: z.start=1 and z.start>1
          if(g==1){
            lp.N1 <- model$calculate(N.nodes[1])
          }else{
            if(g==2){
              lp.N1.not1 <- model$calculate(N.nodes[1])
            }
            lp.N1 <- lp.N1.not1
          }
          #only recruitment likelihoods before first detection can change
          #lp.N.recruit <- model$calculate(N.recruit.nodes)
          lp.N.recruit <- model$calculate(N.recruit.nodes[1:(first.det-1)])
          #only observation likelihoods before first detection can change
          #lp.y <- model$calculate(y.B.nodes[i.idx]) + model$calculate(y.L.nodes[i.idx]) + model$calculate(y.R.nodes[i.idx])
          lp.y <- model$calculate(y.B.nodes[y.idx]) + model$calculate(y.L.nodes[y.idx]) + model$calculate(y.R.nodes[y.idx])
          lp.surv <- model$calculate(z.nodes[i])
          # Add the full multinomial coefficient prior log-prob for this proposed configuration
          #entry.counts.prop <- entry.counts.curr
          #z.super always 1 for detected guys
          #entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          #entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          #lp.prior <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop + 1)))
          #after removing this individual, all multinomial coefficient terms common across
          #candidates cancel, leaving log(n.g+1) for candidate entry class g
          lp.prior <- log(entry.counts.minus[g]+1)
          lp.start[g] <- lp.N1 + lp.N.recruit + lp.y + lp.surv + lp.prior
        }
        maxlp <- max(lp.start) #deal with overflow
        prop.probs <- exp(lp.start-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)
        z.start.prop <- rcat(1,prop.probs)
        model$z.start[i] <<- z.start.curr #set back to original
        if(model$z.start[i]!=z.start.prop){#if proposal is same as current, no need to replace anything
          model$z.start[i] <<- z.start.prop
          z.prop <- rep(0,n.primary)
          z.prop[model$z.start[i]:first.det] <- 1 #must be alive until first detection
          if(first.det < n.primary){
            z.prop[(first.det+1):n.primary] <- z.curr[(first.det+1):n.primary] #fill in remaining current z values, keeping death event the same
          }
          model$z[i,] <<- z.prop
          model$N <<- N.curr - z.curr + z.prop
          model$N.recruit <<- N.recruit.curr #set back to original first
          if(z.start.curr > 1){ #if wasn't in pop in primary occasion 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- N.recruit.curr[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in primary occasion 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- N.recruit.curr[z.start.prop-1] + 1
          }
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          #model$calculate(ER.nodes) #update ER when N updated
          model$calculate(ER.nodes[1:(first.det-1)])
          #model$calculate(d2.nodes[i.idx]); model$calculate(pd.B.nodes[i.idx]); model$calculate(pd.L.nodes[i.idx]); model$calculate(pd.R.nodes[i.idx]) #update d2 and pd nodes
          model$calculate(d2.nodes[y.idx])
          model$calculate(pd.B.nodes[y.idx])
          model$calculate(pd.L.nodes[y.idx])
          model$calculate(pd.R.nodes[y.idx])
          #update these logProbs
          #model$calculate(y.B.nodes[i.idx]) + model$calculate(y.L.nodes[i.idx]) + model$calculate(y.R.nodes[i.idx])
          model$calculate(y.B.nodes[y.idx])
          model$calculate(y.L.nodes[y.idx])
          model$calculate(y.R.nodes[y.idx])
          model$calculate(N.nodes[1])
          #model$calculate(N.recruit.nodes)
          model$calculate(N.recruit.nodes[1:(first.det-1)])
          model$calculate(z.nodes[i])
          mvSaved["z.start",1][i] <<- model[["z.start"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["N.recruit",1] <<- model[["N.recruit"]]
          mvSaved["ER",1] <<- model[["ER"]]
          #only d2 and pd nodes before first detection changed
          #for(g in 1:n.primary){
          for(g2 in 1:(first.det-1)){
            for(j in 1:J[g2]){
              mvSaved["d2",1][i,g2,j] <<- model[["d2"]][i,g2,j]
              mvSaved["pd.B",1][i,g2,j] <<- model[["pd.B"]][i,g2,j]
              mvSaved["pd.L",1][i,g2,j] <<- model[["pd.L"]][i,g2,j]
              mvSaved["pd.R",1][i,g2,j] <<- model[["pd.R"]][i,g2,j]
            }
          }
          #recompute entry counts
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          entry.counts.curr <- entry.counts.prop
        }else{
          model[["z.start"]][i] <<- mvSaved["z.start",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["N.recruit"]] <<- mvSaved["N.recruit",1]
          model[["ER"]] <<- mvSaved["ER",1]
          #only d2 and pd nodes before first detection changed
          #for(g in 1:n.primary){
          for(g2 in 1:(first.det-1)){
            for(j in 1:J[g2]){
              model[["d2"]][i,g2,j] <<- mvSaved["d2",1][i,g2,j]
              model[["pd.B"]][i,g2,j] <<- mvSaved["pd.B",1][i,g2,j]
              model[["pd.L"]][i,g2,j] <<- mvSaved["pd.L",1][i,g2,j]
              model[["pd.R"]][i,g2,j] <<- mvSaved["pd.R",1][i,g2,j]
            }
          }
          #set these logProbs back
          model$calculate(N.nodes[1])
          #model$calculate(N.recruit.nodes)
          model$calculate(N.recruit.nodes[1:(first.det-1)])
          #model$calculate(y.B.nodes[i.idx]) + model$calculate(y.L.nodes[i.idx]) + model$calculate(y.R.nodes[i.idx])
          model$calculate(y.B.nodes[y.idx])
          model$calculate(y.L.nodes[y.idx])
          model$calculate(y.R.nodes[y.idx])
          model$calculate(z.nodes[i])
        }
      }
    }
    #1b) z stop update (z.start update above): Gibbs, compute full conditional
    for(i in 1:M){
      if(z.obs[i]==1&y2D[i,n.primary]==0){ #for detected guys, skip if observed in final primary occasion
        z.curr <- model$z[i,]
        z.stop.curr <- model$z.stop[i]
        N.curr <- model$N
        dets <- which(y2D[i,]>0)
        last.det <- max(dets)
        lp.stop <- rep(-Inf,n.primary)
        i.idx <- seq(i,M*n.primary,M) #used to reference correct y and pd nodes
        #only y, d2, and pd nodes after last detection can change across z.stop candidates
        y.idx <- i.idx[(last.det+1):n.primary]
        #Here, we are looping over all valid z.stops and storing the logProb for each
        #entry prior does not change in z.stop update
        for(g in (last.det):n.primary){ #can't die on or before primary occasion of last detection
          model$z.stop[i] <<- g
          z.prop <- rep(0,n.primary)
          z.prop[last.det:g] <- 1 #must be alive between last detection and this z.stop
          z.prop[1:(last.det)] <- z.curr[1:(last.det)] #fill in remaining current z values, keeping death event the same
          model$z[i,] <<- z.prop
          #update N, number of recruits does not change going backwards
          model$N <<- N.curr - z.curr + z.prop
          #only ER nodes after last detection can change across z.stop candidates
          #model$calculate(ER.nodes) #update ER when N updated
          if(last.det < n.primary-1){
            model$calculate(ER.nodes[(last.det+1):(n.primary-1)])
          }
          #only d2 and pd nodes after last detection can change across z.stop candidates
          #model$calculate(d2.nodes[i.idx]); model$calculate(pd.B.nodes[i.idx]); model$calculate(pd.L.nodes[i.idx]); model$calculate(pd.R.nodes[i.idx]) #update d2 and pd nodes when z changes
          model$calculate(d2.nodes[y.idx])
          model$calculate(pd.B.nodes[y.idx])
          model$calculate(pd.L.nodes[y.idx])
          model$calculate(pd.R.nodes[y.idx])
          #get these logProbs
          #lp.N1 <- model$calculate(N.nodes[1])
          #N[1] cannot change in a z.stop update, so this is constant across candidates and cancels
          #lp.N.recruit <- model$calculate(N.recruit.nodes)
          if(last.det < n.primary-1){
            lp.N.recruit <- model$calculate(N.recruit.nodes[(last.det+1):(n.primary-1)])
          }else{
            lp.N.recruit <- 0
          }
          #lp.y <- model$calculate(y.B.nodes[i.idx]) + model$calculate(y.L.nodes[i.idx]) + model$calculate(y.R.nodes[i.idx])
          lp.y <- model$calculate(y.B.nodes[y.idx]) + model$calculate(y.L.nodes[y.idx]) + model$calculate(y.R.nodes[y.idx])
          lp.surv <- model$calculate(z.nodes[i])
          #no prior term, z.stop update does not change it
          #lp.stop[g] <- lp.N1 + lp.N.recruit + lp.y + lp.surv
          lp.stop[g] <- lp.N.recruit + lp.y + lp.surv
        }
        maxlp <- max(lp.stop) #deal with overflow
        prop.probs <- exp(lp.stop-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)
        z.stop.prop <- rcat(1,prop.probs)
        model$z.stop[i] <<- z.stop.curr #set back to original
        if(model$z.stop[i]!=z.stop.prop){#if proposal differs from current
          model$z.stop[i] <<- z.stop.prop
          z.prop <- rep(0,n.primary)
          z.prop[last.det:model$z.stop[i]] <- 1 #must be alive between last detection and this z.stop
          z.prop[1:(last.det)] <- z.curr[1:(last.det)] #fill in remaining current z values, keeping death event the same
          model$z[i,] <<- z.prop
          model$N <<- N.curr - z.curr + z.prop
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          #model$calculate(ER.nodes) #update ER when N updated
          if(last.det < n.primary-1){
            model$calculate(ER.nodes[(last.det+1):(n.primary-1)])
          }
          #model$calculate(d2.nodes[i.idx]); model$calculate(pd.B.nodes[i.idx]); model$calculate(pd.L.nodes[i.idx]); model$calculate(pd.R.nodes[i.idx]) #update d2 and pd nodes when z changes
          model$calculate(d2.nodes[y.idx])
          model$calculate(pd.B.nodes[y.idx])
          model$calculate(pd.L.nodes[y.idx])
          model$calculate(pd.R.nodes[y.idx])
          #update these logProbs
          #model$calculate(N.nodes[1]) #N[1] does not change in a z.stop update
          #model$calculate(N.recruit.nodes)
          if(last.det < n.primary-1){
            model$calculate(N.recruit.nodes[(last.det+1):(n.primary-1)])
          }
          #model$calculate(y.B.nodes[i.idx]) + model$calculate(y.L.nodes[i.idx]) + model$calculate(y.R.nodes[i.idx])
          model$calculate(y.B.nodes[y.idx])
          model$calculate(y.L.nodes[y.idx])
          model$calculate(y.R.nodes[y.idx])
          model$calculate(z.nodes[i])
          mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["ER",1] <<- model[["ER"]]
          #only d2 and pd nodes after last detection changed
          #for(g in 1:n.primary){
          for(g2 in (last.det+1):n.primary){
            for(j in 1:J[g2]){
              mvSaved["d2",1][i,g2,j] <<- model[["d2"]][i,g2,j]
              mvSaved["pd.B",1][i,g2,j] <<- model[["pd.B"]][i,g2,j]
              mvSaved["pd.L",1][i,g2,j] <<- model[["pd.L"]][i,g2,j]
              mvSaved["pd.R",1][i,g2,j] <<- model[["pd.R"]][i,g2,j]
            }
          }
        }else{
          model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["ER"]] <<- mvSaved["ER",1]
          #only d2 and pd nodes after last detection changed
          #for(g in 1:n.primary){
          for(g2 in (last.det+1):n.primary){
            for(j in 1:J[g2]){
              model[["d2"]][i,g2,j] <<- mvSaved["d2",1][i,g2,j]
              model[["pd.B"]][i,g2,j] <<- mvSaved["pd.B",1][i,g2,j]
              model[["pd.L"]][i,g2,j] <<- mvSaved["pd.L",1][i,g2,j]
              model[["pd.R"]][i,g2,j] <<- mvSaved["pd.R",1][i,g2,j]
            }
          }
          #set these logProbs back
          #model$calculate(N.nodes[1]) #N[1] does not change in a z.stop update
          #model$calculate(N.recruit.nodes)
          if(last.det < n.primary-1){
            model$calculate(N.recruit.nodes[(last.det+1):(n.primary-1)])
          }
          #model$calculate(y.B.nodes[i.idx]) + model$calculate(y.L.nodes[i.idx]) + model$calculate(y.R.nodes[i.idx])
          model$calculate(y.B.nodes[y.idx])
          model$calculate(y.L.nodes[y.idx])
          model$calculate(y.R.nodes[y.idx])
          model$calculate(z.nodes[i])
        }
      }
    }
    #2) undetected guy update. Only if in the superpopulation. 
    # Metropolis-Hastings, Propose z vectors from priors
    #can try Gibbs instead. not sure which is more efficient in which conditions
    #MH is cheap, but wont mix as well, Gibbs is expensive, but mixes better.
    #entry counts current after z.start update
    for(i in 1:M){
      if(z.obs[i]==0&model$z.super[i]==1){
        z.curr <- model$z[i,]
        z.start.curr <- model$z.start[i]
        z.stop.curr <- model$z.stop[i]
        i.idx <- seq(i,M*n.primary,M) #used to reference correct y and pd nodes
        #get forwards recruitment probabilities
        recruit.probs.for <- c(model$lambda.y1,model$ER)
        recruit.probs.for <- recruit.probs.for/sum(recruit.probs.for)
        #track proposal probs
        #survival proposal probabilities cancel exactly with the survival likelihood because
        #the survival history is proposed from the same survival model used in the target
        log.prop.for <- log.prop.back <- 0
        #simulate recruitment
        z.start.prop <- rcat(1,recruit.probs.for)
        z.prop <- rep(0,n.primary)
        z.prop[z.start.prop] <- 1
        log.prop.for <- log.prop.for + log(recruit.probs.for[z.start.prop])
        #simulate survival
        #once the individual dies, remaining z's are already 0 so no more simulation is needed
        z.stop.prop <- z.start.prop
        if(z.start.prop < n.primary){#if you don't recruit in final primary occasion
          for(g in (z.start.prop+1):n.primary){
            if(z.prop[g-1]==1){
              z.prop[g] <- rbinom(1,1,model$phi[i,g-1]*z.prop[g-1])
              if(z.prop[g]==1){
                z.stop.prop <- g
              }
            }
          }
        }
        #if the proposed history is the current history, there is nothing to calculate or update
        if(z.start.prop!=z.start.curr|z.stop.prop!=z.stop.curr){
          #primary occasions where the proposed alive state differs from current
          z.changed <- which(z.prop!=z.curr)
          y.idx.changed <- i.idx[z.changed]
          #get initial logProbs
          lp.initial.entry <- model$getLogProb(N.nodes[1])
          lp.initial.entry <- lp.initial.entry + model$getLogProb(N.recruit.nodes)
          #only observation likelihoods in primary occasions where z changed
          lp.initial.y <- model$getLogProb(y.B.nodes[y.idx.changed]) + model$getLogProb(y.L.nodes[y.idx.changed]) + model$getLogProb(y.R.nodes[y.idx.changed])
          #lp.initial.surv <- model$getLogProb(z.nodes[i]) #cancels exactly with backwards survival proposal probability
          #log.prior.curr <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr + 1))) #full multinomial coefficient calculation replaced by exact ratio below
          model$z[i,] <<- z.prop
          model$z.start[i] <<- z.start.prop
          model$z.stop[i] <<- z.stop.prop
          #update N, N.recruit, N.survive only if individual is in superpopulation
          #1) Update N
          model$N <<- model$N - z.curr + z.prop
          #2) Update N.recruit
          if(z.start.curr > 1){ #if wasn't in pop in primary occasion 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in primary occasion 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          model$calculate(ER.nodes) #update ER when N updated
          #only d2 and pd nodes in primary occasions where z changed
          model$calculate(d2.nodes[y.idx.changed])
          model$calculate(pd.B.nodes[y.idx.changed])
          model$calculate(pd.L.nodes[y.idx.changed])
          model$calculate(pd.R.nodes[y.idx.changed])
          #get proposed logProbs
          lp.proposed.entry <- model$calculate(N.nodes[1])
          lp.proposed.entry <- lp.proposed.entry + model$calculate(N.recruit.nodes)
          #only observation likelihoods in primary occasions where z changed
          lp.proposed.y <- model$calculate(y.B.nodes[y.idx.changed]) + model$calculate(y.L.nodes[y.idx.changed]) + model$calculate(y.R.nodes[y.idx.changed])
          #lp.proposed.surv <- model$calculate(z.nodes[i]) #cancels exactly with forwards survival proposal probability
          # Full multinomial coefficient prior for proposed configuration
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          #log.prior.prop <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop + 1)))
          #In log.prior.prop-log.prior.curr, -lgamma(M+1) and all unchanged entry-class
          #lgamma terms cancel. If entry changes from class a to b, the remaining ratio is
          #log(n.b+1)-log(n.a). If entry does not change, the ratio is 0.
          if(z.start.prop!=z.start.curr){
            log.prior.ratio <- log(entry.counts.curr[z.start.prop]+1)-log(entry.counts.curr[z.start.curr])
          }else{
            log.prior.ratio <- 0
          }
          #get backwards proposal probs
          recruit.probs.back <- c(model$lambda.y1,model$ER)
          recruit.probs.back <- recruit.probs.back/sum(recruit.probs.back)
          log.prop.back <- log.prop.back + log(recruit.probs.back[z.start.curr])
          #survival proposal probabilities are not calculated because they cancel exactly
          #with the survival likelihood ratio in the MH ratio
          lp.initial.total <- lp.initial.entry + lp.initial.y
          lp.proposed.total <- lp.proposed.entry + lp.proposed.y
          #MH step
          log_MH_ratio <- (lp.proposed.total + log.prior.ratio + log.prop.back) - (lp.initial.total + log.prop.for)
          accept <- decide(log_MH_ratio)
          if(accept){
            #update survival logProb once for accepted history; it was not needed when evaluating MH ratio
            model$calculate(z.nodes[i])
            mvSaved["z.start",1][i] <<- model[["z.start"]][i]
            mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
            mvSaved["z",1][i,] <<- model[["z"]][i,]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["ER",1] <<- model[["ER"]]
            #only d2 and pd nodes in changed primary occasions need to be saved
            for(g3 in 1:length(z.changed)){
              g2 <- z.changed[g3]
              for(j in 1:J[g2]){
                mvSaved["d2",1][i,g2,j] <<- model[["d2"]][i,g2,j]
                mvSaved["pd.B",1][i,g2,j] <<- model[["pd.B"]][i,g2,j]
                mvSaved["pd.L",1][i,g2,j] <<- model[["pd.L"]][i,g2,j]
                mvSaved["pd.R",1][i,g2,j] <<- model[["pd.R"]][i,g2,j]
              }
            }
            entry.counts.curr <- entry.counts.prop
          }else{
            model[["z.start"]][i] <<- mvSaved["z.start",1][i]
            model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
            model[["z"]][i,] <<- mvSaved["z",1][i,]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["ER"]] <<- mvSaved["ER",1]
            #restore d2 and pd nodes only in primary occasions where z changed
            for(g3 in 1:length(z.changed)){
              g2 <- z.changed[g3]
              for(j in 1:J[g2]){
                model[["d2"]][i,g2,j] <<- mvSaved["d2",1][i,g2,j]
                model[["pd.B"]][i,g2,j] <<- mvSaved["pd.B",1][i,g2,j]
                model[["pd.L"]][i,g2,j] <<- mvSaved["pd.L",1][i,g2,j]
                model[["pd.R"]][i,g2,j] <<- mvSaved["pd.R",1][i,g2,j]
              }
            }
            #set these logProbs back
            #only observation logProbs in primary occasions where z changed need restoration
            model$calculate(y.B.nodes[y.idx.changed])
            model$calculate(y.L.nodes[y.idx.changed])
            model$calculate(y.R.nodes[y.idx.changed])
            model$calculate(N.nodes[1])
            model$calculate(N.recruit.nodes)
            #model$calculate(z.nodes[i]) #not needed because survival logProb was never recalculated for the proposal
          }
        }
      }
    }
    #3) update z.super: Metropolis-Hastings
    #entry counts current coming out of undetected ind update
    #make lists of currently on/off undetected guys once, then update after accepted proposals
    z.on <- rep(0,M)
    z.off <- rep(0,M)
    non.curr <- 0
    noff.curr <- 0
    for(i in 1:M){
      if(z.obs[i]==0){
        if(model$z.super[i]==1){
          non.curr <- non.curr+1
          z.on[non.curr] <- i
        }else{
          noff.curr <- noff.curr+1
          z.off[noff.curr] <- i
        }
      }
    }
    for(up in 1:z.super.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      if(updown==0){#subtract
        #find all z's currently on and undetected
        non.init <- non.curr
        if(non.init>0){
          pick.pos <- rcat(1,rep(1/non.init,non.init))
          pick <- z.on[pick.pos]
          z.start.curr <- model$z.start[pick]
          z.curr <- model$z[pick,]
          #p select on guy
          log.p.select.for <- log(1/non.init)
          #log multinomial coefficient prior
          #log.z.prior.for <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr+1)))
          pick.idx <- seq(pick,M*n.primary,M) #used to reference correct y and pd nodes
          #get initial logProbs (survival logProb does not change)
          lp.initial.N <- model$getLogProb(N.nodes[1])
          lp.initial.N.recruit <- model$getLogProb(N.recruit.nodes)
          lp.initial.y <- model$getLogProb(y.B.nodes[pick.idx]) + model$getLogProb(y.L.nodes[pick.idx]) + model$getLogProb(y.R.nodes[pick.idx])
          #lp.initial.surv <- model$getLogProb(z.nodes[pick]) #survival likelihood cancels exactly with reverse survival proposal probability
          #propose new N.super/z.super/z.start/z.stop
          model$N.super <<-  model$N.super - 1
          model$z.super[pick] <<- 0
          model$z.start[pick] <<- 0
          model$z.stop[pick] <<- 0
          model$z[pick,] <<- rep(0,n.primary)
          #update N, N.recruit, N.survive
          #1) Update N
          model$N <<- model$N - z.curr
          #2) Update N.recruit
          if(z.start.curr > 1){ #if wasn't in pop in primary occasion 1
            model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(d2.nodes[pick.idx]) #update pd nodes when z changes
          model$calculate(pd.B.nodes[pick.idx]) #update pd nodes when z changes
          model$calculate(pd.L.nodes[pick.idx]) #update pd nodes when z changes
          model$calculate(pd.R.nodes[pick.idx]) #update pd nodes when z changes
          #Reverse proposal probs
          recruit.probs.back <- c(model$lambda.y1, model$ER)
          recruit.probs.back <- recruit.probs.back / sum(recruit.probs.back)
          log.prop.back <- log(recruit.probs.back[z.start.curr])
          #survival proposal probability cancels exactly with current survival likelihood
          #if(z.start.curr < n.primary){
          #  for(g in (z.start.curr+1):n.primary){
          #    log.prop.back <- log.prop.back + dbinom(z.curr[g],1,model$phi[pick,g-1]*z.curr[g-1],log=TRUE)
          #  }
          #}
          #get proposed logProbs for N, N.recruit, and y
          lp.proposed.N <- model$calculate(N.nodes[1])
          lp.proposed.N.recruit <- model$calculate(N.recruit.nodes)
          lp.proposed.y <- model$calculate(y.B.nodes[pick.idx]) + model$calculate(y.L.nodes[pick.idx]) + model$calculate(y.R.nodes[pick.idx]) #will always be 0
          #lp.proposed.surv <- model$calculate(z.nodes[pick]) #survival likelihood cancels exactly with reverse survival proposal probability
          #survival target/proposal terms cancel exactly, so they are omitted from the MH totals
          #lp.initial.total <- lp.initial.N + lp.initial.y + lp.initial.N.recruit + lp.initial.surv
          #lp.proposed.total <- lp.proposed.N + lp.proposed.y + lp.proposed.N.recruit + lp.proposed.surv
          lp.initial.total <- lp.initial.N + lp.initial.y + lp.initial.N.recruit
          lp.proposed.total <- lp.proposed.N + lp.proposed.y + lp.proposed.N.recruit
          #backwards prior and select probs
          #move from class z.start.curr in z.super==0 to class g in z.super==1
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          entry.counts.prop[n.primary + 1] <- entry.counts.prop[n.primary + 1] + 1
          #p select off guy
          noff.back <- noff.curr+1
          log.p.select.back <- log(1/noff.back)
          #log multinomial coefficient prior
          #log.z.prior.back <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop+1)))
          #In log.z.prior.back-log.z.prior.for, -lgamma(M+1) and all unchanged
          #entry-class lgamma terms cancel. Only the old entry class and off class change:
          #lgamma(n.entry)-lgamma(n.entry+1)=-log(n.entry)
          #lgamma(n.off+2)-lgamma(n.off+1)=log(n.off+1)
          log.z.prior.ratio <- log(entry.counts.curr[n.primary+1]+1)-log(entry.counts.curr[z.start.curr])
          log.prop.for <- 0
          #MH step
          log_MH_ratio <- (lp.proposed.total + log.z.prior.ratio + log.p.select.back + log.prop.back) -
            (lp.initial.total + log.p.select.for + log.prop.for)
          accept <- decide(log_MH_ratio)
          if(accept){
            #survival logProb was omitted from the MH calculation because it cancels with the proposal;
            #calculate it once now to synchronize the accepted model state
            model$calculate(z.nodes[pick])
            mvSaved["z.start",1][pick] <<- model[["z.start"]][pick]
            mvSaved["z.stop",1][pick] <<- model[["z.stop"]][pick]
            mvSaved["z",1][pick,] <<- model[["z"]][pick,]
            mvSaved["z.super",1][pick] <<- model[["z.super"]][pick]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["N.super",1][1] <<- model[["N.super"]]
            mvSaved["ER",1] <<- model[["ER"]]
            for(g2 in 1:n.primary){
              for(j in 1:J[g2]){
                mvSaved["d2",1][pick,g2,j] <<- model[["d2"]][pick,g2,j]
                mvSaved["pd.B",1][pick,g2,j] <<- model[["pd.B"]][pick,g2,j]
                mvSaved["pd.L",1][pick,g2,j] <<- model[["pd.L"]][pick,g2,j]
                mvSaved["pd.R",1][pick,g2,j] <<- model[["pd.R"]][pick,g2,j]
              }
            }
            entry.counts.curr <- entry.counts.prop
            #move guy from on list to off list
            z.on[pick.pos] <- z.on[non.curr]
            z.on[non.curr] <- 0
            non.curr <- non.curr-1
            noff.curr <- noff.curr+1
            z.off[noff.curr] <- pick
          }else{
            model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
            model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
            model[["z"]][pick,] <<- mvSaved["z",1][pick,]
            model[["z.super"]][pick] <<- mvSaved["z.super",1][pick]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["N.super"]] <<- mvSaved["N.super",1][1]
            model[["ER"]] <<- mvSaved["ER",1]
            for(g2 in 1:n.primary){
              for(j in 1:J[g2]){
                model[["d2"]][pick,g2,j] <<- mvSaved["d2",1][pick,g2,j]
                model[["pd.B"]][pick,g2,j] <<- mvSaved["pd.B",1][pick,g2,j]
                model[["pd.L"]][pick,g2,j] <<- mvSaved["pd.L",1][pick,g2,j]
                model[["pd.R"]][pick,g2,j] <<- mvSaved["pd.R",1][pick,g2,j]
              }
            }
            #set these logProbs back
            model$calculate(y.B.nodes[pick.idx])
            model$calculate(y.L.nodes[pick.idx])
            model$calculate(y.R.nodes[pick.idx])
            #model$calculate(z.nodes[pick]) #not needed because survival logProb was never recalculated for the proposal
            model$calculate(N.nodes[1])
            model$calculate(N.recruit.nodes)
          }
        }
      }else{#add
        if(model$N.super[1] < M){ #cannot update if z.super maxed out. Need to raise M
          #find all z's currently off and undetected
          noff.init <- noff.curr
          if(noff.init>0){
            pick.pos <- rcat(1,rep(1/noff.init,noff.init)) #select one of these individuals
            pick <- z.off[pick.pos]
            pick.idx <- seq(pick,M*n.primary,M)
            #p select off guy
            log.p.select.for <- log(1/noff.init)
            #log multinomial coefficient prior
            #log.z.prior.for <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr+1)))
            #get initial logProbs (survival logProb does not change)
            lp.initial.N <- model$getLogProb(N.nodes[1])
            lp.initial.N.recruit <- model$getLogProb(N.recruit.nodes)
            lp.initial.y <- model$getLogProb(y.B.nodes[pick.idx]) + model$getLogProb(y.L.nodes[pick.idx]) + model$getLogProb(y.R.nodes[pick.idx]) #will always be 0
            #lp.initial.surv <- model$getLogProb(z.nodes[pick]) #survival likelihood cancels exactly with forward survival proposal probability
            # Propose new z.start for the new on individual
            recruit.probs.for <- c(model$lambda.y1, model$ER)
            recruit.probs.for <- recruit.probs.for / sum(recruit.probs.for)
            z.start.prop <- rcat(1, recruit.probs.for)  # propose entry cohort
            log.prop.for <- log(recruit.probs.for[z.start.prop])
            model$z.start[pick] <<- z.start.prop
            #Simulate survival path
            model$z[pick,] <<- 0 # initialize to 0
            model$z[pick, z.start.prop] <<- 1
            if(z.start.prop < n.primary){
              for(g in (z.start.prop+1):n.primary){
                model$z[pick, g] <<- rbinom(1, 1, model$phi[pick, g-1] * model$z[pick, g-1])
                #log.prop.for <- log.prop.for + dbinom(model$z[pick, g], 1, model$phi[pick, g-1] * model$z[pick, g-1], log=TRUE) #survival proposal probability cancels exactly with proposed survival likelihood
              }
            }
            #Update z.stop
            z.on.prop <- which(model$z[pick,] == 1)
            z.stop.prop <- max(z.on.prop)
            model$z.stop[pick] <<- z.stop.prop
            #propose new N/z
            model$N.super <<-  model$N.super + 1
            model$z.super[pick] <<- 1
            #update N, N.recruit, N.survive
            #1) Update N
            model$N <<- model$N + model$z[pick,]
            #2) Update N.recruit
            if(model$z.start[pick] > 1){ #if wasn't in pop in primary occasion 1
              model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
            }
            #3) Update N.survive
            model$N.survive <<- model$N[2:n.primary] - model$N.recruit #survivors are guys alive in primary occasion g-1 minus recruits in this primary occasion g
            model$calculate(ER.nodes) #update ER when N updated
            model$calculate(d2.nodes[pick.idx]) #update pd nodes when z changes
            model$calculate(pd.B.nodes[pick.idx]) #update pd nodes when z changes
            model$calculate(pd.L.nodes[pick.idx]) #update pd nodes when z changes
            model$calculate(pd.R.nodes[pick.idx]) #update pd nodes when z changes
            #get proposed logprobs for N and y
            lp.proposed.N <- model$calculate(N.nodes[1])
            lp.proposed.N.recruit <- model$calculate(N.recruit.nodes)
            lp.proposed.y <- model$calculate(y.B.nodes[pick.idx]) + model$calculate(y.L.nodes[pick.idx]) + model$calculate(y.R.nodes[pick.idx])
            #lp.proposed.surv <- model$calculate(z.nodes[pick]) #survival likelihood cancels exactly with forward survival proposal probability
            #survival target/proposal terms cancel exactly, so they are omitted from the MH totals
            #lp.initial.total <- lp.initial.N + lp.initial.y + lp.initial.N.recruit + lp.initial.surv
            #lp.proposed.total <- lp.proposed.N + lp.proposed.y + lp.proposed.N.recruit + lp.proposed.surv
            lp.initial.total <- lp.initial.N + lp.initial.y + lp.initial.N.recruit
            lp.proposed.total <- lp.proposed.N + lp.proposed.y + lp.proposed.N.recruit
            #backwards prior and select probs
            #move from class g in z.super==0 to class g in z.super==1
            entry.counts.prop <- entry.counts.curr
            entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
            entry.counts.prop[n.primary + 1] <- entry.counts.prop[n.primary + 1] - 1
            #p select on guy
            non.back <- non.curr+1
            log.p.select.back <- log(1/non.back)
            #log multinomial coefficient prior
            #log.z.prior.back <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop+1)))
            #In log.z.prior.back-log.z.prior.for, -lgamma(M+1) and all unchanged
            #entry-class lgamma terms cancel. Only the new entry class and off class change:
            #lgamma(n.entry+2)-lgamma(n.entry+1)=log(n.entry+1)
            #lgamma(n.off)-lgamma(n.off+1)=-log(n.off)
            log.z.prior.ratio <- log(entry.counts.curr[z.start.prop]+1)-log(entry.counts.curr[n.primary+1])
            log.prop.back <- 0
            #MH step
            log_MH_ratio <- (lp.proposed.total + log.z.prior.ratio + log.p.select.back + log.prop.back) -
              (lp.initial.total + log.p.select.for + log.prop.for)
            accept <- decide(log_MH_ratio)
            if(accept){
              #survival logProb was omitted from the MH calculation because it cancels with the proposal;
              #calculate it once now to synchronize the accepted model state
              model$calculate(z.nodes[pick])
              mvSaved["z.start",1][pick] <<- model[["z.start"]][pick]
              mvSaved["z.stop",1][pick] <<- model[["z.stop"]][pick]
              mvSaved["z",1][pick,] <<- model[["z"]][pick,]
              mvSaved["z.super",1][pick] <<- model[["z.super"]][pick]
              mvSaved["N",1] <<- model[["N"]]
              mvSaved["N.survive",1] <<- model[["N.survive"]]
              mvSaved["N.recruit",1] <<- model[["N.recruit"]]
              mvSaved["N.super",1][1] <<- model[["N.super"]]
              mvSaved["ER",1] <<- model[["ER"]]
              for(g2 in 1:n.primary){
                for(j in 1:J[g2]){
                  mvSaved["d2",1][pick,g2,j] <<- model[["d2"]][pick,g2,j]
                  mvSaved["pd.B",1][pick,g2,j] <<- model[["pd.B"]][pick,g2,j]
                  mvSaved["pd.L",1][pick,g2,j] <<- model[["pd.L"]][pick,g2,j]
                  mvSaved["pd.R",1][pick,g2,j] <<- model[["pd.R"]][pick,g2,j]
                }
              }
              entry.counts.curr <- entry.counts.prop
              #move guy from off list to on list
              z.off[pick.pos] <- z.off[noff.curr]
              z.off[noff.curr] <- 0
              noff.curr <- noff.curr-1
              non.curr <- non.curr+1
              z.on[non.curr] <- pick
            }else{
              model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
              model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
              model[["z"]][pick,] <<- mvSaved["z",1][pick,]
              model[["z.super"]][pick] <<- mvSaved["z.super",1][pick]
              model[["N"]] <<- mvSaved["N",1]
              model[["N.survive"]] <<- mvSaved["N.survive",1]
              model[["N.recruit"]] <<- mvSaved["N.recruit",1]
              model[["N.super"]] <<- mvSaved["N.super",1][1]
              model[["ER"]] <<- mvSaved["ER",1]
              for(g2 in 1:n.primary){
                for(j in 1:J[g2]){
                  model[["d2"]][pick,g2,j] <<- mvSaved["d2",1][pick,g2,j]
                  model[["pd.B"]][pick,g2,j] <<- mvSaved["pd.B",1][pick,g2,j]
                  model[["pd.L"]][pick,g2,j] <<- mvSaved["pd.L",1][pick,g2,j]
                  model[["pd.R"]][pick,g2,j] <<- mvSaved["pd.R",1][pick,g2,j]
                }
              }
              #set these logProbs back
              model$calculate(y.B.nodes[pick.idx])
              model$calculate(y.L.nodes[pick.idx])
              model$calculate(y.R.nodes[pick.idx])
              #model$calculate(z.nodes[pick]) #not needed because survival logProb was never recalculated for the proposal
              model$calculate(N.nodes[1])
              model$calculate(N.recruit.nodes)
            }
          }
        }
      }
    }
    #copy back to mySaved to update logProbs.
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

#Left flank update
IDLSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K2D <- control$K2D
    J.cams <- control$J.cams
    n.L <- control$n.L
    n.fixed <- control$n.fixed
    n.primary <- control$n.primary
    M <- control$M
    J <- control$J
    K <- control$K
    y.L.start <- control$y.L.start
    y.L.stop <- control$y.L.stop
    prop.scale <- control$prop.scale #can scale the distance proposal, but a value of 1 should be a good choice.
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    if(n.L>0){ #skip if no lefts
      z <- model$z
      z.super <- model$z.super
      s <- model$s
      y.L.true <- model$y.L.true
      ID.L <- model$ID.L
      pd.L <- model$pd.L
      sigma <- model$sigma
      #how to propose with z's not consistent? Just propose them that way and reject for now
      ll.y.L <- model$logProb_y.L.true
      ll.y.L.cand <- ll.y.L
      
      # z.start <- model$z.start
      # z.stop <- model$z.stop
      
      for(l in (n.fixed+1):n.L){
        y.L.true.cand <- y.L.true
        this.i <- ID.L[l]
        #other flanks we can't swap in to ID.L[l]
        conflict.idx <- which(y.L.start<model$z.start[this.i]|y.L.stop>model$z.stop[this.i])
        #convert to their individual number these flanks currently assigned to
        conflict.idx2 <- ID.L[conflict.idx]
        #individuals we can swap the focal flank to
        valid.swap <- 1*(model$z.start<=y.L.start[l]&model$z.stop>=y.L.stop[l]&model$z.super==1)
        #remove conflicts from above
        for(l2 in 1:length(conflict.idx2)){
          valid.swap[conflict.idx2[l2]] <- 0
        }
        if(n.fixed>0){
          valid.swap[1:n.fixed] <- 0 #cannot select an individual with fixed flanks
        }
        valid.swap[this.i] <- 0 #don't choose focal
        propprobs <- rep(0,M)
        for(i in 1:M){
          if(valid.swap[i]){
            d2 <- (s[this.i,1]-s[i,1])^2+(s[this.i,2]-s[i,2])^2
            propprobs[i] <- exp(-d2/(2*(prop.scale*sigma[1])^2)) #distance-based proposal. must reference index for sigma
          }
        }
        sum.propprobs <- sum(propprobs)
        if(sum.propprobs>0){ #prereject if no valid swaps
          propprobs <- propprobs/sum.propprobs
          cand.i <- rcat(1,prob=propprobs)
          
          cand.ID.idx <- which(ID.L==cand.i)#which ID.L index is this individual. Will be empty if not in ID.L
          for(g in 1:n.primary){ #"can't do math in >2 dimensions" loop
            for(j in 1:J[g]){
              y.L.true.cand[this.i,g,j] <- y.L.true[cand.i,g,j]
              y.L.true.cand[cand.i,g,j] <- y.L.true[this.i,g,j]
            }
          }
          #update ll.y.L
          for(g in 1:n.primary){
            ll.y.L.cand[this.i,g,1] <- dBernoulliVectorSingle(y.L.true.cand[this.i,g,1:J[g]],
                                                              pd.L[this.i,g,1:J[g]],
                                                              K2D=K2D[g,1:J[g]],
                                                              J.cams=J.cams[g,1:J[g]],
                                                              z=z[this.i,g],
                                                              z.super=z.super[this.i],log=TRUE)
            ll.y.L.cand[cand.i,g,1] <- dBernoulliVectorSingle(y.L.true.cand[cand.i,g,1:J[g]],
                                                              pd.L[cand.i,g,1:J[g]],
                                                              K2D=K2D[g,1:J[g]],
                                                              J.cams=J.cams[g,1:J[g]],
                                                              z=z[cand.i,g],
                                                              z.super=z.super[cand.i],log=TRUE)
          }
          
          #calculate backwards proposal probs
          #other flanks we can't swap in to cand.i
          conflict.idx <- which(y.L.start<model$z.start[cand.i]|y.L.stop>model$z.stop[cand.i])
          #convert to their individual number these flanks currently assigned to
          conflict.idx2 <- ID.L[conflict.idx]
          #individuals we can swap the focal flank to
          valid.swap <- 1*(model$z.start<=y.L.start[l]&model$z.stop>=y.L.stop[l]&model$z.super==1)
          #remove conflicts from above
          for(l2 in 1:length(conflict.idx2)){
            valid.swap[conflict.idx2[l2]] <- 0
          }
          if(n.fixed>0){
            valid.swap[1:n.fixed] <- 0 #cannot select an individual with fixed flanks
          }
          valid.swap[cand.i] <- 0 #don't choose focal
          backprobs <- rep(0,M)
          for(i in 1:M){
            if(valid.swap[i]){
              d2 <- (s[cand.i,1]-s[i,1])^2+(s[cand.i,2]-s[i,2])^2
              backprobs[i] <- exp(-d2/(2*(prop.scale*sigma[1])^2)) #distance-based proposal. must reference index for sigma
            }
          }
          backprobs <- backprobs/sum(backprobs)
          # lp_initial <- sum(ll.y.L[this.i,,]) + sum(ll.y.L[cand.i,,])
          # lp_proposed <- sum(ll.y.L.cand[this.i,,]) + sum(ll.y.L.cand[cand.i,,])
          #"can't do math in >2 dimensions" loop
          lp_initial <- lp_proposed <- 0
          for(g in 1:n.primary){
            lp_initial <- lp_initial + ll.y.L[this.i,g,1] + ll.y.L[cand.i,g,1]
            lp_proposed <- lp_proposed + ll.y.L.cand[this.i,g,1] + ll.y.L.cand[cand.i,g,1]
          }
          log_MH_ratio <- (lp_proposed+log(backprobs[this.i])) - (lp_initial+log(propprobs[cand.i]))
          accept <- decide(log_MH_ratio)
          if(accept){
            for(g in 1:n.primary){ #"can't do math in >2 dimensions" loop
              for(j in 1:J[g]){
                y.L.true[this.i,g,j] <- y.L.true.cand[this.i,g,j]
                y.L.true[cand.i,g,j] <- y.L.true.cand[cand.i,g,j]
              }
              ll.y.L[this.i,g,1] <- ll.y.L.cand[this.i,g,1]
              ll.y.L[cand.i,g,1] <- ll.y.L.cand[cand.i,g,1]
            }
            #swap flank indices
            ID.L[l] <- cand.i
            #if cand.i was in ID.L, update this ID index
            tmp <- nimDim(cand.ID.idx)[1]
            if(tmp>0){
              ID.L[cand.ID.idx] <- this.i
            }
          }
        }
      }
      
      #put everything back into model$stuff
      model$y.L.true <<- y.L.true
      model$ID.L <<- ID.L
      model$calculate(calcNodes) #update logprob
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    }
  },
  methods = list( reset = function () {} )
)

#Right flank update
IDRSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    K2D <- control$K2D
    J.cams <- control$J.cams
    n.R <- control$n.R
    n.fixed <- control$n.fixed
    n.primary <- control$n.primary
    M <- control$M
    J <- control$J
    K <- control$K
    y.R.start <- control$y.R.start
    y.R.stop <- control$y.R.stop
    prop.scale <- control$prop.scale #can scale the distance proposal, but a value of 1 should be a good choice.
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    if(n.R>0){ #skip if no rights
      z <- model$z
      z.super <- model$z.super
      s <- model$s
      y.R.true <- model$y.R.true
      ID.R <- model$ID.R
      pd.R <- model$pd.R
      sigma <- model$sigma
      #how to propose with z's not consistent? Just propose them that way and reject for now
      ll.y.R <- model$logProb_y.R.true
      ll.y.R.cand <- ll.y.R
      for(l in (n.fixed+1):n.R){
        y.R.true.cand <- y.R.true
        this.i <- ID.R[l]
        #other flanks we can't swap in to ID.L[l]
        conflict.idx <- which(y.R.start<model$z.start[this.i]|y.R.stop>model$z.stop[this.i])
        #convert to their individual number these flanks currently assigned to
        conflict.idx2 <- ID.R[conflict.idx]
        #individuals we can swap the focal flank to
        valid.swap <- 1*(model$z.start<=y.R.start[l]&model$z.stop>=y.R.stop[l]&model$z.super==1)
        #remove conflicts from above
        for(l2 in 1:length(conflict.idx2)){
          valid.swap[conflict.idx2[l2]] <- 0
        }
        if(n.fixed>0){
          valid.swap[1:n.fixed] <- 0 #cannot select an individual with fixed flanks
        }
        valid.swap[this.i] <- 0 #don't choose focal
        propprobs <- rep(0,M)
        for(i in 1:M){
          if(valid.swap[i]){
            d2 <- (s[this.i,1]-s[i,1])^2+(s[this.i,2]-s[i,2])^2
            propprobs[i] <- exp(-d2/(2*(prop.scale*sigma[1])^2)) #distance-based proposal. must reference index for sigma
          }
        }
        sum.propprobs <- sum(propprobs)
        if(sum.propprobs>0){ #prereject if no valid swaps
          propprobs <- propprobs/sum.propprobs
          cand.i <- rcat(1,prob=propprobs)
          cand.ID.idx <- which(ID.R==cand.i)#which ID.R index is this individual. Will be empty if not in ID.R
          for(g in 1:n.primary){ #"can't do math in >2 dimensions" loop
            for(j in 1:J[g]){
              y.R.true.cand[this.i,g,j] <- y.R.true[cand.i,g,j]
              y.R.true.cand[cand.i,g,j] <- y.R.true[this.i,g,j]
            }
          }
          #update ll.y.R
          for(g in 1:n.primary){
            ll.y.R.cand[this.i,g,1] <- dBernoulliVectorSingle(y.R.true.cand[this.i,g,1:J[g]],
                                                              pd.R[this.i,g,1:J[g]],
                                                              K2D=K2D[g,1:J[g]],
                                                              J.cams=J.cams[g,1:J[g]],
                                                              z=z[this.i,g],
                                                              z.super=z.super[this.i],log=TRUE)
            ll.y.R.cand[cand.i,g,1] <- dBernoulliVectorSingle(y.R.true.cand[cand.i,g,1:J[g]],
                                                              pd.R[cand.i,g,1:J[g]],
                                                              K2D=K2D[g,1:J[g]],
                                                              J.cams=J.cams[g,1:J[g]],
                                                              z=z[cand.i,g],
                                                              z.super=z.super[cand.i],log=TRUE)
          }
          #calculate backwards proposal probs
          #other flanks we can't swap in to cand.i
          conflict.idx <- which(y.R.start<model$z.start[cand.i]|y.R.stop>model$z.stop[cand.i])
          #convert to their individual number these flanks currently assigned to
          conflict.idx2 <- ID.R[conflict.idx]
          #individuals we can swap the focal flank to
          valid.swap <- 1*(model$z.start<=y.R.start[l]&model$z.stop>=y.R.stop[l]&model$z.super==1)
          #remove conflicts from above
          for(l2 in 1:length(conflict.idx2)){
            valid.swap[conflict.idx2[l2]] <- 0
          }
          if(n.fixed>0){
            valid.swap[1:n.fixed] <- 0 #cannot select an individual with fixed flanks
          }
          valid.swap[cand.i] <- 0 #don't choose focal
          backprobs <- rep(0,M)
          for(i in 1:M){
            if(valid.swap[i]){
              d2 <- (s[cand.i,1]-s[i,1])^2+(s[cand.i,2]-s[i,2])^2
              backprobs[i] <- exp(-d2/(2*(prop.scale*sigma[1])^2)) #distance-based proposal. must reference index for sigma
            }
          }
          backprobs <- backprobs/sum(backprobs)
          # lp_initial <- sum(ll.y.R[this.i,,])+sum(ll.y.R[cand.i,,])
          # lp_proposed <- sum(ll.y.R.cand[this.i,,])+sum(ll.y.R.cand[cand.i,,])
          #"can't do math in >2 dimensions" loop
          lp_initial <- lp_proposed <- 0
          for(g in 1:n.primary){
            lp_initial <- lp_initial + ll.y.R[this.i,g,1] + ll.y.R[cand.i,g,1]
            lp_proposed <- lp_proposed + ll.y.R.cand[this.i,g,1] + ll.y.R.cand[cand.i,g,1]
          }
          log_MH_ratio <- (lp_proposed+log(backprobs[this.i])) - (lp_initial+log(propprobs[cand.i]))
          accept <- decide(log_MH_ratio)
          
          if(accept){
            for(g in 1:n.primary){ #"can't do math in >2 dimensions" loop
              for(j in 1:J[g]){
                y.R.true[this.i,g,j] <- y.R.true.cand[this.i,g,j]
                y.R.true[cand.i,g,j] <- y.R.true.cand[cand.i,g,j]
              }
              ll.y.R[this.i,g,1] <- ll.y.R.cand[this.i,g,1]
              ll.y.R[cand.i,g,1] <- ll.y.R.cand[cand.i,g,1]
            }
            #swap flank indices
            ID.R[l] <- cand.i
            #if cand.i was in ID.R, update this ID index
            tmp <- nimDim(cand.ID.idx)[1]
            if(tmp>0){
              ID.R[cand.ID.idx] <- this.i
            }
          }
        }
      }
      #put everything back into model$stuff
      model$y.R.true <<- y.R.true
      model$ID.R <<- ID.R
      model$calculate(calcNodes) #update logprob
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    }
  },
  methods = list( reset = function () {} )
)

truncGammaPoisSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model,mvSaved,target,control){
    calcNodes <- model$getDependencies(target)
    upper <- model$getBound(target,"upper")
    tau <- control$tau
    if(target=="lambda.y1"){
      is.lambda <- TRUE
      is.fixed.gamma <- FALSE
      g <- 1
      n.recruit <- 1
    }else if(target=="gamma"){
      is.lambda <- FALSE
      is.fixed.gamma <- TRUE
      g <- 1
      recruitNodes <- grep("^N.recruit\\[",model$getNodeNames(stochOnly=TRUE),value=TRUE)
      n.recruit <- length(recruitNodes)
    }else{
      is.lambda <- FALSE
      is.fixed.gamma <- FALSE
      g <- as.integer(gsub("[^0-9]","",target))
      n.recruit <- 1
    }
  },
  run = function(){
    if(is.lambda){
      count <- model$N[1]
      rate <- 1
    }else if(is.fixed.gamma){
      count <- 0
      rate <- 0
      for(j in 1:n.recruit){
        count <- count+model$N.recruit[j]
        rate <- rate+model$N[j]*tau[j]
      }
    }else{
      count <- model$N.recruit[g]
      rate <- model$N[g]*tau[g]
    }
    if(rate>0){
      shape <- count+1
      p.upper <- pgamma(upper,shape=shape,rate=rate)
      model[[target]] <<- qgamma(runif(1,0,p.upper),shape=shape,rate=rate)
    }else{
      model[[target]] <<- runif(1,0,upper)
    }
    model$calculate(calcNodes)
    copy(from=model,to=mvSaved,row=1,nodes=calcNodes,logProb=TRUE)
  },
  methods=list(reset=function(){})
)
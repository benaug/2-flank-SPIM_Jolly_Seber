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
      return(rep(0,J)) #skip calculation if not is superpop, or in superpop, but not alive in this year
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
      return(rep(0,J)) #skip calculation if not is superpop, or in superpop, but not alive in this year
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
    if(z.super==0 | z.super==1&z==0){#skip calculation if not is superpop, or in superpop, but not alive in this year
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
    if(z.super==0 | z.super==1&z==0){#skip calculation if not is superpop, or in superpop, but not alive in this year
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

#this is used to restrict likelihood evaluation to only the years relevant for survival for each individual
dSurvival <- nimbleFunction(
  run = function(x = double(1), phi = double(1), z.start = double(0), z.stop = double(0), z.super = double(0),
                 log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    if(z.super==1){
      n.year <- length(phi)+1
      #extract first and last survival event years
      surv.start <- z.start+1
      surv.stop <- z.stop+1 #count death events, first z[i,]=0
      if(surv.start <= n.year){ #if surv.start beyond last year, no survival events, logProb=0
        if(surv.stop > n.year){ #but can't survive past n.year
          surv.stop <- n.year 
        }
        for(g in surv.start:surv.stop){ #sum logprob over survival event years
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
    n.year <- length(phi)
    return(rep(0,n.year))
  }
)

#all z updates live here
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    M <- control$M
    J <- control$J
    z.super.ups <- control$z.super.ups
    n.year <- control$n.year
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
    y2D <- matrix(0,M,n.year)
    z.obs <- rep(0,M)
    for(i in 1:M){
      for(g in 1:n.year){
        y2D[i,g] <- sum(model$y.B.true[i,g,]) + sum(model$y.L.true[i,g,]) + sum(model$y.R.true[i,g,])
      }
      if(sum(y2D[i,])>0){
        z.obs[i] <- 1
      }
    }
    #precompute entry counts
    entry.counts.curr <- rep(0, n.year+1)
    for(g in 1:n.year){
      entry.counts.curr[g] <- sum(model$z.start==g & model$z.super==1)
    }
    entry.counts.curr[n.year + 1] <- sum(model$z.super==0)
    
    #1) Detected guy updates: z.start, z.stop
    # 1a) z start update (z.stop update below): Gibbs, compute full conditional
    for(i in 1:M){
      if(z.obs[i]==1&y2D[i,1]==0){ #for detected guys, skip if observed 1st year
        z.curr <- model$z[i,]
        z.start.curr <- model$z.start[i]
        N.curr <- model$N
        N.recruit.curr <- model$N.recruit
        dets <- which(y2D[i,]>0)
        first.det <- min(dets)
        lp.start <- rep(-Inf,n.year)
        i.idx <- seq(i,M*n.year,M) #used to reference correct y and pd nodes
        for(g in 1:first.det){ #must be recruited in year with first detection or before
          z.start.prop <- g
          model$z.start[i] <<- z.start.prop
          z.prop <- rep(0,n.year)
          z.prop[g:first.det] <- 1 #must be alive until first detection
          if(first.det < n.year){
            z.prop[(first.det+1):n.year] <- z.curr[(first.det+1):n.year] #fill in remaining current z values, keeping death event the same
          }
          model$z[i,] <<- z.prop
          
          #update N, N.recruit, N.survive. These individuals always in superpopulation
          #1) Update N
          model$N <<- N.curr - z.curr + z.prop
          #2) Update N.recruit
          model$N.recruit <<- N.recruit.curr #set back to original first
          if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- N.recruit.curr[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- N.recruit.curr[z.start.prop-1] + 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.year]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(d2.nodes[i.idx]) #update d2 nodes when a z changes
          model$calculate(pd.B.nodes[i.idx]) #update pd nodes when a z changes
          model$calculate(pd.L.nodes[i.idx]) #update pd nodes when a z changes
          model$calculate(pd.R.nodes[i.idx]) #update pd nodes when a z changes
          #get these logProbs
          lp.N1 <- model$calculate(N.nodes[1])
          lp.N.recruit <- model$calculate(N.recruit.nodes)
          lp.y <- model$calculate(y.B.nodes[i.idx]) + 
            model$calculate(y.L.nodes[i.idx]) + 
            model$calculate(y.R.nodes[i.idx])
          lp.surv <- model$calculate(z.nodes[i])
          # Add the full multinomial coefficient prior log-prob for this proposed configuration
          entry.counts.prop <- entry.counts.curr
          #z.super always 1 for detected guys
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          lp.prior <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop + 1)))
          lp.start[g] <- lp.N1 + lp.N.recruit + lp.y + lp.surv + lp.prior
        }
        maxlp <- max(lp.start) #deal with overflow
        prop.probs <- exp(lp.start-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)
        
        z.start.prop <- rcat(1,prop.probs)
        model$z.start[i] <<- z.start.curr #set back to original
        
        if(model$z.start[i]!=z.start.prop){#if proposal is same as current, no need to replace anything
          model$z.start[i] <<- z.start.prop
          z.prop <- rep(0,n.year)
          z.prop[model$z.start[i]:first.det] <- 1 #must be alive until first detection
          if(first.det < n.year){
            z.prop[(first.det+1):n.year] <- z.curr[(first.det+1):n.year] #fill in remaining current z values, keeping death event the same
          }
          model$z[i,] <<- z.prop
          model$N <<- N.curr - z.curr + z.prop
          model$N.recruit <<- N.recruit.curr #set back to original first
          if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- N.recruit.curr[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- N.recruit.curr[z.start.prop-1] + 1
          }
          model$N.survive <<- model$N[2:n.year]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(d2.nodes[i.idx]) #update d2 nodes when a z changes
          model$calculate(pd.B.nodes[i.idx]) #update pd nodes
          model$calculate(pd.L.nodes[i.idx]) #update pd nodes
          model$calculate(pd.R.nodes[i.idx]) #update pd nodes
          #update these logProbs
          model$calculate(y.B.nodes[i.idx])
          model$calculate(y.L.nodes[i.idx])
          model$calculate(y.R.nodes[i.idx])
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
          model$calculate(z.nodes[i])
          mvSaved["z.start",1][i] <<- model[["z.start"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["N.recruit",1] <<- model[["N.recruit"]]
          mvSaved["ER",1] <<- model[["ER"]]
          for(g in 1:n.year){
            for(j in 1:J[g]){
              mvSaved["d2",1][i,g,j] <<- model[["d2"]][i,g,j]
              mvSaved["pd.B",1][i,g,j] <<- model[["pd.B"]][i,g,j]
              mvSaved["pd.L",1][i,g,j] <<- model[["pd.L"]][i,g,j]
              mvSaved["pd.R",1][i,g,j] <<- model[["pd.R"]][i,g,j]
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
          for(g in 1:n.year){
            for(j in 1:J[g]){
              model[["d2"]][i,g,j] <<- mvSaved["d2",1][i,g,j]
              model[["pd.B"]][i,g,j] <<- mvSaved["pd.B",1][i,g,j]
              model[["pd.L"]][i,g,j] <<- mvSaved["pd.L",1][i,g,j]
              model[["pd.R"]][i,g,j] <<- mvSaved["pd.R",1][i,g,j]
            }
          }
          #set these logProbs back
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
          model$calculate(y.B.nodes[i.idx])
          model$calculate(y.L.nodes[i.idx])
          model$calculate(y.R.nodes[i.idx])
          model$calculate(z.nodes[i])
        }
      }
    }
    
    #1b) z stop update (z.start update above): Gibbs, compute full conditional
    for(i in 1:M){
      if(z.obs[i]==1&y2D[i,n.year]==0){ #for detected guys, skip if observed in final year
        z.curr <- model$z[i,]
        z.stop.curr <- model$z.stop[i]
        N.curr <- model$N
        dets <- which(y2D[i,]>0)
        last.det <- max(dets)
        lp.stop <- rep(-Inf,n.year)
        i.idx <- seq(i,M*n.year,M) #used to reference correct y and pd nodes
        for(g in (last.det):n.year){ #can't die on or before year of last detection
          model$z.stop[i] <<- g
          z.prop <- rep(0,n.year)
          z.prop[last.det:g] <- 1 #must be alive between last detection and this z.stop
          z.prop[1:(last.det)] <- z.curr[1:(last.det)] #fill in remaining current z values, keeping death event the same
          model$z[i,] <<- z.prop
          #update N, number of recruits does not change going backwards
          model$N <<- N.curr - z.curr + z.prop
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(d2.nodes[i.idx]) #update d2 nodes when a z changes
          model$calculate(pd.B.nodes[i.idx]) #update pd nodes when z changes
          model$calculate(pd.L.nodes[i.idx]) #update pd nodes when z changes
          model$calculate(pd.R.nodes[i.idx]) #update pd nodes when z changes
          #get these logProbs
          lp.N1 <- model$calculate(N.nodes[1])
          lp.N.recruit <- model$calculate(N.recruit.nodes)
          lp.y <- model$calculate(y.B.nodes[i.idx]) + 
            model$calculate(y.L.nodes[i.idx]) + 
            model$calculate(y.R.nodes[i.idx])
          lp.surv <- model$calculate(z.nodes[i])
          #no prior term, z.stop update does not change it
          lp.stop[g] <- lp.N1 + lp.N.recruit + lp.y + lp.surv
        }
        maxlp <- max(lp.stop) #deal with overflow
        prop.probs <- exp(lp.stop-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)
        z.stop.prop <- rcat(1,prop.probs)
        model$z.stop[i] <<- z.stop.curr #set back to original
        if(model$z.stop[i]!=z.stop.prop){#if proposal differs from current
          model$z.stop[i] <<- z.stop.prop
          z.prop <- rep(0,n.year)
          z.prop[last.det:model$z.stop[i]] <- 1 #must be alive between last detection and this z.stop
          z.prop[1:(last.det)] <- z.curr[1:(last.det)] #fill in remaining current z values, keeping death event the same
          model$z[i,] <<- z.prop
          model$N <<- N.curr - z.curr + z.prop
          model$N.survive <<- model$N[2:n.year]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(d2.nodes[i.idx]) #update d2 nodes when a z changes
          model$calculate(pd.B.nodes[i.idx]) #update pd nodes when z changes
          model$calculate(pd.L.nodes[i.idx]) #update pd nodes when z changes
          model$calculate(pd.R.nodes[i.idx]) #update pd nodes when z changes
          #update these logProbs
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
          model$calculate(y.B.nodes[i.idx])
          model$calculate(y.L.nodes[i.idx])
          model$calculate(y.R.nodes[i.idx])
          model$calculate(z.nodes[i])
          mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["ER",1] <<- model[["ER"]]
          for(g in 1:n.year){
            for(j in 1:J[g]){
              mvSaved["d2",1][i,g,j] <<- model[["d2"]][i,g,j]
              mvSaved["pd.B",1][i,g,j] <<- model[["pd.B"]][i,g,j]
              mvSaved["pd.L",1][i,g,j] <<- model[["pd.L"]][i,g,j]
              mvSaved["pd.R",1][i,g,j] <<- model[["pd.R"]][i,g,j]
            }
          }
        }else{
          model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["ER"]] <<- mvSaved["ER",1]
          for(g in 1:n.year){
            for(j in 1:J[g]){
              model[["d2"]][i,g,j] <<- mvSaved["d2",1][i,g,j]
              model[["pd.B"]][i,g,j] <<- mvSaved["pd.B",1][i,g,j]
              model[["pd.L"]][i,g,j] <<- mvSaved["pd.L",1][i,g,j]
              model[["pd.R"]][i,g,j] <<- mvSaved["pd.R",1][i,g,j]
            }
          }
          #set these logProbs back
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
          model$calculate(y.B.nodes[i.idx])
          model$calculate(y.L.nodes[i.idx])
          model$calculate(y.R.nodes[i.idx])
          model$calculate(z.nodes[i])
        }
      }
    }
    #2) undetected guy update. Only if in the superpopulation. 
    # Metropolis-Hastings, Propose z vectors from priors
    #entry counts current after z.start update
    for(i in 1:M){
      if(z.obs[i]==0&model$z.super[i]==1){
        z.curr <- model$z[i,]
        z.start.curr <- model$z.start[i]
        z.stop.curr <- model$z.stop[i]
        i.idx <- seq(i,M*n.year,M) #used to reference correct y and pd nodes
        #get forwards recruitment probabilities
        recruit.probs.for <- c(model$lambda.y1,model$ER)
        recruit.probs.for <- recruit.probs.for/sum(recruit.probs.for)
        #get initial logProbs
        lp.initial.entry <- model$getLogProb(N.nodes[1])
        lp.initial.entry <- lp.initial.entry + model$getLogProb(N.recruit.nodes)
        lp.initial.y <- model$getLogProb(y.B.nodes[i.idx]) + 
          model$getLogProb(y.L.nodes[i.idx]) + 
          model$getLogProb(y.R.nodes[i.idx])
        lp.initial.surv <- model$getLogProb(z.nodes[i])
        log.prior.curr <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr + 1)))
        
        #track proposal probs - survival is symmetric, but not recruitment and detection
        log.prop.for <- log.prop.back <- 0
        
        #simulate recruitment
        z.start.prop <- rcat(1,recruit.probs.for)
        z.prop <- rep(0,n.year)
        z.prop[z.start.prop] <- 1
        log.prop.for <- log.prop.for + log(recruit.probs.for[z.start.prop])
        
        #simulate survival
        if(z.start.prop < n.year){#if you don't recruit in final year
          for(g in (z.start.prop+1):n.year){
            z.prop[g] <- rbinom(1,1,model$phi[i,g-1]*z.prop[g-1])
            log.prop.for <- log.prop.for + dbinom(z.prop[g],1,model$phi[i,g-1]*z.prop[g-1],log=TRUE)
          }
        }
        z.on.prop <- which(z.prop==1)
        z.stop.prop <- max(z.on.prop)
        model$z[i,] <<- z.prop
        model$z.start[i] <<- z.start.prop
        model$z.stop[i] <<- z.stop.prop
        
        #update N, N.recruit, N.survive only if individual is in superpopulation
        #1) Update N
        model$N <<- model$N - z.curr + z.prop
        #2) Update N.recruit
        if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
          model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
        }
        if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
          model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
        }
        #3) Update N.survive
        model$N.survive <<- model$N[2:n.year]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
        
        model$calculate(ER.nodes) #update ER when N updated
        model$calculate(d2.nodes[i.idx]) #update d2 nodes when a z changes
        model$calculate(pd.B.nodes[i.idx]) #update pd nodes when z changes
        model$calculate(pd.L.nodes[i.idx]) #update pd nodes when z changes
        model$calculate(pd.R.nodes[i.idx]) #update pd nodes when z changes
        #get proposed logProbs
        lp.proposed.entry <- model$calculate(N.nodes[1])
        lp.proposed.entry <- lp.proposed.entry + model$calculate(N.recruit.nodes)
        lp.proposed.y <- model$calculate(y.B.nodes[i.idx]) + 
          model$calculate(y.L.nodes[i.idx]) + 
          model$calculate(y.R.nodes[i.idx])
        lp.proposed.surv <- model$calculate(z.nodes[i])
        
        # Full multinomial coefficient prior for proposed configuration
        entry.counts.prop <- entry.counts.curr
        entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
        entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
        log.prior.prop <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop + 1)))
        
        #get backwards proposal probs
        recruit.probs.back <- c(model$lambda.y1,model$ER)
        recruit.probs.back <- recruit.probs.back/sum(recruit.probs.back)
        log.prop.back <- log.prop.back + log(recruit.probs.back[z.start.curr])
        if(z.start.curr < n.year){#if you don't recruit in final year
          for(g in (z.start.curr+1):n.year){
            log.prop.back <- log.prop.back + dbinom(z.curr[g],1,model$phi[i,g-1]*z.curr[g-1],log=TRUE)
          }
        }
        lp.initial.total <- lp.initial.entry + lp.initial.y + lp.initial.surv + log.prior.curr
        lp.proposed.total <- lp.proposed.entry + lp.proposed.y + lp.proposed.surv + log.prior.prop
        
        #MH step
        log_MH_ratio <- (lp.proposed.total + log.prop.back) - (lp.initial.total + log.prop.for)
        # log_MH_ratio <- (lp.proposed) - (lp.initial)
        accept <- decide(log_MH_ratio)
        
        if(accept){
          mvSaved["z.start",1][i] <<- model[["z.start"]][i]
          mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["N.recruit",1] <<- model[["N.recruit"]]
          mvSaved["ER",1] <<- model[["ER"]]
          for(g in 1:n.year){
            for(j in 1:J[g]){
              mvSaved["d2",1][i,g,j] <<- model[["d2"]][i,g,j]
              mvSaved["pd.B",1][i,g,j] <<- model[["pd.B"]][i,g,j]
              mvSaved["pd.L",1][i,g,j] <<- model[["pd.L"]][i,g,j]
              mvSaved["pd.R",1][i,g,j] <<- model[["pd.R"]][i,g,j]
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
          for(g in 1:n.year){
            for(j in 1:J[g]){
              model[["d2"]][i,g,j] <<- mvSaved["d2",1][i,g,j]
              model[["pd.B"]][i,g,j] <<- mvSaved["pd.B",1][i,g,j]
              model[["pd.L"]][i,g,j] <<- mvSaved["pd.L",1][i,g,j]
              model[["pd.R"]][i,g,j] <<- mvSaved["pd.R",1][i,g,j]
            }
          }
          #set these logProbs back
          model$calculate(N.recruit.nodes)
          model$calculate(N.nodes[1])
          model$calculate(y.B.nodes[i.idx])
          model$calculate(y.L.nodes[i.idx])
          model$calculate(y.R.nodes[i.idx])
          model$calculate(z.nodes[i])
        }
      }
    }
    #3) update z.super: Metropolis-Hastings
    #entry counts current coming out of undetected ind update
    for(up in 1:z.super.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected individual
      if(updown==0){#subtract
      #find all z's currently on
      z.on <- which(model$z.super==1)
      non.init <- length(z.on)
      pick <- rcat(1,rep(1/non.init,non.init))
      pick <- z.on[pick]
      if(z.obs[pick]==1){ #is this individual detected?
        reject <- TRUE #if so, we reject (could never select these inds, but then need to account for asymmetric proposal)
      }
      if(!reject){
        z.start.curr <- model$z.start[pick]
        z.curr <- model$z[pick,]
        
        #p select off guy
        log.p.select.for <- log(1/non.init)
        #log multinomial coefficient prior
        log.z.prior.for <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr+1)))
        pick.idx <- seq(pick,M*n.year,M) #used to reference correct y nodes
        
        #get initial logProbs (survival logProb does not change)
        lp.initial.N <- model$getLogProb(N.nodes[1])
        lp.initial.N.recruit <- model$getLogProb(N.recruit.nodes)
        lp.initial.y <- model$getLogProb(y.B.nodes[pick.idx]) + 
          model$getLogProb(y.L.nodes[pick.idx]) + 
          model$getLogProb(y.R.nodes[pick.idx])
        lp.initial.surv <- model$getLogProb(z.nodes[pick])
        
        # propose new N.super/z.super/z.start/z.stop
        model$N.super <<-  model$N.super - 1
        model$z.super[pick] <<- 0
        model$z.start[pick] <<- 0
        model$z.stop[pick] <<- 0
        model$z[pick,] <<- rep(0,n.year)
        
        #update N, N.recruit, N.survive
        #1) Update N
        model$N <<- model$N - z.curr
        #2) Update N.recruit
        if(z.start.curr > 1){ #if wasn't in pop in year 1
          model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
        }
        #3) Update N.survive
        model$N.survive <<- model$N[2:n.year]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
        model$calculate(ER.nodes) #update ER when N updated
        model$calculate(d2.nodes[pick.idx]) #update d2 nodes when a z changes
        model$calculate(pd.B.nodes[pick.idx]) #update pd nodes when z changes
        model$calculate(pd.L.nodes[pick.idx]) #update pd nodes when z changes
        model$calculate(pd.R.nodes[pick.idx]) #update pd nodes when z changes
        
        #Reverse proposal probs
        recruit.probs.back <- c(model$lambda.y1, model$ER)
        recruit.probs.back <- recruit.probs.back / sum(recruit.probs.back)
        log.prop.back <- log(recruit.probs.back[z.start.curr])
        if(z.start.curr < n.year){
          for(g in (z.start.curr+1):n.year){
            log.prop.back <- log.prop.back + dbinom(z.curr[g],1,model$phi[pick,g-1]*z.curr[g-1],log=TRUE)
          }
        }
        
        #get proposed logProbs for N, N.recruit, and y
        lp.proposed.N <- model$calculate(N.nodes[1])
        lp.proposed.N.recruit <- model$calculate(N.recruit.nodes)
        lp.proposed.y <- model$calculate(y.B.nodes[pick.idx]) + #will always be 0
          model$calculate(y.L.nodes[pick.idx]) + 
          model$calculate(y.R.nodes[pick.idx])
        lp.proposed.surv <- model$calculate(z.nodes[pick]) #will always be 0
        
        lp.initial.total <- lp.initial.N + lp.initial.y + lp.initial.N.recruit + lp.initial.surv
        lp.proposed.total <- lp.proposed.N + lp.proposed.y + lp.proposed.N.recruit + lp.proposed.surv
        
        #backwards prior and select probs
        #move from class z.start.curr in z.super==0 to class g in z.super==1
        entry.counts.prop <- entry.counts.curr
        entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
        entry.counts.prop[n.year + 1] <- entry.counts.prop[n.year + 1] + 1
        
        #p select on guy
        noff.back <- sum(model$z.super == 0)
        log.p.select.back <- log(1/noff.back)
        #log multinomial coefficient prior
        log.z.prior.back <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop+1)))
        log.prop.for <- 0
        #MH step
        log_MH_ratio <- (lp.proposed.total + log.z.prior.back + log.p.select.back + log.prop.back) -
          (lp.initial.total + log.z.prior.for + log.p.select.for + log.prop.for)
        
        accept <- decide(log_MH_ratio)
        if(accept){
          mvSaved["z.start",1][pick] <<- model[["z.start"]][pick]
          mvSaved["z.stop",1][pick] <<- model[["z.stop"]][pick]
          mvSaved["z",1][pick,] <<- model[["z"]][pick,]
          mvSaved["z.super",1] <<- model[["z.super"]]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["N.recruit",1] <<- model[["N.recruit"]]
          mvSaved["N.super",1][1] <<- model[["N.super"]]
          mvSaved["ER",1] <<- model[["ER"]]
          for(g in 1:n.year){
            for(j in 1:J[g]){
              mvSaved["d2",1][pick,g,j] <<- model[["d2"]][pick,g,j]
              mvSaved["pd.B",1][pick,g,j] <<- model[["pd.B"]][pick,g,j]
              mvSaved["pd.L",1][pick,g,j] <<- model[["pd.L"]][pick,g,j]
              mvSaved["pd.R",1][pick,g,j] <<- model[["pd.R"]][pick,g,j]
            }
          }
          entry.counts.curr <- entry.counts.prop
        }else{
          model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
          model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
          model[["z"]][pick,] <<- mvSaved["z",1][pick,]
          model[["z.super"]] <<- mvSaved["z.super",1]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["N.recruit"]] <<- mvSaved["N.recruit",1]
          model[["N.super"]] <<- mvSaved["N.super",1][1]
          model[["ER"]] <<- mvSaved["ER",1]
          for(g in 1:n.year){
            for(j in 1:J[g]){
              model[["d2"]][pick,g,j] <<- mvSaved["d2",1][pick,g,j]
              model[["pd.B"]][pick,g,j] <<- mvSaved["pd.B",1][pick,g,j]
              model[["pd.L"]][pick,g,j] <<- mvSaved["pd.L",1][pick,g,j]
              model[["pd.R"]][pick,g,j] <<- mvSaved["pd.R",1][pick,g,j]
            }
          }
          #set these logProbs back
          model$calculate(y.B.nodes[pick.idx])
          model$calculate(y.L.nodes[pick.idx])
          model$calculate(y.R.nodes[pick.idx])
          model$calculate(z.nodes[pick])
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
        }
      }
    }else{#add
      if(model$N.super[1] < M){ #cannot update if z.super maxed out. Need to raise M
        z.off <- which(model$z.super==0)
        noff.init <- length(z.off)
        pick <- rcat(1,rep(1/noff.init,noff.init)) #select one of these individuals
        pick <- z.off[pick]
        pick.idx <- seq(pick,M*n.year,M)
        
        non.init <- sum(model$z.super == 1)
        
        #p select off guy
        log.p.select.for <- log(1/noff.init)
        
        #log multinomial coefficient prior
        log.z.prior.for <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr+1)))
        
        #get initial logProbs (survival logProb does not change)
        lp.initial.N <- model$getLogProb(N.nodes[1])
        lp.initial.N.recruit <- model$getLogProb(N.recruit.nodes)
        lp.initial.y <- model$getLogProb(y.B.nodes[pick.idx]) + #will always be 0
          model$getLogProb(y.L.nodes[pick.idx]) +
          model$getLogProb(y.R.nodes[pick.idx])
        lp.initial.surv <- model$getLogProb(z.nodes[pick]) #will always be 0
        
        # Propose new z.start for the new on individual
        recruit.probs.for <- c(model$lambda.y1, model$ER)
        recruit.probs.for <- recruit.probs.for / sum(recruit.probs.for)
        z.start.prop <- rcat(1, recruit.probs.for)  # propose entry cohort
        log.prop.for <- log(recruit.probs.for[z.start.prop])
        model$z.start[pick] <<- z.start.prop
        
        # Simulate survival path
        model$z[pick,] <<- 0 # initialize to 0
        model$z[pick, z.start.prop] <<- 1
        if(z.start.prop < n.year){
          for(g in (z.start.prop+1):n.year){
            model$z[pick, g] <<- rbinom(1, 1, model$phi[pick, g-1] * model$z[pick, g-1])
            log.prop.for <- log.prop.for + dbinom(model$z[pick, g], 1, model$phi[pick, g-1] * model$z[pick, g-1], log=TRUE)
          }
        }
        # Update z.stop
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
        if(model$z.start[pick] > 1){ #if wasn't in pop in year 1
          model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
        }
        #3) Update N.survive
        model$N.survive <<- model$N[2:n.year] - model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
        model$calculate(ER.nodes) #update ER when N updated
        model$calculate(d2.nodes[pick.idx]) #update d2 nodes when a z changes
        model$calculate(pd.B.nodes[pick.idx]) #update pd nodes when z changes
        model$calculate(pd.L.nodes[pick.idx]) #update pd nodes when z changes
        model$calculate(pd.R.nodes[pick.idx]) #update pd nodes when z changes
        #get proposed logprobs for N and y
        lp.proposed.N <- model$calculate(N.nodes[1])
        lp.proposed.N.recruit <- model$calculate(N.recruit.nodes)
        lp.proposed.y <- model$calculate(y.B.nodes[pick.idx]) +  #will always be 0
          model$calculate(y.L.nodes[pick.idx]) + 
          model$calculate(y.R.nodes[pick.idx])
        lp.proposed.surv <- model$calculate(z.nodes[pick])
        
        lp.initial.total <- lp.initial.N + lp.initial.y + lp.initial.N.recruit + lp.initial.surv
        lp.proposed.total <- lp.proposed.N + lp.proposed.y + lp.proposed.N.recruit + lp.proposed.surv
        
        #backwards prior and select probs
        #move from class g in z.super==0 to class g in z.super==1
        entry.counts.prop <- entry.counts.curr
        entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
        entry.counts.prop[n.year + 1] <- entry.counts.prop[n.year + 1] - 1
        
        #p select on guy
        non.back <- sum(model$z.super == 1)
        log.p.select.back <- log(1/non.back)
        #log multinomial coefficient prior
        log.z.prior.back <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop+1)))
        log.prop.back <- 0
        
        #MH step
        log_MH_ratio <- (lp.proposed.total + log.z.prior.back + log.p.select.back + log.prop.back) -
          (lp.initial.total + log.z.prior.for + log.p.select.for + log.prop.for)
        
        accept <- decide(log_MH_ratio)
        if(accept){
          mvSaved["z.start",1][pick] <<- model[["z.start"]][pick]
          mvSaved["z.stop",1][pick] <<- model[["z.stop"]][pick]
          mvSaved["z",1][pick,] <<- model[["z"]][pick,]
          mvSaved["z.super",1] <<- model[["z.super"]]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["N.recruit",1] <<- model[["N.recruit"]]
          mvSaved["N.super",1][1] <<- model[["N.super"]]
          mvSaved["ER",1] <<- model[["ER"]]
          for(g in 1:n.year){
            for(j in 1:J[g]){
              mvSaved["d2",1][pick,g,j] <<- model[["d2"]][pick,g,j]
              mvSaved["pd.B",1][pick,g,j] <<- model[["pd.B"]][pick,g,j]
              mvSaved["pd.L",1][pick,g,j] <<- model[["pd.L"]][pick,g,j]
              mvSaved["pd.R",1][pick,g,j] <<- model[["pd.R"]][pick,g,j]
            }
          }
          entry.counts.curr <- entry.counts.prop
        }else{
          model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
          model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
          model[["z"]][pick,] <<- mvSaved["z",1][pick,]
          model[["z.super"]] <<- mvSaved["z.super",1]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["N.recruit"]] <<- mvSaved["N.recruit",1]
          model[["N.super"]] <<- mvSaved["N.super",1][1]
          model[["ER"]] <<- mvSaved["ER",1]
          for(g in 1:n.year){
            for(j in 1:J[g]){
              model[["d2"]][pick,g,j] <<- mvSaved["d2",1][pick,g,j]
              model[["pd.B"]][pick,g,j] <<- mvSaved["pd.B",1][pick,g,j]
              model[["pd.L"]][pick,g,j] <<- mvSaved["pd.L",1][pick,g,j]
              model[["pd.R"]][pick,g,j] <<- mvSaved["pd.R",1][pick,g,j]
            }
          }
          #set these logProbs back
          model$calculate(y.B.nodes[pick.idx])
          model$calculate(y.L.nodes[pick.idx])
          model$calculate(y.R.nodes[pick.idx])
          model$calculate(z.nodes[pick])
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
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
    n.year <- control$n.year
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
          for(g in 1:n.year){ #"can't do math in >2 dimensions" loop
            for(j in 1:J[g]){
              y.L.true.cand[this.i,g,j] <- y.L.true[cand.i,g,j]
              y.L.true.cand[cand.i,g,j] <- y.L.true[this.i,g,j]
            }
          }
          #update ll.y.L
          for(g in 1:n.year){
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
          for(g in 1:n.year){
            lp_initial <- lp_initial + ll.y.L[this.i,g,1] + ll.y.L[cand.i,g,1]
            lp_proposed <- lp_proposed + ll.y.L.cand[this.i,g,1] + ll.y.L.cand[cand.i,g,1]
          }
          log_MH_ratio <- (lp_proposed+log(backprobs[this.i])) - (lp_initial+log(propprobs[cand.i]))
          accept <- decide(log_MH_ratio)
          if(accept){
            for(g in 1:n.year){ #"can't do math in >2 dimensions" loop
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
    n.year <- control$n.year
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
          for(g in 1:n.year){ #"can't do math in >2 dimensions" loop
            for(j in 1:J[g]){
              y.R.true.cand[this.i,g,j] <- y.R.true[cand.i,g,j]
              y.R.true.cand[cand.i,g,j] <- y.R.true[this.i,g,j]
            }
          }
          #update ll.y.R
          for(g in 1:n.year){
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
          for(g in 1:n.year){
            lp_initial <- lp_initial + ll.y.R[this.i,g,1] + ll.y.R[cand.i,g,1]
            lp_proposed <- lp_proposed + ll.y.R.cand[this.i,g,1] + ll.y.R.cand[cand.i,g,1]
          }
          log_MH_ratio <- (lp_proposed+log(backprobs[this.i])) - (lp_initial+log(propprobs[cand.i]))
          accept <- decide(log_MH_ratio)
          
          if(accept){
            for(g in 1:n.year){ #"can't do math in >2 dimensions" loop
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

sSampler <- nimbleFunction(
  # name = 'sampler_RW',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    i <- control$i
    J <- control$J
    n.primary <- control$n.primary
    xlim <- control$xlim
    ylim <- control$ylim
    ## control list extraction
    # logScale            <- extractControlElement(control, 'log',                 FALSE)
    # reflective          <- extractControlElement(control, 'reflective',          FALSE)
    adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
    adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
    adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
    scale               <- extractControlElement(control, 'scale',               1)
    
    ## node list generation
    #s prior nodes
    s.nodes <- model$expandNodeNames(paste0("s[",i,",1:2]"))
    
    #occasion-specific detection nodes
    d2.nodes <- pd.B.nodes <- pd.L.nodes <- pd.R.nodes <- character(0)
    y.B.nodes <- y.L.nodes <- y.R.nodes <- character(0)
    for(g.setup in 1:n.primary){
      d2.nodes <- c(d2.nodes,model$expandNodeNames(paste0("d2[",i,",",g.setup,",1:",J[g.setup],"]")))
      pd.B.nodes <- c(pd.B.nodes,model$expandNodeNames(paste0("pd.B[",i,",",g.setup,",1:",J[g.setup],"]")))
      pd.L.nodes <- c(pd.L.nodes,model$expandNodeNames(paste0("pd.L[",i,",",g.setup,",1:",J[g.setup],"]")))
      pd.R.nodes <- c(pd.R.nodes,model$expandNodeNames(paste0("pd.R[",i,",",g.setup,",1:",J[g.setup],"]")))
      y.B.nodes <- c(y.B.nodes,model$expandNodeNames(paste0("y.B.true[",i,",",g.setup,",1:",J[g.setup],"]")))
      y.L.nodes <- c(y.L.nodes,model$expandNodeNames(paste0("y.L.true[",i,",",g.setup,",1:",J[g.setup],"]")))
      y.R.nodes <- c(y.R.nodes,model$expandNodeNames(paste0("y.R.true[",i,",",g.setup,",1:",J[g.setup],"]")))
    }
    calcNodes <- c(s.nodes,d2.nodes,pd.B.nodes,pd.L.nodes,pd.R.nodes,y.B.nodes,y.L.nodes,y.R.nodes)
    
    ## numeric value generation
    scaleOriginal <- scale
    timesRan      <- 0
    timesAccepted <- 0
    timesAdapted  <- 0
    scaleHistory  <- c(0, 0)   ## scaleHistory
    acceptanceHistory  <- c(0, 0)   ## scaleHistory
    if(nimbleOptions('MCMCsaveHistory')) {
      saveMCMChistory <- TRUE
    } else saveMCMChistory <- FALSE
    optimalAR     <- 0.44
    gamma1        <- 0
    ## checks
    # if(length(targetAsScalar) > 1)   stop('cannot use RW sampler on more than one target; try RW_block sampler')
    # if(model$isDiscrete(target))     stop('cannot use RW sampler on discrete-valued target; try slice sampler')
    # if(logScale & reflective)        stop('cannot use reflective RW sampler on a log scale (i.e. with options log=TRUE and reflective=TRUE')
    if(adaptFactorExponent < 0)      stop('cannot use RW sampler with adaptFactorExponent control parameter less than 0')
    if(scale < 0)                    stop('cannot use RW sampler with scale control parameter less than 0')
  },
  run = function() {
    z.super <- model$z.super[i]
    if(z.super==0){#propose from uniform prior
      model$s[i,1:2] <<- c(runif(1,xlim[1],xlim[2]),runif(1,ylim[1],ylim[2]))
      model$calculate(s.nodes)
      copy(from = model, to = mvSaved, row = 1, nodes = s.nodes, logProb = TRUE)
    }else{#MH
      s.cand <- c(rnorm(1,model$s[i,1],scale),rnorm(1,model$s[i,2],scale))
      inbox <- s.cand[1]<xlim[2]&s.cand[1]>xlim[1]&s.cand[2]<ylim[2]&s.cand[2]>ylim[1]
      if(inbox){
        #initial log probability: s prior plus observation likelihoods while alive
        lp.initial <- model$getLogProb(s.nodes)
        for(g.use in 1:n.primary){
          if(model$z[i,g.use]==1){
            lp.initial <- lp.initial+model$getLogProb(y.B.nodes[g.use])+model$getLogProb(y.L.nodes[g.use])+model$getLogProb(y.R.nodes[g.use])
          }
        }
        #propose activity center
        model$s[i,1:2] <<- s.cand
        lp.proposed <- model$calculate(s.nodes)
        
        #only detection terms in primary occasions where the individual is alive can change
        for(g.use in 1:n.primary){
          if(model$z[i,g.use]==1){
            model$calculate(d2.nodes[g.use])
            model$calculate(pd.B.nodes[g.use])
            model$calculate(pd.L.nodes[g.use])
            model$calculate(pd.R.nodes[g.use])
            lp.proposed <- lp.proposed+model$calculate(y.B.nodes[g.use])+model$calculate(y.L.nodes[g.use])+model$calculate(y.R.nodes[g.use])
          }
        }
        
        log_MH_ratio <- lp.proposed-lp.initial
        accept <- decide(log_MH_ratio)
        if(accept) {
          copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        } else {
          copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        }
        if(adaptive){ #we only tune for z=1 proposals
          adaptiveProcedure(accept)
        }
      }
    }
  },
  methods = list(
    adaptiveProcedure = function(jump = logical()) {
      timesRan <<- timesRan + 1
      if(jump)     timesAccepted <<- timesAccepted + 1
      if(timesRan %% adaptInterval == 0) {
        acceptanceRate <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        if(saveMCMChistory) {
          setSize(scaleHistory, timesAdapted)                 ## scaleHistory
          scaleHistory[timesAdapted] <<- scale                ## scaleHistory
          setSize(acceptanceHistory, timesAdapted)            ## scaleHistory
          acceptanceHistory[timesAdapted] <<- acceptanceRate  ## scaleHistory
        }
        gamma1 <<- 1/((timesAdapted + 3)^adaptFactorExponent)
        gamma2 <- 10 * gamma1
        adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
        scale <<- scale * adaptFactor
        ## If there are upper and lower bounds, enforce a maximum scale of
        ## 0.5 * (upper-lower).  This is arbitrary but reasonable.
        ## Otherwise, for a poorly-informed posterior,
        ## the scale could grow without bound to try to reduce
        ## acceptance probability.  This creates enormous cost of
        ## reflections.
        # if(reflective) {
        #   lower <- model$getBound(target, 'lower')
        #   upper <- model$getBound(target, 'upper')
        #   if(scale >= 0.5*(upper-lower)) {
        #     scale <<- 0.5*(upper-lower)
        #   }
        # }
        timesRan <<- 0
        timesAccepted <<- 0
      }
    },
    getScaleHistory = function() {       ## scaleHistory
      returnType(double(1))
      if(saveMCMChistory) {
        return(scaleHistory)
      } else {
        print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
        return(numeric(1, 0))
      }
    },          
    getAcceptanceHistory = function() {  ## scaleHistory
      returnType(double(1))
      if(saveMCMChistory) {
        return(acceptanceHistory)
      } else {
        print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
        return(numeric(1, 0))
      }
    },
    ##getScaleHistoryExpanded = function() {                                                 ## scaleHistory
    ##    scaleHistoryExpanded <- numeric(timesAdapted*adaptInterval, init=FALSE)            ## scaleHistory
    ##    for(iTA in 1:timesAdapted)                                                         ## scaleHistory
    ##        for(j in 1:adaptInterval)                                                      ## scaleHistory
    ##            scaleHistoryExpanded[(iTA-1)*adaptInterval+j] <- scaleHistory[iTA]         ## scaleHistory
    ##    returnType(double(1)); return(scaleHistoryExpanded) },                             ## scaleHistory
    reset = function() {
      scale <<- scaleOriginal
      timesRan      <<- 0
      timesAccepted <<- 0
      timesAdapted  <<- 0
      if(saveMCMChistory) {
        scaleHistory  <<- c(0, 0)    ## scaleHistory
        acceptanceHistory  <<- c(0, 0)
      }
      gamma1 <<- 0
    }
  )
)
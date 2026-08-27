sSampler <- nimbleFunction(
  # name = 'sampler_RW',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    i <- control$i
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
    # targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes <- model$getDependencies(target)
    
    #old approach: calculate every dependency of s[i,] for every proposal
    #calcNodes <- model$getDependencies(target)
    
    #group s[i,] dependencies by the primary occasion-specific z[i,g] they also depend on.
    #Nodes associated with exactly one z[i,g] can be skipped when z[i,g]=0.
    #Nodes associated with no z[i,g], or with more than one z[i,g], are always recalculated.
    n.primary <- length(model$z[i,])
    dep.count <- integer(length(calcNodes))
    z.dep.nodes <- vector("list",n.primary)
    for(g in 1:n.primary){
      z.dep <- model$getDependencies(paste0("z[",i,",",g,"]"))
      z.dep.nodes[[g]] <- calcNodes[calcNodes %in% z.dep]
      if(length(z.dep.nodes[[g]])>0){
        idx <- match(z.dep.nodes[[g]],calcNodes)
        dep.count[idx] <- dep.count[idx]+1L
      }
    }
    
    #flatten occasion-local node lists for use in compiled run code
    local.nodes <- character(0)
    local.node.start <- integer(n.primary)
    local.node.n <- integer(n.primary)
    for(g in 1:n.primary){
      these.nodes <- z.dep.nodes[[g]]
      if(length(these.nodes)>0){
        these.idx <- match(these.nodes,calcNodes)
        these.nodes <- these.nodes[dep.count[these.idx]==1]
      }
      local.node.start[g] <- length(local.nodes)+1L
      local.node.n[g] <- length(these.nodes)
      if(length(these.nodes)>0){
        local.nodes <- c(local.nodes,these.nodes)
      }
    }
    other.nodes <- calcNodes[dep.count!=1]
    
    # calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    # isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)   ## should be made faster
    # calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
    # calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
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
      model$s[i, 1:2] <<- c(runif(1, xlim[1], xlim[2]), runif(1, ylim[1], ylim[2]))
      
      #old full-dependency calculation
      #model$calculate(calcNodes)
      model$calculate(other.nodes)
      
      #old full-dependency copy
      #copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
      copy(from = model, to = mvSaved, row = 1, nodes = other.nodes, logProb = TRUE)
    }else{#MH
      s.cand <- c(rnorm(1,model$s[i,1],scale), rnorm(1,model$s[i,2],scale))
      inbox <- s.cand[1]< xlim[2] & s.cand[1]> xlim[1] & s.cand[2] < ylim[2] & s.cand[2] > ylim[1]
      if(inbox){
        
        #old full-dependency initial log probability
        #model_lp_initial <- model$getLogProb(calcNodes)
        model_lp_initial <- model$getLogProb(other.nodes)
        for(g in 1:n.primary){
          if(model$z[i,g]==1&local.node.n[g]>0){
            node.start <- local.node.start[g]
            node.end <- node.start+local.node.n[g]-1
            model_lp_initial <- model_lp_initial+model$getLogProb(local.nodes[node.start:node.end])
          }
        }
        
        model$s[i, 1:2] <<- s.cand
        
        #old full-dependency proposed calculation
        #model_lp_proposed <- model$calculate(calcNodes)
        model_lp_proposed <- model$calculate(other.nodes)
        for(g in 1:n.primary){
          if(model$z[i,g]==1&local.node.n[g]>0){
            node.start <- local.node.start[g]
            node.end <- node.start+local.node.n[g]-1
            model_lp_proposed <- model_lp_proposed+model$calculate(local.nodes[node.start:node.end])
          }
        }
        
        log_MH_ratio <- model_lp_proposed - model_lp_initial
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

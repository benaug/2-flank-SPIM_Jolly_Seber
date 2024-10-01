#This version uses a 2D trap operation matrix and camera number matrix (year by trap)
#This does not allow the camera number to change through time at the same
#station, but runs faster than 3D version that allows one station to have
#some occasions with 1 camera and others with 2

#this version also considers activity center relocation between primary periods.
#the current algorithm works, but the flank swap updates could potentially be made
#more efficient, with lower rejection rates
library(nimble)
library(coda)
source("sim.2flank.JS.2D.R")
source("init.2flank.JS.mobileAC.2D.R")
source("LRmatch.R")
source("Nimble Model 2-flank-SPIM-JS-mobileAC-2D.R")
source("Nimble Functions 2-flank-SPIM-JS-mobileAC-2D.R") #contains custom distributions and updates
source("sSampler Multi.R") # activity center sampler that proposes from prior when z.super=0.
#this one works for fixed activity centers over years only

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

n.year <- 5 #number of years
lambda.y1 <- 50 #expected N in year 1
gamma <- rep(0.2,n.year-1) #yearly per-capita recruitment
phi <- rep(0.85,n.year-1) #yearly survival, model file set up for fixed
#yearly detection probabilities at activity center. Model file set up for p0.L=p0.R
p0.B <- rep(0.05,n.year) #both-flank detections
p0.L <- rep(0.05,n.year) #left-flank detections
p0.R <- rep(0.05,n.year) #right-flank detections
sigma <- rep(0.5,n.year) #yearly detection function spatial scale
sigma.move <- 1 #movement sigma, fixed over primary periods
K <- rep(10,n.year) #yearly sampling occasions

buff <- 2 #state space buffer. Buffers maximal x and y dimensions of X below across years
X <- vector("list",n.year) #one trapping array per year
for(g in 1:n.year){ #using same trapping array every year here
  X[[g]] <- as.matrix(expand.grid(3:11,3:11))
}
J <- unlist(lapply(X,nrow))
J.max <- max(J)
K.max <- max(K)

#year by trap matrix indicating which traps had 1 vs. 2 cameras in each year
#if you do not include 2 camera stations, you cannot observe both-flank captures
J.cams <- matrix(0,n.year,J.max)
for(g in 1:n.year){
  J.cams[g,1:J[g]] <- 1 #start with all 1 cam stations
  # J.cams[g,seq(1,J[g],2)] <- 2 #set some to 2 (or not)
}
#year by trap matrix of camera operation - assuming all operational for all occasions in each year
K2D <- matrix(0,n.year,J.max)
for(g in 1:n.year){
  K2D[g,1:J[g]] <- K[g]
}

#You may know some flank linkages without both-side captures. If NA, we assume we do not.
# n.fixed <- 5 #minimum number of matched flanks for simulation
n.fixed <- NA #supply NA for only n.B flank matches to be known (if there are any).
#must supply this for simulation so we know to retain left and right known flank individuals with no captures
data <- sim.2flank.JS.2D(lambda.y1=lambda.y1,gamma=gamma,n.year=n.year,
            phi=phi,p0.B=p0.B,p0.L=p0.L,p0.R=p0.R,sigma=sigma,X=X,buff=buff,
            K=K,K2D=K2D,J.cams=J.cams,n.fixed=n.fixed,sigma.move=sigma.move)
#true N.super
data$truth$N[1] + sum(data$truth$N.recruit)

#captures per individual by year for each capture types, are these realistic?
apply(data$y.B,c(1,2),sum)
apply(data$y.L,c(1,2),sum)
apply(data$y.R,c(1,2),sum)

##Initialize##
M <- 200 #data augmentation level. Check N.super posterior to make sure it never hits M

#initialize data. Can init at truth for simulated data sets, useful for development
nimbuild <- init.2flank.JS.mobileAC.2D(data=data,M=M,n.fixed=data$n.fixed,initTrue=FALSE)

#constants for Nimble
constants <- list(n.year=data$n.year,M=M,J=data$J,xlim=data$xlim,ylim=data$ylim,
                  K2D=data$K2D,J.cams=data$J.cams,
                  n.L=nimbuild$n.L,n.R=nimbuild$n.R)
#inits for Nimble. includes left and right flank histories that are partially or fully latent
Niminits <- list(N=nimbuild$N,lambda.y1=nimbuild$N[1],
                 N.survive=nimbuild$N.survive,N.recruit=nimbuild$N.recruit,
                 ER=nimbuild$N.recruit,N.super=nimbuild$N.super,z.super=nimbuild$z.super,
                 z=nimbuild$z,z.start=nimbuild$z.start,z.stop=nimbuild$z.stop,
                 s=nimbuild$s,ID.L=nimbuild$ID.L,ID.R=nimbuild$ID.R,
                 y.L.true=nimbuild$y.true[,,,2],
                 y.R.true=nimbuild$y.true[,,,3],
                 sigma=2, #set sigma higher than expected so that starting logProb is finite
                 p0.B=rep(0.1,data$n.year),p0.S=rep(0.1,data$n.year))

#data for Nimble
Nimdata <- list(y.B.true=nimbuild$y.true[,,,1],X=nimbuild$X.nim) #both flank data is fully observed

# set parameters to monitor
parameters <- c('N','gamma.fixed','N.recruit','N.survive','N.super',
                'lambda.y1','phi.fixed','p0.B','p0.S','sigma','sigma.move')
parameters2 <- c("ID.L","ID.R") #monitor these with reduced thinning rate

nt <- 1 #thinning rate 1
nt2 <- 10 #thinning rate 2
# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)

#OK! what are we doing here? If you just let nimble configure as normal, it will assign incorrect samplers
#to z and N objects. We could then remove them and replace them, but it is much faster to not let nimble
#make the assignments in the first place. So! put all terms with priors in config.nodes here except for
#N.recruit. If you change the model parameters, you will need to make the same changes here. Finally, 
#we have to tell nimble which nodes to assign samplers for for the individual covariate when manually
#instructing nimble which samplers to assign.
config.nodes <- c('phi.fixed','gamma.fixed','lambda.y1','p0.B','p0.S','sigma','sigma.move')
# config.nodes <- c()
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      monitors2=parameters2,thin2=nt2,
                      nodes=config.nodes,useConjugacy = FALSE)

###*required* sampler replacements###
z.super.ups <- round(M*0.25) #how many z.super update proposals per iteration?
#20% of M seems reasonable, but optimal will depend on data set
#loop here bc potentially different numbers of traps to vectorize in each year
y.B.nodes <- pd.B.nodes <- c()
y.L.nodes <- pd.L.nodes <- c()
y.R.nodes <- pd.R.nodes <- c()
d2.nodes <- c()
for(g in 1:data$n.year){
  #if you change y structure, change here
  d2.nodes <- c(d2.nodes,Rmodel$expandNodeNames(paste0("d2[1:",M,",",g,",1:",data$J[g],"]")))
  y.B.nodes <- c(y.B.nodes,Rmodel$expandNodeNames(paste0("y.B.true[1:",M,",",g,",1:",data$J[g],"]")))
  pd.B.nodes <- c(pd.B.nodes,Rmodel$expandNodeNames(paste0("pd.B[1:",M,",",g,",1:",data$J[g],"]")))
  y.L.nodes <- c(y.L.nodes,Rmodel$expandNodeNames(paste0("y.L.true[1:",M,",",g,",1:",data$J[g],"]")))
  pd.L.nodes <- c(pd.L.nodes,Rmodel$expandNodeNames(paste0("pd.L[1:",M,",",g,",1:",data$J[g],"]")))
  y.R.nodes <- c(y.R.nodes,Rmodel$expandNodeNames(paste0("y.R.true[1:",M,",",g,",1:",data$J[g],"]")))
  pd.R.nodes <- c(pd.R.nodes,Rmodel$expandNodeNames(paste0("pd.R[1:",M,",",g,",1:",data$J[g],"]")))
}
N.nodes <- Rmodel$expandNodeNames(paste0("N"))
N.survive.nodes <- Rmodel$expandNodeNames(paste0("N.survive[1:",data$n.year-1,"]"))
N.recruit.nodes <- Rmodel$expandNodeNames(paste0("N.recruit[1:",data$n.year-1,"]"))
ER.nodes <- Rmodel$expandNodeNames(paste0("ER[1:",data$n.year-1,"]"))
z.nodes <- Rmodel$expandNodeNames(paste0("z[1:",M,",1]"))
calcNodes <- c(N.nodes,N.recruit.nodes,y.B.nodes,y.L.nodes,y.R.nodes,z.nodes) #the ones that need likelihoods updated in mvSaved
conf$addSampler(target = c("z"),
                type = 'zSampler',control = list(M=M,n.year=data$n.year,J=data$J,
                                                 z.super.ups=z.super.ups,
                                                 y.B.nodes=y.B.nodes,pd.B.nodes=pd.B.nodes,
                                                 y.L.nodes=y.L.nodes,pd.L.nodes=pd.L.nodes,
                                                 y.R.nodes=y.R.nodes,pd.R.nodes=pd.R.nodes,
                                                 N.nodes=N.nodes,
                                                 z.nodes=z.nodes,ER.nodes=ER.nodes,
                                                 N.survive.nodes=N.survive.nodes,
                                                 N.recruit.nodes=N.recruit.nodes,
                                                 calcNodes=calcNodes), silent = TRUE)

#left and right flank updates
J.max <- max(data$J)
conf$addSampler(target = paste0("y.L.true[1:",M,",1:",data$n.year,",1:",J.max,"]"),
                type = 'IDLSampler',control = list(K2D=data$K2D,J.cams=data$J.cams,n.fixed=data$n.fixed,
                                                   n.year=data$n.year,M=nimbuild$M,J=data$J,K=data$K,
                                                   n.L=nimbuild$n.L,prop.scale=1),silent = TRUE)
conf$addSampler(target = paste0("y.R.true[1:",M,",1:",data$n.year,",1:",J.max,"]"),
                type = 'IDRSampler',control = list(K2D=data$K2D,J.cams=data$J.cams,n.fixed=data$n.fixed,
                                                   n.year=data$n.year,M=nimbuild$M,J=data$J,K=data$K,
                                                   n.R=nimbuild$n.R,prop.scale=1),silent = TRUE)

#activity center sampler. There are 3 samplers here for these cases
#1) z.super=1 and z=1, sSampler1 uses Metropolis-Hastings
#2) z.super=1 and z=0, sSampler2 uses Metropolis-Hastings with proposal sd tuned separately from above
#3) z.super=0, sSampler3 simulates entire new trajectory from priors.
for(i in 1:M){
  for(g in 1:data$n.year){
    conf$addSampler(target = paste0("s[",i,",",g,",1:2]"),
                    type = 'sSampler1',control=list(i=i,g=g,xlim=data$xlim,ylim=data$ylim,scale=1),silent = TRUE)
    conf$addSampler(target = paste0("s[",i,",",g,",1:2]"),
                    type = 'sSampler2',control=list(i=i,g=g,xlim=data$xlim,ylim=data$ylim,scale=1),silent = TRUE)
    #scale parameter here is just the starting scale. It will be tuned.
  }
  conf$addSampler(target = paste0("s[",i,",1:",data$n.year,",1:2]"),
                  type = 'sSampler3',control=list(i=i,xlim=data$xlim,ylim=data$ylim,
                                                  n.year=data$n.year,scale=1),silent = TRUE)
}

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=2) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(2000,reset=FALSE) #can extend run by rerunning this line
end.time <- Sys.time()
time1 <- end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2 <- end.time-start.time2 # post-compilation run time

mvSamples <-  as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[-c(1:250),]))

#reminder what the targets are
data$N
data$N.recruit
data$N.survive
data$N[1] + sum(data$N.recruit) #N.super

#posterior sample match probs. not removing burnin here
n.L <- nimbuild$n.L
n.R <- nimbuild$n.R
mvSamples2 <- as.matrix(Cmcmc$mvSamples2)
postL <- mvSamples2[,1:n.L] #which 1:M individual is left flank matched to on each iteration>
postR <- mvSamples2[,(n.L+1):(n.R+n.L)] #same for right flanks

#Calculate posterior prob that left flank l matches right flank r
postprobs <- matrix(NA,nrow=n.L,ncol=n.R)
for(l in 1:n.L){
  for(r in 1:n.R){
    tmp <- postL[,l]==postR[,r]
    postprobs[l,r] <- mean(tmp[-1])
  }
}

#Posterior match probability that left flanks are matched to correct right flanks (for simulated data)
#Higher match probs indicate more certainty in matches and less precision lost relative to true data
#**match probs should be 1 for 1:n.fixed individuals. If not, something went wrong.**
for(i in 1:n.L){
  idx <- which(data$ID.R==data$ID.L[i])
  if(length(idx)>0){
    print(postprobs[i,idx])
  }else{
    print(NA) #was no matching right flank
  }
}
#Posterior match probability that right flanks are matched to correct left flanks (for simulated data)
for(i in 1:n.R){
  idx <- which(data$ID.L==data$ID.R[i])
  if(length(idx)>0){
    print(postprobs[idx,i])
  }else{
    print(NA) #was no matching left flank
  }
}

#posterior for number of matching flank pairs and number of unique individuals captured
n.match <- rep(NA,nrow(postL)-1) #number of matching flank pairs
n.cap <-  rep(NA,nrow(postL)-1) #number of unique individuals captured
for(i in 2:nrow(postL)){
  n.match[i] <- sum(postL[i,]%in%postR[i,])
  n.cap[i] <- length(unique(c(postL[i,],postR[i,])))
}
tmp <- mcmc(n.match[-1])
varnames(tmp) <- "LRmatches"
plot(tmp)
sum(data$ID.L%in%data$ID.R) #truth
tmp2 <- mcmc(n.cap[-1])
varnames(tmp2) <- "n.cap"
plot(tmp2)
length(unique(c(data$ID.L,data$ID.R)))

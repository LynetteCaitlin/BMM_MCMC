
# Metropolis-Hastings algorithm for inferring
#   the overall scaled mutation rate, the mutation bias, and the strength of directional and quadratic selection
#   from the joint site frequency spectrum of two populations, as well as the split time between the two
#########################################################################################################

# load required packages
library(ggmcmc)
library(gridExtra)
library(matrixStats)

# calculate the approximate transition rate matrix of the boundary-mutation Moran model
#   with biased mutations, linear, and quadratic selection
#     (see C.Vogl and L.C.Mikula, Theoretical Population Biology 139 (2021), p.1-17 -  Eqs.(27-32))
#   - note the transition rate matrix is assumed to be quadratic (dimension NxN here)

TransitionMatrix_BG<-function(N,beta,theta,B1=0,B2=0){
  if (!is.numeric(c(N,beta,theta,B1,B2))) {
    stop("TransitionMatrix_BG(): input variables must be numeric", call. = FALSE)
  }
  if (length(c(N,beta,theta,B1,B2))!=5) {
    stop("TransitionMatrix_BG(): incorrect number of input varibles", call. = FALSE)
  }
  if (N<3) {
    stop("TransitionMatrix_BG(): N too small", call. = FALSE)
  }
  if (any(c(beta,theta)<0)) {
    stop("MetrHast(): beta,theta must be positive", call. = FALSE)
  }  
  matT=matrix(rep(0,(N+1)*(N+1)),nrow=N+1)
  vi=1:(N-1)
  C0=exp((B1+B2)/(2*N))/(1-beta*theta*sum(1/vi*exp(B1*vi/N+B2*vi*(N-vi)/(N*N))))
  matT[1,2]=theta*beta/N*C0
  matT[1,1]=1-matT[1,2]
  for(x in 1:(N-1)){
    temp=x*(N-x)/(N*N)
    matT[x+1,x]=exp(-B1/(2*N)+B2*(2*x-N)/(2*N*N))*temp
    matT[x+1,x+2]=exp(B1/(2*N)-B2*(2*x-N)/(2*N*N))*temp
    matT[x+1,x+1]=1-matT[x+1,x]-matT[x+1,x+2]
  }
  C1=exp((-B1+B2)/(2*N))/(1-(1-beta)*theta*sum(1/(N-vi)*exp(-B1*(N-vi)/N+B2*vi*(N-vi)/(N*N))))  
  matT[N+1,N]=(1-beta)*theta/N*C1
  matT[N+1,N+1]=1-matT[N+1,N]
  matT
}


# helper function: 
#   hypergeometric sampling of J haploid individuals from a population of N haploid individuals...
getH_=function(N,J){
  if (J>N) {
    stop("getH_(): J (or the dimensions of jSFS) must not be larger than N", call. = FALSE)
  }  
  matH=matrix(0,nrow=N+1,ncol=J+1)
  for(x in 0:N){
    matH[x+1,]=dhyper(0:J,x,N-x,J)
  }
  return(matH)
} 


# calculate the log-likelihood of an observed joint site-frequency spectrum 
#     given the underlying model and current parameter values
#   - note: works via spectral decomposition
LogLikelihood_=function(matP,matHJ,matHL,t,jSFS){
  eigenSystem=eigen(t(matP))
  eigenValues=eigenSystem$values
  eigenForw=t(eigenSystem$vectors)
  pi=eigenForw[1,]/sum(eigenForw[1,])
  eigenBacJ=solve(eigenForw)
  matPt=eigenBacJ%*%diag(exp(t*log(eigenValues)))%*%eigenForw
  matML=t(matHJ)%*%t(matPt)%*%diag(pi)%*%matPt%*%matHL
  sum(log(matML)*jSFS)
}

# MAIN FUNCTION
# Metropolis-Hastings algorithm (see accompanying paper for details)
MetrHast <- function(N=NULL,beta,theta,B1=0,B2=0,tau,jSFS,maxIT= 20000,burnin=NULL,fix.background=FALSE,
                     fix.B1=TRUE,fix.B2=TRUE,tune=0.2,verbose=TRUE){
  
  if(!is.null(N)){
    if(!is.numeric(c(N,beta,theta,B1,B2,tau,maxIT))) {
      stop("MetrHast(): input variables must be numeric", call. = FALSE)
    }
    if(length(c(N,beta,theta,B1,B2,tau))!=6) {
      stop("MetrHast(): incorrect number of input varibles", call. = FALSE)
    }
    if(N<3) {
      stop("MetrHast(): N too small", call. = FALSE)
    }
  }
  if(is.null(N)){
    if(!is.numeric(c(beta,theta,B1,B2,tau,maxIT))) {
      stop("MetrHast(): input variables must be numeric", call. = FALSE)
    }
    if(length(c(beta,theta,B1,B2,tau))!=5) {
      stop("MetrHast(): incorrect number of input varibles", call. = FALSE)
    }
    N <- max(nrow(jSFS),ncol(jSFS))
  }
  
  if(any(c(beta,theta,tau)<0)) {
    stop("MetrHast(): beta,theta,t must be positive", call. = FALSE)
  }
  if(!is.matrix(jSFS)) {
    stop("MetrHast(): jSFS must be a matrix", call. = FALSE)
  }
  if(sum(jSFS)<((N^2)*100)) {
    stop("SimJSFS(): size too small", call. = FALSE)
  }
  if(!is.numeric(burnin)) {
   burnin <- maxIT*0.4
  }
  
  # set up results matrix
  matResults <- matrix(rep(0,maxIT*6),ncol=6)
  
  # calculate initial log-likelihood
  J=length(jSFS[,1])-1
  L=length(jSFS[1,])-1
  matHJ=getH_(N,J)
  matHL=getH_(N,L)
  
  t <- tau*(N^2)
  
  matP=TransitionMatrix_BG(N,beta,theta,B1,B2)
  LogLike=LogLikelihood_(matP,matHJ,matHL,t,jSFS)  
  matResults[1,] <- c(LogLike,beta,theta,B1,B2,t/(N^2))
  
  # set an upper bound for the divergence time
  maxT=2*N*N
 
  # set tuning factors for different variables
  x1=tune*sum(jSFS) 
  x2=100/x1 
  x3=100000/x1  

  # start iterations!
  #   so: in each iteration, a new value is proposed for each parameter
  #         the candidate value is accepted if the acceptance ratio (of the new log likelihood vs the old) is greater than a uniform random number 
  for(it in 2:maxIT){
    
    # print progress
    xx <- NULL
      if (it %% round(maxIT / 10) == 0) {
        message(paste0(round(floor(it / maxIT * 100), -1), "%"))
      }
    
    # update split time: propose a new time near the current one 
    #   i.e. the change in time has a mean of zero and variance of N/4 (variance of the binomial distribution)
    
    if(!fix.background){
    newt=t-N/2+rbinom(1,size=N,p=0.5)
    if(3*N<=newt & newt<=maxT){
      matP=TransitionMatrix_BG(N,beta,theta,B1,B2)
      newLogLike=LogLikelihood_(matP,matHJ,matHL,newt,jSFS)
      lr=newLogLike-LogLike
      if(lr>=log(runif(1))){
        t=newt
        LogLike=newLogLike
      } 
    }
    
    # update mutation bias: propose a new value for the bias near the current one 
    #     i.e. the new bias is drawn from a beta distribution with the mean equal to the old value,
    #           and variance determined by the first smoothing parameter
    newBeta=rbeta(1,beta*x1,(1-beta)*x1)
    matP=TransitionMatrix_BG(N,newBeta,theta,B1,B2)
    newLogLike=LogLikelihood_(matP,matHJ,matHL,t,jSFS)
    lr=newLogLike-LogLike-dbeta(newBeta,beta*x1,(1-beta)*x1,log=TRUE)+dbeta(beta,newBeta*x1,(1-newBeta)*x1,log=TRUE)
    if(lr>=log(runif(1))){
      beta=newBeta
      LogLike=newLogLike
    }
    
    if(!fix.B1){
    # mutation bias and directional selection are highly correlated
    #   therefore, a joint updating step is required
    newBeta=rbeta(1,beta*(x1*0.1),(1-beta)*(x1*0.1))
    c=beta*exp(B1)/(1-beta+beta*exp(B1))
    newB1=log(c*(1-newBeta)/((1-c)*newBeta))
    matP=TransitionMatrix_BG(N,newBeta,theta,newB1,B2)
    newLogLike=LogLikelihood_(matP,matHJ,matHL,t,jSFS)
    lr=newLogLike-LogLike-dbeta(newBeta,beta*x1,(1-beta)*x1,log=TRUE)+dbeta(beta,newBeta*x1,(1-newBeta)*x1,log=TRUE)
    if(lr>=log(runif(1))){
      beta=newBeta
      B1=newB1
      LogLike=newLogLike
      }
    }
    
    }
    
    # update directional selection strength: propose a new value near the current one 
    #     i.e. the new value is drawn from a normal distribution with the mean equal to the old value,
    #           and variance determined by the second smoothing parameter
    if(!fix.B1){
      newB1=rnorm(1,B1,x2)
      matP=TransitionMatrix_BG(N,beta,theta,newB1,B2)
      newLogLike=LogLikelihood_(matP,matHJ,matHL,t,jSFS)
      lr=newLogLike-LogLike
      if(lr>=log(runif(1))){
        B1=newB1
        LogLike=newLogLike
      }
    }
    
    # update quadratic selection strength: propose a new value near the current one 
    #     i.e. the new value is drawn from a normal distribution with the mean equal to the old value,
    #           and variance determined by the second smoothing parameter
    if(!fix.B2){
      newB2=rnorm(1,B2,x3)
      matP=TransitionMatrix_BG(N,beta,theta,B1,newB2)
      newLogLike=LogLikelihood_(matP,matHJ,matHL,t,jSFS)
      lr=newLogLike-LogLike
      if(lr>=log(runif(1))){
        B2=newB2
        LogLike=newLogLike
      }
    }
    
    # update mutation rate: propose a new rate near the current one 
    #     i.e. the new rate is drawn from a gamma distribution with the mean equal to the old value,
    #           and variance determined by the first smoothing parameter
    
    if(!fix.background){
    newTheta=rgamma(1,theta*x1,x1)
    matP=TransitionMatrix_BG(N,beta,newTheta,B1,B2)
    newLogLike=LogLikelihood_(matP,matHJ,matHL,t,jSFS)
    lr=newLogLike-LogLike-dgamma(newTheta,theta*x1,x1,log=TRUE)+dgamma(theta,newTheta*x1,x1,log=TRUE)
    if(lr>=log(runif(1))){
      theta=newTheta
      LogLike=newLogLike
    }
    }
    # store the results for each iteration
    matResults[it,]=c(LogLike,beta,theta,B1,B2,t/(N^2))
    colnames(matResults) <- c("log-likelihood","mutation bias","mutation rate","directional selection","quadratic selection","split time")
    
    if(verbose){
      if(it%%100==0){
        print(c(paste0("iter=",it),paste0("log-likelihood=",LogLike),paste0("beta=",beta),paste0("theta=",theta),paste0("B1=",B1),paste0("B2=",B2),paste0("t=",t/(N^2))))
      }
    }
    
  } 
  
  #estimates
  est <- c(mean(matResults[((burnin+1):maxIT),1]),mean(matResults[((burnin+1):maxIT),2]),mean(matResults[((burnin+1):maxIT),3]), mean(matResults[((burnin+1):maxIT),4]),mean(matResults[((burnin+1):maxIT),5]),mean(matResults[((burnin+1):maxIT),6]))
  names(est) <- c("log-likelihood","mutation bias","mutation rate","directional selection","quadratic selection","split time")
  # store output
  res <- list(samples=matResults,
              estimates=est)
  
  # re-format results
  res2 <- reshape2::melt(matResults[((burnin+1):maxIT),])
  res2$Chain <- 1
  colnames(res2) <- c("Iteration", "Parameter", "value", "Chain")
  attr(res2, "nChains") <- 1
  attr(res2, "nParameters") <- 1
  attr(res2, "nIterations") <- maxIT-burnin
  attr(res2, "nBurnin") <- 0
  attr(res2, "nThin") <- 1
  
  # below, we plot the traces of all estimated parameters to assess convergence
  #     - burnin is excluded
  #     - means (after burnin) are added in horizontal white  ines
  
  if(!fix.B2){
  B2info <- res2[res2$Parameter=="quadratic selection",]
  B2trace <- ggmcmc::ggs_traceplot(B2info)+geom_hline(yintercept=mean(as.numeric(unlist(B2info[3]))),colour="white")
  grid.arrange(B2trace,nrow=1)
  }
  
  if(!fix.B1){
    B1info <- res2[res2$Parameter=="directional selection",]
    B1trace <- ggmcmc::ggs_traceplot(B1info)+geom_hline(yintercept=mean(as.numeric(unlist(B1info[3]))),colour="white")
    grid.arrange(B1trace,nrow=1)
  }
  
  if(!fix.background){
    beinfo <- res2[res2$Parameter=="mutation bias",]
    betrace <- ggmcmc::ggs_traceplot(beinfo)+geom_hline(yintercept=mean(as.numeric(unlist(beinfo[3]))),colour="white")
    
    muinfo <- res2[res2$Parameter=="mutation rate",]
    mutrace <- ggmcmc::ggs_traceplot(muinfo)+geom_hline(yintercept=mean(as.numeric(unlist(muinfo[3]))),colour="white")  
    
    tinfo <- res2[res2$Parameter=="split time",]
    ttrace <- ggmcmc::ggs_traceplot(tinfo)+geom_hline(yintercept=mean(as.numeric(unlist(tinfo[3]))),colour="white")
    
    grid.arrange(mutrace,betrace, ttrace, ncol=1)
    }
  
  llinfo <- res2[res2$Parameter=="log-likelihood",]
  lltrace <- ggmcmc::ggs_traceplot(llinfo)+geom_hline(yintercept=mean(as.numeric(unlist(llinfo[3]))),colour="white")
  grid.arrange(lltrace,nrow=1)

  return(res)
  
} 

# hypergeometrically downsample the columns of the matrix mMarg to length newM
downSampleM <- function(mMarg,newM){
  mMargDown <- matrix(0,nrow=nrow(mMarg),ncol=newM)
  oldM <- ncol(mMarg)-1
  for(i in 1:nrow(mMarg)){
    for(j in 0:oldM){
      mMargDown[i,] <- mMargDown[i,]+table(factor(rhyper(mMarg[i,j+1],j,oldM-j,newM-1),levels = 0:(newM-1)))
    }
  }
  return(mMargDown)
}


# simulate an equilibrium joint site frequency spectrum (jSFS) of two populations of dimension LxK with a specified split time
#   this can include biased mutations, linear, and quadratic selection
# the output will be a matrix of dimension LxK, but transitions are calculated based on an NxN (quadratic) matrix
SimJSFS=function(J,L,beta,theta,B1,B2,tau,size,N=NULL){
  
  
  if (!is.numeric(c(J,L,beta,theta,B1,B2,tau,size))) {
    stop("SimJSFS(): input variables must be numeric", call. = FALSE)
  }
  if (length(c(J,L,beta,theta,B1,B2))!=6) {
    stop("SimJSFS(): incorrect number of input varibles", call. = FALSE)
  }
  if (any(c(beta,theta,tau)<0)) {
    stop("SimJSFS(): beta,theta,t must be positive", call. = FALSE)
  }

  if (size<(J*L*100)) {
    stop("SimJSFS(): size too small", call. = FALSE)
  }
  
  J <- J-1
  L <- L-1
  N <- max(J,L)+1
  if(!is.null(N)){
    t <- tau*(N^2)
  }
  # calculate transition matrix and sampling matrices
  matP <- TransitionMatrix_BG(N,beta,theta,B1,B2)
  matHK=getH_(N,J)
  matHL=getH_(N,L)
  
  # calculate the marginal likelihood (sum of columns of the maximum likelihood matrix)
  eigenSystem=eigen(t(matP))
  eigenValues=eigenSystem$values
  eigenForw=t(eigenSystem$vectors)
  pi=eigenForw[1,]/sum(eigenForw[1,])
  eigenBack=solve(eigenForw)
  matPt=eigenBack%*%diag(exp(t*log(eigenValues)))%*%eigenForw
  matML=t(matHK)%*%t(matPt)%*%diag(pi)%*%matPt%*%matHL
  vMarg=matML[1,] 
  
  # generate a multinomial sample - so a single spectrum -with probabilities based on the marginal likelihoods
  for(i in 2:length(matML[,1])){
    vMarg=vMarg+matML[i,]  
  }
  vNDotJ=rmultinom(1,size=size,p=vMarg)
  
  # multinomially distribute the single spectrum onto a joint spectrum based on the marginal likelihoods 
  mSim=matrix(nrow=(J+1),ncol=(L+1))
  for(j in 1:ncol(mSim)){
    mSim[,j]=rmultinom(1,size=vNDotJ[j],p=matML[,j])
  }
  return(mSim)
}



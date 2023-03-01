# little programme for inferring:
#     the overall scaled mutation rate, the mutation bias, and the directional components from a single site frequency spectrum 
#     of a population distributed according to the boundary-mutation Moran model
#     (see C.Vogl and L.C.Mikula, Theoretical Population Biology 139 (2021))
#########################################################################################################


# hypergeometrically downsample a vector vMarg to size newM
downSampleV=function(vMarg,newM){
  vMargDown=rep(0,newM+1)
  oldM=length(vMarg)-1
  for(j in 0:oldM){
    vMargDown=vMargDown+table(factor(rhyper(vMarg[j+1],j,oldM-j,newM),levels = 0:newM))
  }
  vMargDown
}

# hypergeometrically downsample the columns of the matrix mMarg to newM
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


# helper function: local search for the maximum-likelihood estimate of the directional component 
getArgMaxGamma_ <- function(mPoly, upper=10){
  M <- length(mPoly)+1
  y <- 1:(M-1)
  vGamma <- seq(0,upper,0.001)
  vlogL <- rep(0,length(vGamma))
  for(i in 1:length(vGamma)){
    temp <- exp(vGamma[i]*y/M)*(1/y+1/(M-y))
    temp <- temp/sum(temp)
    vlogL[i] <- sum(mPoly*log(temp))
  }
  maxgamma <- vGamma[which(vlogL==max(vlogL))]
  return(maxgamma)
}

# helper function: defines the harmonic mean for the directional component 
GenHarmNumM1_=function(gamma,M){
  y <- 1:(M-1)
  return(sum(exp(gamma*y/M)*1/y))
}

# calculates the estimators
getParms <- function(mData){
  if(!is.null(dim(mData))){
    M <- ncol(mData)-1
    ret_est <- matrix(0,nrow=nrow(mData),ncol=5)
    colnames(ret_est) <- c("gamma","rho","vartheta","beta","theta")
    rownames(ret_est) <- rownames(mData)
    for(i in 1:nrow(ret_est)){
      ret_est[i,1] <- getArgMaxGamma_(mData[i,2:M])
      g0 <- GenHarmNumM1_(ret_est[i,1],M)
      g1 <- exp(ret_est[i,1])*GenHarmNumM1_(-ret_est[i,1],M)
      c0 <- g0/(g0+g1)
      ret_est[i,2] <- (mData[i,M+1]+c0*sum(mData[i,2:M]))/sum(mData[i,])
      ret_est[i,3] <- sum(mData[i,2:M])/((g0+g1)*sum(mData[i,]))
      ret_est[i,4] <- ret_est[i,2]/(exp(ret_est[i,1])*(1-ret_est[i,2])+ret_est[i,2]) 
      ret_est[i,5] <- ret_est[i,3]/(ret_est[i,4]*(1-ret_est[i,4]))
    }
    return(ret_est)
  }
  else{
    M <- length(mData)-1
    ret_est <- matrix(0,nrow=1,ncol=5)
    colnames(ret_est) <- c("gamma","rho","vartheta","beta","theta")
     ret_est[1] <- getArgMaxGamma_(mData[2:M])
      g0 <- GenHarmNumM1_(ret_est[1],M)
      g1 <- exp(ret_est[1])*GenHarmNumM1_(-ret_est[1],M)
      c0 <- g0/(g0+g1)
      ret_est[2] <- (mData[M+1]+c0*sum(mData[2:M]))/sum(mData)
      ret_est[3] <- sum(mData[2:M])/((g0+g1)*sum(mData))
      ret_est[4] <- ret_est[2]/(exp(ret_est[1])*(1-ret_est[2])+ret_est[2]) 
      ret_est[5] <- ret_est[3]/(ret_est[4]*(1-ret_est[4]))
    return(ret_est)
  }
}


# calculates the equilibrium distribution from the above estimators
equilDistr_dir <- function(gammaHat,rhoHat,varThetaHat,M){
  M <- M-1
  y <-1:(M-1)
  vals <- matrix(0,nrow=length(gammaHat),ncol=M+1)
  for(i in 1:length(gammaHat)){
    vPoly <- varThetaHat[i]*exp(gammaHat[i]*y/M)*(1/y+1/(M-y))
    vals[i,] <- c(1-rhoHat[i]-varThetaHat[i]*GenHarmNumM1_(gammaHat[i],M),vPoly,rhoHat[i]-varThetaHat[i]*exp(gammaHat[i])*GenHarmNumM1_(-gammaHat[i],M))
  }
  
  return(vals)
}


# Metropolis - Hastings algorithm 
# INPUT: joint site frequency spectrum of two populations assuming equilibrium
#        should be biallelic spectra such as those obtained by contrasting (A+T) vs (C+G) alleles
# OUTPUT: estimators for mutation rate, mutation bias, directional and quadratic selection, and split time between the populations 


# load the source file:
#      source('(INSERT FILE LOCATION)/PopGenJSFS_MCMC_FunctionFile.R', chdir = TRUE)
#############################

J <- 22           # haploid population size of population 1; single (marginal) spectrum is a column vector
L <- 36           # haploid population size of population 2; single (marginal) spectrum is a row vector
theta <- 0.002    # overall scaled mutation rate (scaled by haploid population size)
beta <- 0.33      # mutation bias
B1 <- 1.11        # linear component: directional selection, or biased gene conversion,...
B2 <- 0           # quadratic component: quadratic selection,...
tau <- 0.8        # split time, on the diffusion time scale 
size <- 40000000  # number of sites in the spectrum; NOTE THAT THIS MAY HAVE TO BE QUITE LARGE TO ENSURE INFERENCE ACCURACY

# joint spectrum is simulated :
simMat <- SimJSFS(J=J,L=L,beta=beta,theta=theta,B1=B1,B2=B2,tau=tau,size=size)
# joint spectrum is downsampled to a square matrix 
#   (this step is not necessary, the sampler does it for you if necessary)
simDown <- downSampleM(simMat,23)

# run the Metropolis-Hastings algorithm
#   N - dimension of the square joint spectrum that is the basis of the calculations; when = NULL N=max(J,L)+1 
#   beta - initial value for mutation bias
#   theta - initial value for overall mutation rate
#   B1 - initial value for linear component, default is 0
#   B2 -  initial value for quadratic component, default is 0
#   tau - initial value for split time
#   jSFS - observed spectrum
#   maxIT - number of iterations, default is 20000
#   burnin - length of burnin, default is 40% of maxIT
#   fix.background - keeps beta, theta, tau constant to infer B1 and/or B2 given known background variation, default is FALSE
#   fix.B1 - fixes B1 to infer the background parameters (beta, theta, tau), default is TRUE
#   fix.B2 - fixes B2 to infer the background parameters (beta, theta, tau), default is TRUE
#   tune - tune*(size of the joint spectrum) scale the rate and scale parameters for the proposal distributions for beta (beta distribution) and theta (gamma distribution),
#           while 100/tune and 100000/tune are the variances for the proposal distributions for B1 and B2 (normal)
#           --- note that while the default is 0.2, playing with this can improve convergence behaviour!
#   verbose - when TRUE (default), percent of completed iterations is printed in 10% steps
# -> traces should pop up in the plot window when the iterations are complete

InferParms <-MetrHast(beta=0.2,theta=0.0008,B1=0.6,tau=1.5,jSFS=simDown,fix.B1=FALSE,verbose=FALSE)
# view the results
InferParms$estimates


############################################################################################
# 'intermission' to infer estimators from a site frequency spectrum from a single population 
#       (note quadratic selection and split time estimation is not implemented here!)
# load the source file:
#      source('(INSERT FILE LOCATION)/PopGenSFS_FunctionFile.R', chdir = TRUE)

gamma <- 1.105        # corresponds to B1 in the MCMC
rho <- 0.59           # = beta*exp(gamma) / ((1-beta)+beta*exp(gamma))
vartheta <- 0.000255  # = (beta*(1-beta)*theta) / ((1-beta)+beta*exp(gamma))
M <- 24               # sample size

# simulate an equilibrium spectrum
simVec <- equilDistr_dir(gamma,rho,vartheta,M)
# infer parameters
getParms(simVec)

# end 'intermission'
############################################################################################


############################################################################################
# MCMC example from the Appendix of accompanying paper

J <- 20
L <- 22
beta <- 0.4
theta <- 0.001
B1 <- -(log(0.4))
B2<- 1
tau <- 0.65
size <- 10^{8}

bigSim <- SimJSFS(J=J,L=L,beta=beta,theta=theta,B1=B1,B2=B2,tau=tau,size=size)
compSim1<- SimJSFS(J=J,L=L,beta=0.1,theta=0.002,B1=0.9,B2=0,tau=0.3,size=size)
plot(colSums(bigSim)[2:21],pch="+")
lines(colSums(compSim1)[2:21],pch="x",col="grey")


Infer_bigSim <- MetrHast(beta=0.1,theta=0.002,B1=0.9,tau=0.3,jSFS=bigSim,fix.B1=FALSE,fix.B2=FALSE,verbose=FALSE)
Infer_bigSim$estimates




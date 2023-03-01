## Inference based on the boundary-mutation Moran Model (BMM)

These methods use bi-allelic sample site frequency spectra from two populations (A.) or a single (B.) population.
Assuming equilibrium, mutation and selection parameters are estimated. 

The site frequencies are assumed to evolve according to the boundary-mutation Moran model:
This is an approximation to the bi-allic Moran model with reversible mutations for small scaled overall mutation rates.

Specifically, the methods are based on the approximate equilibrium transition rates in:
Equations  27-32 of C.Vogl and L.C.Mikula, Theoretical Population Biology 139 (2021), p.1-17
<https://www.sciencedirect.com/science/article/pii/S0040580921000289>


Please follow the examples provided in 
<https://LynetteCaitlin.github.io/BMM_MCMC/PopGenExamples.R>
and read the comments for a full understanding of how to apply the methods.

The examples pertain to two inference techniques:

## A. BMM_MCMC 

This is a Metropolis-Hastings algorithm for inferring the overall scaled mutation rate, the mutation bias, and the strength of directional selection (or biased gene conversion), and the strength of quadratic selection from the joint site frequency spectrum of two populations, as well as the population split time. 

The source file for this is:
<https://LynetteCaitlin.github.io/BMM_MCMC/PopGenJSFS_MCMC_FunctionFile.R>

Details about the algorithm can be found in the article (...tba when available...)

## B. Simple Method of Moment Estimators

Basic code for inferring the overall scaled mutation rate, the mutation bias, and the strength of directional selection (or biased gene conversion) from the joint site frequency spectrum of an individual population.

The source file for this is:
<https://LynetteCaitlin.github.io/BMM_MCMC/PopGenSFS_FunctionFile.R>



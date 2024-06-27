# UnSteadyErosion calculator

This toolbox has two main purposes: (1) provide forward models to calculate
nuclide concentrations for unsteady erosion scenarios, and (2) provide inversion
schemes to analyse data for these scenarios.

The current erosion scenarios are:
* a multi-step change in erosion rates. For multiple samples, the step change
occurs at the same time. All erosion rates are allowed to vary between basins

* a same-magnitude multi-step change in erosion. This means that in all basins
and for every step change, erosion rates change by the same magnitude.

## Production rates
Production rates are calculated based on the online calculators formerly known as
CRONUS-Earth online calculators v3. 

## Example scripts
* 'Test_MCMC_inversion' Shows how to use the multi-step change erosion model for 
a parameter inversion. It generates test data and tries to recover the initial
data to illustrate what parameter/sample combinations can theoretically be re-
solved

* 'Test_MCMC_same_magnitude' Shows the same for the same-magnitude change model

* 'WC_MCMC_inversion' shows the use with Crete 14C-10Be data

## Inversion sampler
This toolbox uses the 'gwmcmc' Bayesian Ensemble Sampler (Goodman and Weare, 2010).
This algorithm converges faster on a solution than traditional MCMC samplers in 
high-dimensional parameters space, as it is unaffected by affine tranformations 
of space (linear transformations etc.).

Goodman, J., & Weare, J. (2010). Ensemble samplers with affine invariance. Communications in Applied Mathematics and Computational Science, 5(1), 65-80. https://doi.org/10.2140/camcos.2010.5.65
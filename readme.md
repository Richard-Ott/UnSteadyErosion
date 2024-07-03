# UnSteadyErosion calculator

This toolbox has two main purposes: (1) provide forward models to calculate
nuclide concentrations for unsteady erosion scenarios, and (2) provide inversion
schemes to analyze data for these scenarios. The code is currently adapted to work
for 10Be-14C paired nuclide measurements. However, the forward models can be used for
single nuclide analysis and other nuclides such as Al can be easily added. All models
can be run for a single sample as well as a suite of data. However, the inversion
cannot constrain meaningful results for a single sample, because there are less
data than parameters.

Currently 8 erosion scenarios are supported. There are two main scenario types:
step-changes in erosion and erosion spikes. In step change models, the erosion rate
changes at certain times, whereas in spike models a sudden soil loss is simulated
by removing the upper part of the production profile. For these main scenarios, there
are four sub-scenarios with varying degrees of freedom. All models allow multiple
changes in erosion through time. However, the number of parameters should be smaller
than the number of data points (your nuclide measurements). Therefore, if you want
to resolve a model with more than one change in erosion through time, choose a simpler
erosion scenario. In all scenarios, the timing of change is assumed to be uniform
across all catchments. The scenarios are:

* 'step': a (multi-)step-change in erosion rates at one or more times, where erosion
rates are allowed to vary between catchments. 
* 'samestep': a (multi-)step-change in erosion rates, where background erosion varies 
between catchments, but erosion increases/decreases by a common factor.
* 'samebackground_step': a (multi-)step-change in erosion rates, where all background
erosion is the same, but catchments have different erosion rate changes
* 'samebackground_samestep': a (multi-)step-change in erosion rates, where all background
erosion is the same and all change factors.
* 'spike': a (multi-)spike in soil loss at one or more times, where erosion
rates are allowed to vary between catchments, as well as, the soil loss heights. 
* 'samespike': a (multi-)spike soil loss, where background erosion varies 
between catchments, but the amount of soil loss is the same.
* 'samebackground_spike': a (multi-)spike soil loss, where all background
erosion is the same, but catchments have soil loss.
* 'samebackground_samespike': a (multi-)spike soil loss, where all background
erosion is the same and all soil loss heights.


## Production rates
Production rates are calculated based on the online calculators formerly known as
CRONUS-Earth online calculators v3. 

## Example scripts
* 'Test_MCMC' Script that can generate test data for all erosion scenarios
and then performs the parameter inversion. The inversion tries to recover the initial
data to illustrate what parameter/sample combinations can theoretically be re-
solved.

* 'WC_MCMC_inversion' shows the use with Crete 14C-10Be data

## Inversion sampler
This toolbox uses the 'gwmcmc' Bayesian Ensemble Sampler (Goodman and Weare, 2010).
This algorithm converges faster on a solution than traditional MCMC samplers in 
high-dimensional parameters space, as it is unaffected by affine tranformations 
of space (linear transformations etc.).

Goodman, J., & Weare, J. (2010). Ensemble samplers with affine invariance. Communications in Applied Mathematics and Computational Science, 5(1), 65-80. https://doi.org/10.2140/camcos.2010.5.65
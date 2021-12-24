# Three Phase Analysis
This R code implements a generalized raking and multiple imputation approaches to be applied in 3-phase studies.

There are 2 files:
- 3ph_Weibull_withRakingMI_timevarying: This is the main function. It runs 1,000 Monte Carlo simulations to compare the following estimators: IPW, 2- and 3-phase generalized, 2- and 3-phase MI, and 2- and 3-phase estimators that combined generalized raking and MI.
- functions_simulation: A set of functions that are used to generate data and apply all estimators listed above. This function should be loaded before running the main function

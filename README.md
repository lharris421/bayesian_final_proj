# Bayesian Final Project

This repository contains code for the final project titled "Evaluating the association between CF carrier status and severity of chronic sinusitis using Divide and Conquer Algorithms"

We have structured the repository as follows:

- The folder `logistic_regression` contains all the relevant code for the current project
- Within the `logistic_regression` folder we have the following sub-folders and scripts:
    1) `normal_approx`
        - `normal_approx_full.R`: script for non-D&C approach
        - `normal_approx_dnc.R`: script for D&C w recentering approach
        - `normal_approx_dnc_WASP.R`: script for D&C w WASP approach
        - `normal_approximation_sim_K_plot.R`: script for simulation plots assessing bias and computation time over different values of k
    2) `adaptive_MH`
        - `adaptive_MH_full.R`: script for non-D&C approach
        - `adaptive_MH_dnc.R`: script for D&C w recentering approach
        - `adaptive_MH_dnc_WASP.R`: script for D&C w WASP approach
    3) `MH`
        - `mh_full.R`: script for non-D&C approach
        - `mh_dnc.R`: script for D&C w recentering approach
        - `mh_dnc_timing.R`: script for obtaining run times for D&C approaches
        - `submit_mh_full.R`:
        - `submit_stochastic_logistic_regression.R`
        - `mh_dnc.job`
        - `mh_full.job`
        - `mh_sim_K_plot.R` : script for simulation plots assessing bias and time over different values of k
    4) `partitions`
        - `partition.R`: script to partition the dataset into k partitions
        - `partition.job`: job scirpt to run `partition.R` on HPC
    5) `misc`
        - `likelihoods.md`: markdown of our non-D&C and D&C likelihoods used in the analysis
        - `likelihoods.Rmd`: Rmarkdown to generate `likelihoods.md`
        - `glm_full.R`: Run the frequentist GLM model
        - `time_compare_plot.R`: script to generate plot of time comparisons for the different approaches
    6) `archive`
        - contains old scripts of various different iterations of initial attempts





#!/bin/bash
#$ -N logistic_regression
#$ -q UI
#$ -pe smp 2
#$ -m a
#$ -M loharris@uiowa.edu
#$ -V
#$ -j y

module load stack/2020.1
module load r-stanheaders/2.18.1-10_gcc-9.2.0
module load r-rstan/2.19.2_gcc-9.2.0

R CMD BATCH --no-save --no-restore submit_logistic_regression.R out/logistic_regression.rout
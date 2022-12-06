
library(tidyverse)
library(magrittr)

source('https://myweb.uiowa.edu/pbreheny/7110/f21/notes/fun.R')

# Binomial Binary trial is independently repeated n=20 times with x=6 successes:
n <- 20
x <- 6
theta <- seq(0, 1, len=499)
L <- plotL(theta, dbinom(x, n, theta))
plotL(theta, dbinom(x, n, theta))
plotL(theta, (dbinom(x, n, theta)^10), add = T, col = "red")
legend("topright", 
       col = c("blue", "red"),
       lty = 1,
       lwd = 3,
       legend = c("Likelihood", "Power Likelihood"))

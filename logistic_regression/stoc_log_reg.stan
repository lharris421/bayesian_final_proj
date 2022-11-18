functions {
  real stoc_binomial_lpmf (int y, int ntrail, real p, real nrep) {
    return (nrep * binomial_lpmf(y | ntrail, p));
  }
}

data {
  int<lower=0> n;
  int<lower=0> p;
  int<lower=1> ntrail;
  matrix[n, p] x;
  int<lower=0, upper=ntrail> y[n];
  real nrep;
}

parameters {
  vector[p] beta;
}

transformed parameters {
  vector[n] odds;
  vector[n] prob;

  odds = exp(x * beta);

  for (ii in 1:n) {
    prob[ii] = odds[ii] / (1 + odds[ii]);
  }  
}

model {
  
  for (ii in 1:n) {
    y[ii] ~ stoc_binomial(ntrail, prob[ii], nrep);
  }  
}
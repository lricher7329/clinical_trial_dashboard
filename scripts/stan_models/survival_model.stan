// Weibull Survival Model with Treatment Effect
data {
  int<lower=0> N;             // number of observations
  int<lower=0> K;             // number of predictors
  matrix[N, K] X;             // predictor matrix
  vector<lower=0>[N] time;    // observed event times
  vector<lower=0, upper=1>[N] event;  // censoring indicator (1=event, 0=censored)
}

parameters {
  vector[K] beta;             // coefficients for predictors
  real<lower=0> alpha;        // shape parameter for Weibull distribution
}

model {
  // Priors
  beta ~ normal(0, 2.5);
  alpha ~ gamma(2, 0.1);
  
  // Likelihood
  for (i in 1:N) {
    // Linear predictor
    real eta = X[i] * beta;
    
    // Weibull log-likelihood with censoring
    if (event[i] == 1)
      target += weibull_lpdf(time[i] | alpha, exp(-eta/alpha));
    else
      target += weibull_lccdf(time[i] | alpha, exp(-eta/alpha));
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N] hazard_ratio;
  
  for (i in 1:N) {
    real eta = X[i] * beta;
    real lambda = exp(-eta/alpha);
    
    // Log-likelihood for LOO-CV
    if (event[i] == 1)
      log_lik[i] = weibull_lpdf(time[i] | alpha, lambda);
    else
      log_lik[i] = weibull_lccdf(time[i] | alpha, lambda);
      
    // Hazard ratio for treatment effect
    hazard_ratio[i] = exp(beta[2]); // Assuming treatment is 2nd predictor
  }
}

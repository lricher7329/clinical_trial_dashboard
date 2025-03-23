// Longitudinal mixed effects model for biomarker data
data {
  int<lower=0> N;                  // number of observations
  int<lower=0> J;                  // number of subjects
  int<lower=0> K;                  // number of fixed effect predictors
  array[N] int<lower=1, upper=J> subject; // subject ID for each observation
  matrix[N, K] X;                  // fixed effect design matrix
  vector[N] y;                     // outcome variable (biomarker)
  vector[N] time;                  // time variable
}

parameters {
  vector[K] beta;                  // fixed effect coefficients
  real<lower=0> sigma;             // residual standard deviation
  
  vector[J] u_intercept;           // subject-specific random intercepts
  vector[J] u_slope;               // subject-specific random slopes
  
  vector<lower=0>[2] tau;          // random effect standard deviations
  real<lower=-1, upper=1> rho;     // correlation between random effects
}

transformed parameters {
  matrix[2, 2] Omega;              // random effects correlation matrix
  matrix[2, 2] Sigma;              // random effects covariance matrix
  
  // Construct correlation matrix
  Omega[1, 1] = 1;
  Omega[2, 2] = 1;
  Omega[1, 2] = rho;
  Omega[2, 1] = rho;
  
  // Construct covariance matrix
  Sigma = diag_matrix(tau) * Omega * diag_matrix(tau);
}

model {
  // Priors for fixed effects
  beta ~ normal(0, 5);
  
  // Priors for variance components
  sigma ~ cauchy(0, 2.5);
  tau ~ cauchy(0, 2.5);
  rho ~ normal(0, 0.5);
  
  // Random effects distribution
  for (j in 1:J) {
    [u_intercept[j], u_slope[j]]' ~ multi_normal(rep_vector(0, 2), Sigma);
  }
  
  // Likelihood
  {
    vector[N] mu;
    for (i in 1:N) {
      mu[i] = X[i] * beta + u_intercept[subject[i]] + u_slope[subject[i]] * time[i];
    }
    y ~ normal(mu, sigma);
  }
}

generated quantities {
  // Predicted values
  vector[N] y_pred;
  vector[N] log_lik;
  
  // Treatment effect at different time points
  real effect_baseline = beta[2]; // Assuming treatment is the 2nd predictor
  real effect_endpoint = beta[2] + beta[4]; // Assuming treatment*time is the 4th predictor
  
  for (i in 1:N) {
    real mu = X[i] * beta + u_intercept[subject[i]] + u_slope[subject[i]] * time[i];
    y_pred[i] = normal_rng(mu, sigma);
    log_lik[i] = normal_lpdf(y[i] | mu, sigma);
  }
}

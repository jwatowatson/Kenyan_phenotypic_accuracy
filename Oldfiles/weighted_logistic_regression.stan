data {
  int<lower=0> N;
  int<lower=2> K;
  int<lower=0,upper=1> case_control[N];
  vector<lower=0,upper=2>[N] g1;
  vector<lower=0,upper=1>[N] ws;
  matrix[N,K] xx;
}
parameters {
  real alpha;
  vector[K] beta;
  real beta0;
}
model {
  alpha ~ normal(0,1);
  beta ~ normal(0,1);

  for(i in 1:N){
    target += bernoulli_logit_lpmf(case_control[i] | xx[i] * beta + alpha*g1[i] + beta0)*ws[i];
  }
}


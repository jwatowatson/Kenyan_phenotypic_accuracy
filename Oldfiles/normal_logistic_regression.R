require(rstan)

standard_logit_mod = '
data {
  int<lower=0> N;
  int<lower=0,upper=1> case_control[N];
  vector[N] g;
  int<lower=0> K;
  matrix[N,K] x;
}
parameters {
  real alpha;
  real beta_g;
  vector[K] beta_x;
}
model {
  //alpha ~ normal(0,1);
  //beta_g ~ normal(0,1);
  //beta_x ~ normal(0,1);
  case_control ~ bernoulli_logit(beta_g*g + x*beta_x + alpha);
}'

weighted_logit_mod = '
data {
  int<lower=0> N;
  int<lower=0,upper=1> case_control[N];
  vector[N] g;
  vector[N] ws;
  int<lower=0> K;
  matrix[N,K] x;
}
parameters {
  real alpha;
  real beta_g;
  vector[K] beta_x;
}
model {
  //alpha ~ normal(0,1);
  //beta_g ~ normal(0,1);
  //beta_x ~ normal(0,1);
  for(n in 1:N){
    target += bernoulli_logit_lpmf(case_control[n] | beta_g*g[n] + x[n]*beta_x + alpha)*ws[n];
  }
}
'

normal_glm = stan_model(model_code = standard_logit_mod)


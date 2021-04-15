data {
  int<lower=0> Ncases;
  int<lower=0> Ncontrols;
  int<lower=2> K;
  int<lower=0,upper=2> g_cases[Ncases];
  int<lower=0,upper=2> g_controls[Ncontrols];
  vector<lower=0,upper=1>[Ncases] ps;
  matrix[Ncases,K] x_cases;
  matrix[Ncontrols,K] x_controls;
}
transformed data {
  int y_controls[Ncontrols];
  for(i in 1:Ncontrols){
    y_controls[i] = 0;
  }
}
parameters {
  real alpha_P1;
  real alpha_P2;
  vector[K] beta;
}
model {
  alpha_P1 ~ normal(0,1);
  alpha_P2 ~ normal(0,1);
  beta ~ normal(0,1);
  g_controls ~ binomial_logit(2, x_controls * beta);
  for(n in 1:Ncases){
    target += log_sum_exp(log(ps[n]) +
    binomial_logit_lpmf(g_cases[n] | 2, x_cases[n] * beta + alpha_P1),
    log1m(ps[n]) + binomial_logit_lpmf(g_cases[n] | 2, x_cases[n] * beta + alpha_P2));
  }
}


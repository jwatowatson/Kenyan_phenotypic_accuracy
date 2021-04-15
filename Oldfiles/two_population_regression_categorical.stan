data {
  int<lower=0> Ncases;
  int<lower=0> Ncontrols;
  int<lower=2> K;
  vector<lower=0,upper=2>[Ncases] g_cases;
  vector<lower=0,upper=2>[Ncontrols] g_controls;
  vector<lower=0,upper=1>[Ncases] ps;
  matrix[Ncases,K] x_cases;
  matrix[Ncontrols,K] x_controls;
}
transformed data {
  int y_controls[Ncontrols];
  for(i in 1:Ncontrols){
    y_controls[i] = 1;
  }
}
parameters {
  real alpha;
  real alpha_P1;
  real alpha_P2;
  vector[K] beta;
}
model {
  alpha_P1 ~ normal(0,1);
  alpha_P2 ~ normal(0,1);
  beta ~ normal(0,1);
  y_controls ~ categorical_logit(x_controls * beta + alpha*g_controls);
  for(n in 1:Ncases){
    target += log_sum_exp(log(ps[n]) +
    categorical_logit_lpmf(1 | x_cases[n] * beta + alpha_P1),
    log1m(ps[n]) + categorical_logit_lpmf(2 | x_cases[n] * beta + alpha_P2));
  }
}


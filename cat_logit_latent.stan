data {
  int Ncases;
  int Ncontrols;
  int D;
  vector<lower=0,upper=1>[Ncases] ps;
  matrix[Ncases, D] x_cases;
  matrix[Ncontrols, D] x_controls;
}
transformed data {
  vector[D] zeros = rep_vector(0, D);
}
parameters {
  matrix[D, 2] beta_raw;
}
transformed parameters {
  matrix[D, 3] beta;
  beta = append_col(beta_raw, zeros);
}
model {
  matrix[Ncases, 3] x_beta_cases = x_cases * beta;
  matrix[Ncontrols, 3] x_beta_controls = x_controls * beta;

  to_vector(beta_raw) ~ normal(0, 5);

  for(n in 1:Ncontrols){
    target += categorical_logit_lpmf(3 | x_beta_controls[n]');
  }

  for(n in 1:Ncases){
    target += log_sum_exp(log(ps[n]) + categorical_logit_lpmf(1 | x_beta_cases[n]'),
                          log1m(ps[n]) + categorical_logit_lpmf(2 | x_beta_cases[n]'));
  }
}

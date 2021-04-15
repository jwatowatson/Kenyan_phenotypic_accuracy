require(rstan)
mixture_model_student = '
  data {
      int<lower=1> D;                             // number of dimensions
      int<lower=1> N_Training;                    // Training data size
      int<lower=1> N_Testing;                     // Testing data size
      int<lower=1> K_notSM;
      vector[D] y_Training[N_Training];           // data
      vector[D] y_Testing[N_Testing];             // data with unobserved labels
      vector[D] mu_priorSM;
      vector[D] mu_prior_notSM;
      vector[D] sd_priorSM;
      vector[D] sd_prior_notSM;
      real<lower=0> prior_beta1;
      real<lower=0> prior_beta2;
  }
  transformed data {
    vector[K_notSM] alpha;
    for(k in 1:K_notSM){
      alpha[k] = 0.1;
    }
  }
  parameters {
    real<lower=0,upper=1> theta_SM;
    simplex[K_notSM] theta_notSM;               // mixture proportions Training data
    vector[D] mu[1+K_notSM];
    cholesky_factor_corr[D] L[1+K_notSM];    // cholesky factor of covariance
    vector<lower=0>[D] Sigma_scale[1+K_notSM];
  }
  transformed parameters {
    matrix[D, D] Sigma[1+K_notSM];
    // covariance matrices
    for(k in 1:(1+K_notSM)){
      Sigma[k] = crossprod(diag_pre_multiply(Sigma_scale[k], L[k]));
    }
  }
  model {
      // prior
      mu[1] ~ normal(mu_priorSM, sd_priorSM);
      L[1] ~ lkj_corr_cholesky(1);
      Sigma_scale[1] ~ normal(0, .25);

      for(k in 2:(1+K_notSM)){
        mu[k] ~ normal(mu_prior_notSM, sd_prior_notSM);
        L[k] ~ lkj_corr_cholesky(1);
        Sigma_scale[k] ~ normal(0, .25);
      }
      theta_notSM ~ dirichlet(alpha);
      theta_SM ~ beta(prior_beta1, prior_beta2);

      // likelihood - training data with known labels
      y_Training ~ multi_student_t(7, mu[1], Sigma[1]);

      // likelihood - testing data with unknown labels
      for (n in 1:N_Testing){
        vector[1+K_notSM] lps;
        lps[1] = log(theta_SM)+multi_student_t_lpdf(y_Testing[n] | 7, mu[1], Sigma[1]);
        for(k in 2:(1+K_notSM)){
          lps[k] = log1m(theta_SM) + log(theta_notSM[k-1]) + multi_student_t_lpdf(y_Testing[n] | 7, mu[k], Sigma[k]);
        }
        target += log_sum_exp(lps);
      }

  }
  generated quantities {
    vector[N_Testing] prob_SM;
    // This computes respective densities of each mixture component
    for(n in 1:N_Testing){
      vector[1+K_notSM] probs_raw;
      probs_raw[1] = theta_SM*exp(multi_student_t_lpdf(y_Testing[n] | 7, mu[1], Sigma[1]));
      for(k in 2:(1+K_notSM)){
          probs_raw[k] = (1-theta_SM)*theta_notSM[k-1]*exp(multi_student_t_lpdf(y_Testing[n] | 7, mu[k], Sigma[k]));
      }
      // normalise to 1
      prob_SM[n] =  probs_raw[1]/sum(probs_raw);
    }
  }

'


mix_model_student = stan_model(model_code = mixture_model_student)

require(rstan)
mixture_model3 = '
  data {
      int<lower=1> D;                             // number of dimensions
      int<lower=1> N_Training;                    // Training data size
      int<lower=1> N_Testing;                     // Testing data size
      int<lower=1> K_notSM;
      int<lower=1> K_SM;
      vector[D] y_Training[N_Training];           // data
      vector[D] y_Testing[N_Testing];             // data with unobserved labels
      vector[D] mu_prior[K_SM+K_notSM];
      real<lower=0> prior_beta1;
      real<lower=0> prior_beta2;
  }
  parameters {
      simplex[K_SM+K_notSM] theta_Testing;        // mixture proportions Kenya
      simplex[K_SM] theta_Training;               // mixture proportions Training data
      vector[D] mu[K_SM+K_notSM];
      cholesky_factor_corr[D] L[K_SM+K_notSM];    // cholesky factor of covariance
      vector<lower=0>[D] Sigma_scale[K_SM+K_notSM];
  }
  transformed parameters {
    matrix[D, D] Sigma[K_SM+K_notSM];
    // covariance matrices
    for(k in 1:(K_SM+K_notSM)){
      Sigma[k] = crossprod(diag_pre_multiply(Sigma_scale[k], L[k]));
    }
  }
  model {
      vector[K_SM+K_notSM] log_theta_Testing = log(theta_Testing);  // cache log calculation
      vector[K_SM] log_theta_Training = log(theta_Training);  // cache log calculation
      // prior
      for(k in 1:(K_SM+K_notSM)){
        mu[k] ~ normal(mu_prior[k], .1);
        L[k] ~ lkj_corr_cholesky(1);
        Sigma_scale[k] ~ normal(0, .25);
      }
      theta_Testing[1] ~ beta(prior_beta1, prior_beta2);

      // likelihood - training data with known labels
      if(K_SM==1){
        y_Training ~ multi_normal(mu[1], Sigma[1]);
      } else {
        for (n in 1:N_Training){
          vector[K_SM] lps = log_theta_Training;
          for(k in 1:K_SM){
            lps[k] += multi_normal_lpdf(y_Training[n] | mu[k], Sigma[k]);
          }
          target += log_sum_exp(lps);
        }
      }

      // likelihood - testing data with unknown labels
      for (n in 1:N_Testing){
        vector[K_SM+K_notSM] lps = log_theta_Testing;
        for(k in 1:(K_SM+K_notSM)){
          lps[k] += multi_normal_lpdf(y_Testing[n] | mu[k], Sigma[k]);
        }
        target += log_sum_exp(lps);
      }

  }
  generated quantities {
    vector[K_SM+K_notSM] prob_SM[N_Testing];
    // This computes respective densities of each mixture component
    for(n in 1:N_Testing){
      vector[K_SM+K_notSM] probs_raw;
      for(k in 1:(K_SM+K_notSM)){
          probs_raw[k] = exp(theta_Testing[k])*exp(multi_normal_lpdf(y_Testing[n] | mu[k], Sigma[k]));
      }
      // normalise to 1
      prob_SM[n] =  probs_raw/sum(probs_raw);
    }
  }

'


mix_model3 = stan_model(model_code = mixture_model3)

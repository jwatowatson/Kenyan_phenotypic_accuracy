require(rstan)
mixture_model2 = '
  data {
      int D;                            // number of dimensions
      int N_SM;                         // number of data with unobserved labels
      int N_latent;                     // number of data with observed labels
      vector[D] y_SM[N_SM];             // data
      vector[D] y_latent[N_latent];     // data with unobserved labels
      vector[D] mu_prior_SM;
      vector[D] mu_prior_not_SM;
      real prior_beta1;
      real prior_beta2;
   //   real theta;
  }
  parameters {
      real<lower=0,upper=1> theta;      // probability SM
      vector[D] mu_SM;
      vector[D] mu_not_SM;
      cholesky_factor_corr[D] L_SM;     // cholesky factor of covariance
      cholesky_factor_corr[D] L_not_SM; // cholesky factor of covariance
      vector<lower=0>[D] Sigma_scale_SM;
      vector<lower=0>[D] Sigma_scale_not_SM;
  }
  transformed parameters {
    matrix[D, D] Sigma_SM;
    matrix[D, D] Sigma_not_SM;
    // covariance matrix
    Sigma_SM = crossprod(diag_pre_multiply(Sigma_scale_SM, L_SM));
    Sigma_not_SM = crossprod(diag_pre_multiply(Sigma_scale_not_SM, L_not_SM));
  }
  model {
      real ps[2];
      // prior
      mu_SM ~ normal(mu_prior_SM,.1);
      mu_not_SM ~ normal(mu_prior_not_SM, .2);
      L_SM ~ lkj_corr_cholesky(10);
      L_not_SM ~ lkj_corr_cholesky(1);
      Sigma_scale_SM ~ normal(0, .25);
      Sigma_scale_not_SM ~ normal(0, .25);
      theta ~ beta(prior_beta1,prior_beta2);

      // likelihood
      y_SM ~ multi_normal(mu_SM, Sigma_SM);
      for (n in 1:N_latent){
        ps[1] = log(theta) + multi_normal_lpdf(y_latent[n] | mu_SM, Sigma_SM);
        ps[2] = log1m(theta) + multi_normal_lpdf(y_latent[n] | mu_not_SM, Sigma_not_SM);
        target += log_sum_exp(ps);
      }
  }
  generated quantities {
    vector[N_latent] prob_SM;
    // This computes respective densities of each mixture component
    for(n in 1:N_latent){
      vector[2] prob_labels_raw;
      // SM
      prob_labels_raw[1] = theta*exp(multi_normal_lpdf(y_latent[n] | mu_SM, Sigma_SM));
      prob_labels_raw[2] = exp(1-theta)*exp(multi_normal_lpdf(y_latent[n] | mu_not_SM, Sigma_not_SM));
      // normalise to 1
      prob_SM[n] =  prob_labels_raw[1]/sum(prob_labels_raw);
    }
  }

'


mix_model2 = stan_model(model_code = mixture_model2)

require(rstan)
mixture_model_student2 = '
  data {
      int<lower=1> D;                       // number of dimensions
      int<lower=1> N_Training_Africa;       // Training data size
      int<lower=1> N_Training_Asia;         // Training data size
      int<lower=1> N_Testing;               // Testing data size
      int<lower=1> K_notSM;                 // Number of mixture components in notSM
      vector[D] y_Training_Africa[N_Training_Africa];   // data
      vector[D] y_Training_Asia[N_Training_Asia];       // data
      vector[D] y_Testing[N_Testing];                   // data with unobserved labels
      vector[D] mu_priorSM;
      vector[D] mu_prior_notSM;
      vector[D] sd_priorSM;
      vector[D] sd_prior_notSM;
      real<lower=0> prior_beta1;
      real<lower=0> prior_beta2;
  }
  transformed data {
    // make a vector for the Dirichlet prior
    vector[K_notSM] alpha;
    for(k in 1:K_notSM){
      alpha[k] = 1.0 / K_notSM;
    }
  }
  parameters {
    real<lower=0,upper=1> theta_SM;             // mixture proportion in training data
    real<lower=0,upper=1> theta_SMstudy;        // mixture proportions Training data
    simplex[K_notSM] theta_notSM;               // mixture proportions Training data
    vector[D] mu_SM;                            // hierarchical mean in SM
    vector[D] mu_SM_study[2];                   // means in SM in the two training datasets
    vector[D] mu_notSM[K_notSM];
    real<lower=0> sd_muSM;
    cholesky_factor_corr[D] L_SM[2];              // cholesky factor of covariance
    cholesky_factor_corr[D] L_notSM[K_notSM];     // cholesky factor of covariance
    vector<lower=0>[D] Sigma_scale_SM[2];
    vector<lower=0>[D] Sigma_scale_notSM[K_notSM];
  }
  transformed parameters {
    matrix[D, D] Sigma_notSM[K_notSM];
    matrix[D, D] Sigma_SM[2];
    // covariance matrices
    for(k in 1:K_notSM){
      Sigma_notSM[k] = crossprod(diag_pre_multiply(Sigma_scale_notSM[k], L_notSM[k]));
    }
    for(k in 1:2){
      Sigma_SM[k] = crossprod(diag_pre_multiply(Sigma_scale_SM[k], L_SM[k]));
    }
  }
  model {
      // prior
      mu_SM ~ normal(mu_priorSM, sd_priorSM);
      sd_muSM ~ normal(0,.1);
      for(k in 1:2){
        mu_SM_study[k] ~ normal(mu_SM, sd_muSM);
        L_SM[k] ~ lkj_corr_cholesky(1);
        Sigma_scale_SM[k] ~ normal(0, .25);
      }

      for(k in 1:(K_notSM)){
        mu_notSM[k] ~ normal(mu_prior_notSM, sd_prior_notSM);
        L_notSM[k] ~ lkj_corr_cholesky(1);
        Sigma_scale_notSM[k] ~ normal(0, .25);
      }
      theta_notSM ~ dirichlet(alpha);
      theta_SM ~ beta(prior_beta1, prior_beta2);
      theta_SMstudy ~ beta(2,2);

      // likelihood - training data with known labels
      y_Training_Asia ~ multi_student_t(7, mu_SM_study[1], Sigma_SM[1]);
      y_Training_Africa ~ multi_student_t(7, mu_SM_study[2], Sigma_SM[2]);

      // likelihood - testing data with unknown labels
      for (n in 1:N_Testing){
        vector[2+K_notSM] lps;
        lps[1] = log(theta_SM)+log(theta_SMstudy)+multi_student_t_lpdf(y_Testing[n] | 7, mu_SM_study[1], Sigma_SM[1]);
        lps[2] = log(theta_SM)+log1m(theta_SMstudy)+multi_student_t_lpdf(y_Testing[n] | 7, mu_SM_study[2], Sigma_SM[2]);

        for(k in 1:K_notSM){
          lps[k+2] = log1m(theta_SM) + log(theta_notSM[k]) + multi_student_t_lpdf(y_Testing[n] | 7, mu_notSM[k], Sigma_notSM[k]);
        }
        target += log_sum_exp(lps);
      }

  }
  generated quantities {
    vector[N_Testing] prob_SM;
    // This computes respective densities of each mixture component
    for(n in 1:N_Testing){
      vector[2+K_notSM] probs_raw;
      probs_raw[1] = theta_SM*theta_SMstudy*exp(multi_student_t_lpdf(y_Testing[n] | 7, mu_SM_study[1], Sigma_SM[1]));
      probs_raw[2] = theta_SM*(1-theta_SMstudy)*exp(multi_student_t_lpdf(y_Testing[n] | 7, mu_SM_study[2], Sigma_SM[2]));
      for(k in 1:K_notSM){
          probs_raw[k+2] = (1-theta_SM)*theta_notSM[k]*exp(multi_student_t_lpdf(y_Testing[n] | 7, mu_notSM[k], Sigma_notSM[k]));
      }
      // normalise to 1
      prob_SM[n] =  (probs_raw[1]+probs_raw[2])/sum(probs_raw);
    }
  }

'


mix_model_student2 = stan_model(model_code = mixture_model_student2)

require(rstan)
univariate_mixture_model_studentt = '
  data {
      int<lower=1> N_Training_SM;           // SM data size
      int<lower=1> N_Training_UM;           // UM data size
      int<lower=1> N_Unknown;               // SM/UM mix data size
      int<lower=1> K_UM;                    // Number of mixture components in UM
      int<lower=1> K_SM;                    // Number of mixture components in SM
      vector[N_Training_SM] y_Training_SM;  // data
      vector[N_Training_UM] y_Training_UM;  // data
      vector[N_Unknown] y_Unknown;          // data with unobserved labels
      real mu_prior_SM[K_SM];
      real mu_prior_UM[K_UM];
      real sd_prior_SM;
      real sd_prior_UM;
      real<lower=0> prior_beta1;
      real<lower=0> prior_beta2;
  }
  transformed data {
    // make a vector for the Dirichlet prior
    vector[K_UM] alpha_UM;
    vector[K_SM] alpha_SM;
    for(k in 1:K_UM){
      alpha_UM[k] = 1.0;
    }
    for(k in 1:K_SM){
      alpha_SM[k] = 1.0;
    }
  }
  parameters {
    real<lower=0,upper=1> P_SM;             // mixture proportion in unknown data
    simplex[K_UM] theta_UM;                     // mixture proportions UM data
    simplex[K_SM] theta_SM;                     // mixture proportions UM data
    real<lower=0> mu_SM[K_SM];                            // hierarchical mean in SM
    real<lower=0> mu_UM[K_UM];
    real<lower=0> sd_SM[K_SM];
    real<lower=0> sd_UM[K_UM];
  }
  model {
      // prior
      for(k in 1:(K_SM)){
        mu_SM[k] ~ normal(mu_prior_SM[k], sd_prior_SM);
        sd_SM[k] ~ normal(0, 1);
      }

      for(k in 1:(K_UM)){
        mu_UM[k] ~ normal(mu_prior_UM[k], sd_prior_UM);
        sd_UM[k] ~ normal(0, 1);
      }

      theta_UM ~ dirichlet(alpha_UM);
      theta_SM ~ dirichlet(alpha_SM);
      P_SM ~ beta(prior_beta1, prior_beta2);

      // likelihood - training data with known labels
      for(n in 1:N_Training_SM){
        vector[K_SM] lps;
        for(k in 1:K_SM){
          lps[k] = log(theta_SM[k]) + normal_lpdf(y_Training_SM[n] | mu_SM[k], sd_SM[k]);
        }
        target += log_sum_exp(lps);
      }

      for(n in 1:N_Training_UM){
        vector[K_UM] lps;
        for(k in 1:K_UM){
          lps[k] = log(theta_UM[k]) + normal_lpdf(y_Training_UM[n] | mu_UM[k], sd_UM[k]);
        }
        target += log_sum_exp(lps);
      }

      // likelihood - testing data with unknown labels
      for (n in 1:N_Unknown){
        vector[K_SM+K_UM] lps;
        for(k in 1:K_SM){
          lps[k] = log(P_SM) + log(theta_SM[k]) + normal_lpdf(y_Unknown[n] | mu_SM[k], sd_SM[k]);
        }
        for(k in 1:K_UM){
          lps[k+K_SM] = log1m(P_SM) + log(theta_UM[k]) + normal_lpdf(y_Unknown[n] | mu_UM[k], sd_UM[k]);
        }
        target += log_sum_exp(lps);
      }

  }
  generated quantities {
    vector[N_Unknown] prob_SM;
    // This computes respective densities of each mixture component
    for(n in 1:N_Unknown){
      vector[K_SM+K_UM] probs_raw;
      for(k in 1:K_SM){
          probs_raw[k] = P_SM*theta_SM[k]*exp(normal_lpdf(y_Unknown[n] | mu_SM[k], sd_SM[k]));
      }
      for(k in 1:K_UM){
          probs_raw[k+K_SM] = (1-P_SM)*theta_UM[k]*exp(normal_lpdf(y_Unknown[n] | mu_UM[k], sd_UM[k]));
      }
      // normalise to 1
      prob_SM[n] = sum(probs_raw[1:K_SM])/sum(probs_raw);
    }
  }

'


univariate_mix_model_studentt = stan_model(model_code = univariate_mixture_model_studentt)

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> Nunknown;
  int<lower=0> N1;
  int<lower=0> N2;
  int<lower=0> N1_sickle;
  int<lower=0> N2_sickle;
  int<lower=0> Nunknown_sickle;
  real b1;
  real b2;
  real a1;
  real a2;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0,upper=1> q;
  real<lower=0,upper=1> p1;
  real<lower=0,upper=1> p2;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  p1 ~ beta(b1, b2);
  p2 ~ beta(a1, a2);
  q ~ beta(1,1);
  Nunknown_sickle ~ binomial(Nunknown, q*p1 + (1-q)*p2);
  N1_sickle ~ binomial(N1, p1);
  N2_sickle ~ binomial(N2, p2);
}

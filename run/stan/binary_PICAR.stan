/*
  * Binary Data Example
*/
  
  data {
    int N; //the number of observations
    int K; //the number of columns in the model matrix
    int P; //the number of columns in the basis matrix
    matrix[N,K] X; //the model matrix
    matrix[N,P] M; //the basis matrix
    matrix[P,P] MQM; //the precision matrix
    int<lower=0,upper=1> y[N];
  }
  transformed data {
  row_vector[P] zeros;              // create reference level coefs of zero
  zeros = rep_row_vector(0, P);
}

parameters {
  vector[K] beta; //the regression parameters
  vector[P] delta; //the reparameterized random effects
  real<lower=0>  tau; //the precision parameter
}

transformed parameters {
  vector[N] linpred;
  linpred = X*beta+M*delta;
}
model {  
  beta[1] ~ normal(0,100); //prior for the B1 following Gelman 2008
  beta[2] ~ normal(0,100); //prior for the B2 following Gelman 2008
  beta[3] ~ normal(0,100); //prior for the B3 following Gelman 2008
  tau ~ gamma(0.5,2000); //prior for the tau2 for new spatial random effects
  delta ~ multi_normal_prec(zeros, tau * MQM); // MVN prior for the new spatial random effects
  y ~ bernoulli_logit(linpred);
}


/*
  *Count Data with spatially-varying coefficients example
*/
  
  data {
    int N; //the number of observations
    int K; //the number of columns in the model matrix
    int P; //the number of columns in the basis matrix
    matrix[N,K] X; //the model matrix
    matrix[N,P] M; //the basis matrix
    matrix[P,P] MQM; //the basis matrix
    int<lower=0> y[N]; //the response
  }
  transformed data {
  row_vector[P] zeros;              // create reference level coefs of zero
  zeros = rep_row_vector(0, P);
}

parameters {
  vector[K] beta; //the regression parameters
  vector[P] delta; //the reparameterized random effects
  vector[P] delta_B; //Spatially Varying Coefficients
  real<lower=0>  tau; //the precision parameter
  real<lower=0>  tau_B; //the precision parameter for SVCs
}

transformed parameters {
  vector[N] linpred;
  linpred = X*beta+M*delta+ X[,1].*(M*delta_B);
}
model {  
  beta[1] ~ normal(0,100); //prior for the B1 following Gelman 2008
  beta[2] ~ normal(0,100); //prior for the B2 following Gelman 2008
  beta[3] ~ normal(0,100); //prior for the B3 following Gelman 2008
  tau ~ gamma(0.5,2000); //prior for the tau2 for new spatial random effects
  tau_B ~ gamma(0.5,2000); //prior for the tau2 for the spatially varying coefficients
  delta ~ multi_normal_prec(zeros, tau * MQM); // MVN prior for the new spatial random effects
  delta_B ~ multi_normal_prec(zeros, tau_B * MQM); // MVN prior for the spatially varying coefficients
  y ~ poisson_log(linpred);
}

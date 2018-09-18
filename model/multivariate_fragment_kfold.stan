/*Multivariate normal model to multifunctional dataset
*using a cholesky decomposition for the covariance matrix
*with a fragment-level hierarchical intercept
*and some functionalities for cross-validations
*this model was tested and works on simulated data
*/

data {
	int<lower=1> N; //number of rows
	int<lower=1> K; //number of columns in the response matrix
	int<lower=1> J; //number of columns in the X matrix
	int<lower=1> F; //number of fragments

	int<lower=1,upper=F> F_id[N]; //fragment index

	int<lower=0,upper=1> holdout[N]; //index whether the observation should be held out (1) or used (0)

	matrix[N, J] X; //matrix of predictors
	vector[K] y[N]; //matrix of response values
}

parameters {
	matrix[J, K] beta; //regression parameters
	cholesky_factor_corr[K] L_omega; //covariance terms (off-diagonal)
	vector<lower=0>[K] L_sigma; //variance terms (diagonals) 
	real<lower=0> sigma_fragment; //variation in fragment-level intercept one value for all functions
	//vector<lower=0>[K] sigma_fragment; //one value per function

	matrix[F,K] z;
}

transformed parameters {	
	matrix[K,K] L_Sigma; //the cholesky decomposed covariance matrix	
	
	L_Sigma = diag_pre_multiply(L_sigma,L_omega); //re-constructing the covariance matrix	
}

model {	
	matrix[N, K] Xbeta_temp = X * beta; //temporary container to speed up computation
	//matrix[N, K] Xbeta_temp = X * beta; //temporary container to speed up computation
	row_vector[K] mu[N]; //the linear predictor

	
	to_vector(z) ~ normal(0,1); //random normal deviates for fragment-level variation

	//priors
	to_vector(beta) ~ normal(0,5);
	L_omega ~ lkj_corr_cholesky(2);
	L_sigma ~ cauchy(0,2.5);
	sigma_fragment ~ normal(0,1);

	//likelihood on the Kth-fold
	for(n in 1:N){
		mu[n] = Xbeta_temp[n] + z[F_id[n]] * rep_vector(sigma_fragment, K); // compute the matrix of linear predictor

		if(holdout[n] == 0){
			target += multi_normal_cholesky_lpdf(y[n] | mu[n], L_Sigma);
		}
	}
}

generated quantities {
	//matrix[N, K] y_post;//for posterior predictive checks
	vector[N] log_lik;
	matrix[K, K] Sigma; //the variance-covariance matrix
	Sigma = tcrossprod(L_Sigma);
	for (n in 1:N){
		//y_post[n] = to_row_vector(multi_normal_cholesky_rng(to_vector(X[n,] * beta),L_Sigma)); //note that the fragment-level deviation is not taken into account there
                log_lik[n] = multi_normal_lpdf(y[n] | X[n,] * beta, Sigma);
	}
}




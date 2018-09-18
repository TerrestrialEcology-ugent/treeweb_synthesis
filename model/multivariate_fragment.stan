/*Multivariate normal model to multifunctional dataset
*using a cholesky decomposition for the covariance matrix
*with a fragment-level hierarchical intercept
*this model was tested and works on simulated data
*/

functions {

	matrix revv(matrix y_new, vector direction, int N, int K) {
		matrix[K, N] new_matrix;
		//get proper direction (minimize/maximize) so loop through the variables (columns in y_pred)
		//and put them rowise in new_matrix
		for(k in 1:K){
			if(direction[k] == -1)
				new_matrix[k] =  to_row_vector((- col(y_new, k)) + (max(col(y_new, k)) + min(col(y_new, k))));
			else
				new_matrix[k] = to_row_vector(col(y_new, k));
		}
		return new_matrix;

	} 

	vector desirability(matrix new_matrix, vector importance, int N, int K) {		
		vector[N] avg_weighted;
		
		//compute weighted average for each plot
		for(n in 1:N){
			avg_weighted[n] = sum(importance .* col(new_matrix, n)) / sum(importance);
		}
		return avg_weighted;
	}
}


data {
	int<lower=1> N; //number of rows
	int<lower=1> N_pred; //number of rows in the new X matrix used for prediction
	int<lower=1> K; //number of columns in the response matrix
	int<lower=1> J; //number of columns in the X matrix
	int<lower=1> F; //number of fragments

	int<lower=1,upper=F> F_id[N]; //fragment index

	matrix[N, J] X; //matrix of predictors
	matrix[N_pred, J] X_pred; //new matrix for predictions
	vector[K] y[N]; //matrix of response values

	vector<lower=-1,upper=1>[K] direction_manager; //direction to which the function must go (-1 minimized, 1 maximized)
	vector<lower=0,upper=10>[K] importance_manager; //importance weight of the function

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

	for(n in 1:N){
		mu[n] = Xbeta_temp[n] + z[F_id[n]] * rep_vector(sigma_fragment, K);
		//mu[n] = Xbeta_temp[n] + z[F_id[n]] * sigma_fragment;
	}

	//priors
	to_vector(beta) ~ normal(0,5);
	L_omega ~ lkj_corr_cholesky(2);
	L_sigma ~ cauchy(0,2.5);
	sigma_fragment ~ normal(0,1);

	//likelihood
	y ~ multi_normal_cholesky(mu,L_Sigma);
}

generated quantities {
	matrix[N, K] y_post;//for posterior predictive checks
	matrix[N_pred, K] y_pred; //predicted values on a new X matrix
	vector[N] log_lik;
	matrix[K, K] Sigma; //the variance-covariance matrix
	matrix[K, N_pred] new_mat1;

	vector[N_pred] desirability_manager; //the desirability scores for forest managers

	Sigma = tcrossprod(L_Sigma);
	for (n in 1:N){
		y_post[n] = to_row_vector(multi_normal_cholesky_rng(to_vector(X[n,] * beta),L_Sigma)); //note that the fragment-level deviation is not taken into account there
                log_lik[n] = multi_normal_lpdf(y[n] | X[n,] * beta, Sigma);
	}
	for(n in 1:N_pred){
		y_pred[n] = to_row_vector(X_pred[n,] * beta); //note that the fragment-level and the residual covariations deviation is not taken into account there
	}
	new_mat1 = revv(y_pred, direction_manager, N_pred, K);
	
	desirability_manager = desirability(new_mat1, importance_manager, N_pred, K);

}




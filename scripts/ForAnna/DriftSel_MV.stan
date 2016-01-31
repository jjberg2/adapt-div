data {
	int 								N_inds;
	int 								N_obs_data_pts;
	int 								N_unobs_data_pts;
	int 								N_pops;
	int 								N_traits;
	int								obs_data_idx[N_obs_data_pts] ;
	int 								unobs_data_idx[N_unobs_data_pts];
	int								ind_record[N_inds,6];
	vector[N_obs_data_pts] 			obs_trait_vals;
	cholesky_factor_corr[N_pops] 	chol_pop_corr;
	vector[N_pops]					rt_pop_drifts;
	cholesky_factor_corr[N_inds] 	chol_study_corr;
	vector[N_inds]					rt_ind_vars;
}
transformed data {
	int								tot_N_data;
	tot_N_data <- N_obs_data_pts + N_unobs_data_pts;
}
parameters {
	real								lnrm_mu;
	real	<lower=0>					lnrm_sig;
	vector[N_unobs_data_pts] 		x;
	cholesky_factor_corr[N_traits] 	chol_G_corr;
	vector<lower=0>[N_traits]		G_vars;
	vector[N_pops*N_traits]			A_raw; // raw pop level effects
	vector[tot_N_data]				S_raw; // raw individual level effects
	real	<lower=0>					V_e[N_traits]; // environmental variance
	real								mu[N_traits]; // global means
}
transformed parameters {
	// vector of both observed and unobserved data
	real																		comb_trait_vals[N_traits*N_inds]; 
	real																		tmp_mp_A;
	cholesky_factor_corr[N_pops*N_traits]		G_by_pop_chol_corr;
	cholesky_factor_corr[N_inds*N_traits]		G_by_study_chol_corr;
	vector[N_pops*N_traits]									A; // pop level effects
	vector[tot_N_data]											S; // individual level effects	
	vector[tot_N_data]											gen_vals; // vector of all genetic effects
	

	// fill in G %x% Pop matrix
	for ( i in 1:N_traits ) {
		for ( j in 1:N_traits) {
			for ( k in 1:N_pops ) {
				if ( i == N_traits ) {
					// fill in appropriately scaled A vector
					A [ (j-1)*N_pops+k ] <- A_raw [ (j-1)*N_pops+k ]*G_vars[j]*rt_pop_drifts[k];
				}
				for ( l in 1:N_pops) {
					//print ( i , j , k , l );
					G_by_pop_chol_corr [ (i-1)*N_pops+k , (j-1)*N_pops+l ] <- chol_G_corr [i,j]*chol_pop_corr[k,l];
				}
			}
		}
	}
	
	// fill in G %x% theta_b matrix
	for ( i in 1:N_traits ) {
		for ( j in 1:N_traits) {
			for ( k in 1:N_inds ) {
				if ( i == N_traits) {
					// fill in appropriately scaled S vector
					S [ (j-1)*N_inds+k ] <- S_raw [ (j-1)*N_inds+k ]*G_vars[j]*rt_ind_vars[k];
					//print ( (j-1)*N_inds+k );
					//print ( S [ (j-1)*N_inds+k ] );
				}
				for ( l in 1:N_inds) {
					//print ( i , j , k , l );
					G_by_study_chol_corr [ (i-1)*N_inds+k , (j-1)*N_inds+l ] <- chol_G_corr [i,j]*chol_study_corr[k,l];
				}
			}
		}
	}
	
	for ( j in 1:N_traits ) {
		for ( i in 1:N_inds ) { 
			tmp_mp_A <- A [ (j-1)*N_pops + ind_record[i,5]]/2 + A [ (j-1)*N_pops + ind_record[i,6]]/2;
			gen_vals [ (j-1)*N_inds + i ] <- mu[j] + tmp_mp_A + S [ (j-1)*N_inds + i ];
			//print ( A [ (j-1)*N_pops + ind_record[i,5]] );
			//print ( (j-1)*N_pops);
			//print ( (j-1)*N_traits + i );
			//print ( gen_vals[(j-1)*N_inds + i] );
		}
	}
	
	// add observed data to combined vector
	for ( j in  1:N_obs_data_pts ) {
		comb_trait_vals[obs_data_idx[j]] <- obs_trait_vals[j];
		//print ( comb_trait_vals[obs_data_idx[j]] )
	} 
	// add unobserved data to combined vector
	for ( j in  1:N_unobs_data_pts ) {
		comb_trait_vals[unobs_data_idx[j]] <- x[j];
		//print ( comb_trait_vals[obs_data_idx[j]] )
	}
}
model {
	vector[tot_N_data]		expand_V_e;
	for ( j in 1:N_traits ) {
		for ( i in 1:N_inds ) {
			expand_V_e[(j-1)*N_inds+i] <- V_e[j];
		}
	}
	chol_G_corr ~ lkj_corr_cholesky(2);
	A_raw ~ multi_normal_cholesky(rep_vector(0,N_pops*N_traits),G_by_pop_chol_corr);
	S_raw ~ multi_normal_cholesky(rep_vector(0,N_inds*N_traits),G_by_study_chol_corr);	
	G_vars ~ lognormal ( lnrm_mu , lnrm_sig );
	V_e ~ lognormal ( lnrm_mu , lnrm_sig );

	increment_log_prob(normal_log(comb_trait_vals,gen_vals,expand_V_e));
}
data {
	int 													N_inds;
	int 													N_obs_data_pts;
	int 													N_unobs_data_pts;
	int 													N_pops;
	int 													N_traits;
	int														ind_record[N_inds,6];
	int 													obs_pheno_record[N_obs_data_pts,2];
	int 													unobs_pheno_record[N_unobs_data_pts,2];
	vector[N_obs_data_pts] 									obs_trait_vals;
	corr_matrix[N_pops] 									pop_corr;
	vector[N_pops]											rt_pop_drifts;
	corr_matrix[N_inds] 									study_corr;
	vector[N_inds]											rt_ind_vars;
}
transformed data {
	int								tot_N_data;
	matrix[N_pops,N_pops]			inv_pop_corr;
	matrix[N_inds,N_inds]			inv_study_corr;
	real							det_pop_corr_term;
	real							det_study_corr_term;

	inv_pop_corr <- inverse(pop_corr);
	inv_study_corr <- inverse(study_corr);
	det_pop_corr_term <- N_pops*log_determinant(pop_corr)/2 + N_pops*N_traits*log(2*pi())/2;
	det_study_corr_term <- N_inds*log_determinant(study_corr)+ N_pops*N_traits*log(2*pi())/2;
	tot_N_data <- N_obs_data_pts + N_unobs_data_pts;
}
parameters {
	real									lnrm_mu;
	real	<lower=0>						lnrm_sig;
	vector[N_unobs_data_pts] 				x;
	cholesky_factor_corr[N_traits] 			chol_G_corr;
	vector<lower=0>[N_traits]				G_vars;
	matrix[N_pops,N_traits]					A_raw; // raw pop level effects
	matrix[N_inds,N_traits]					S_raw; // raw individual level effects
	real<lower=0>							V_e[N_traits]; // environmental variance
	real									mu[N_traits]; // global means
}
transformed parameters {
	// vector of both observed and unobserved data
	real										comb_trait_vals[N_inds,N_traits]; 
	real										tmp_mp_A;
	matrix[N_pops,N_traits]						A; // pop level effects
	matrix[N_inds,N_traits]						S; // individual level effects	
	matrix[N_inds,N_traits]						gen_vals; // vector of all genetic effects
	

	// fill in G %x% Pop matrix
	for ( j in 1:N_traits) {
		for ( k in 1:N_pops ) {
			A[k,j] <- A_raw [k,j]*G_vars[j]*rt_pop_drifts[k];
		}
	}
	
	// fill in G %x% theta_b matrix
	for ( j in 1:N_traits) {
		for ( k in 1:N_inds ) {
			S [k,j] <- S_raw [k,j]*G_vars[j]*rt_ind_vars[k];
			tmp_mp_A <- A [ind_record[k,5],j]/2 + A[ind_record[k,6],j]/2;
			gen_vals [k,j] <- mu[j] + tmp_mp_A + S [k,j];
		}
	}
	

	for ( k in 1:N_obs_data_pts) {
			comb_trait_vals[obs_pheno_record[k,1],obs_pheno_record[k,2]] <- obs_trait_vals[k];
	}

	for ( k in 1:N_unobs_data_pts) {
			comb_trait_vals[unobs_pheno_record[k,1],unobs_pheno_record[k,2]] <- x[k];
	}

}
model {
	real					pop_effects_log_lik;
	real					study_effects_log_lik;
	//vector[tot_N_data]		expand_V_e;
	//for ( j in 1:N_traits ) {
	//	for ( i in 1:N_inds ) {
	//		expand_V_e[(j-1)*N_inds+i] <- V_e[j];
	//	}
	//}
	chol_G_corr ~ lkj_corr_cholesky(2);
	

	pop_effects_log_lik <- -trace_gen_quad_form(chol_G_corr*chol_G_corr',inv_pop_corr,A_raw)/2;
	pop_effects_log_lik <- pop_effects_log_lik - det_pop_corr_term - N_traits*log_determinant(chol_G_corr);

	study_effects_log_lik <- -trace_gen_quad_form(chol_G_corr*chol_G_corr',inv_study_corr,S_raw)/2;
	study_effects_log_lik <- study_effects_log_lik - det_study_corr_term - N_traits*log_determinant(chol_G_corr);


	//A_raw ~ multi_normal_cholesky(rep_vector(0,N_pops*N_traits),G_by_pop_chol_corr);
	//S_raw ~ multi_normal_cholesky(rep_vector(0,N_inds*N_traits),G_by_study_chol_corr);	
	G_vars ~ lognormal ( lnrm_mu , lnrm_sig );
	V_e ~ lognormal ( lnrm_mu , lnrm_sig );
	increment_log_prob(pop_effects_log_lik);
	increment_log_prob(study_effects_log_lik);
	for ( i in 1:N_traits) {
		increment_log_prob(normal_log(comb_trait_vals[,i],gen_vals[,i],V_e[i]));
	}
}
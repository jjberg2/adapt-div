setwd("~/Documents/Academics/AdaptDiv/Scripts/ForAnna")
source ( "DataGen.R")
if ( FALSE ) test.data <- GenData ( 4 , 3 , 2, 5 , 2 )
library ( rstan )

PrepStanInput <- function ( data.matrix , pop.matrix , ind.var.matrix ) {
	
	#recover()
	NA.idx <- which ( is.na ( data.matrix ) , arr.ind = T )
	has.data.idx <- which ( !is.na ( data.matrix ) , arr.ind = T )
	
	
	#flat.has.data.idx <- which ( !is.na ( c ( data.matrix ) ) )
	obs.data <- data.matrix [ has.data.idx ]
	n.obs.data.pts <- nrow ( has.data.idx )
	n.unobs.data.pts <- nrow ( NA.idx )
	trait.groupings <- matrix ( seq_along ( data.matrix ) , nrow = nrow ( data.matrix ) )
	
	
	# 
	pop.corr <- cov2cor(pop.matrix)
	chol.pop.corr <- t ( chol ( pop.corr ) )
	rt.pop.drifts <- sqrt ( diag ( pop.matrix ) )
	
	ind.corr <- cov2cor ( ind.var.matrix )
	chol.ind.corr <- t ( chol ( ind.corr ) )
	rt.ind.vars <- sqrt ( diag ( ind.var.matrix ) )
	
	return ( list ( obs.data = obs.data , has.data.idx = has.data.idx , NA.idx = NA.idx , n.obs.data.pts = n.obs.data.pts , n.unobs.data.pts = n.unobs.data.pts , trait.groupings = trait.groupings , pop.corr = pop.corr , chol.pop.corr = chol.pop.corr , rt.pop.drifts = rt.pop.drifts , chol.ind.corr = chol.ind.corr , rt.ind.vars = rt.ind.vars ) )
	
}

code <- 
"
data {
	int 								N_inds;
	int 								N_obs_data_pts;
	int 								N_unobs_data_pts;
	int 								N_pops;
	int 								N_traits;
	int								obs_data_idx[N_obs_data_pts,2] ;
	int 								unobs_data_idx[N_unobs_data_pts,2];
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
	real										comb_trait_vals[N_traits,N_inds]; 
	cholesky_factor_corr[N_pops*N_traits]	G_by_pop_chol_corr;
	cholesky_factor_corr[N_inds*N_traits]	G_by_study_chol_corr;
	vector[N_pops*N_traits]					A; // pop level effects
	vector[tot_N_data]						S; // individual level effects	
	vector[tot_N_data]						gen_vals; // vector of all genetic effects
	
	for ( j in 1:N_traits ) {
		for ( i in 1:N_inds ) { 
			
			gen_vals [ (j-1)*N_traits + i ] <- (A [ (j-1)*N_pops + ind_record[i,5]] + A [ (j-1)*N_pops + ind_record[i,5]])/2 + S [ (j-1)*N_traits + i ] + mu[j];
		}
	}
	
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
					S [ (j-1)*N_pops+k ] <- S_raw [ (j-1)*N_inds+k ]*G_vars[j]*rt_ind_vars[k];
				}
				for ( l in 1:N_inds) {
					//print ( i , j , k , l );
					G_by_study_chol_corr [ (i-1)*N_inds+k , (j-1)*N_inds+l ] <- chol_G_corr [i,j]*chol_study_corr[k,l];
				}
			}
		}
	}
	
	
	// add observed data to combined vector
	for ( j in  1:N_obs_data_pts ) {
		comb_trait_vals[obs_data_idx[j,2],obs_data_idx[j,1]] <- obs_trait_vals[j];
	} 
	// add unobserved data to combined vector
	for ( j in  1:N_unobs_data_pts ) {
		comb_trait_vals[unobs_data_idx[j,2],unobs_data_idx[j,1]] <- x[j];
	}
}
model {
	chol_G_corr ~ lkj_corr_cholesky(10);
	G_vars ~ lognormal ( lnrm_mu , lnrm_sig );
	A_raw ~ multi_normal_cholesky(rep_vector(0,N_pops*N_traits),G_by_pop_chol_corr);
	S_raw ~ multi_normal_cholesky(rep_vector(0,N_inds*N_traits),G_by_study_chol_corr);
	
	V_e ~ lognormal ( lnrm_mu , lnrm_sig );
	for ( i in 1:N_traits) {
		increment_log_prob(normal_log(comb_trait_vals[i],0,V_e[i]));
	}
}
"

stan.input <- PrepStanInput ( test.data$obs.trait.vals , test.data$pop.matrix , test.data$theta.b.mat)

the.data <- list ( 	N_inds = nrow ( test.data$inds ) ,
					N_obs_data_pts = stan.input$n.obs.data.pts , 
					N_unobs_data_pts = stan.input$n.unobs.data.pts , 
					N_pops = test.data$n.pops , 
					N_traits = test.data$n.envs * test.data$n.trait , 
					obs_trait_vals = stan.input$obs.data , 
					obs_data_idx = stan.input$has.data.idx ,
					ind_record = test.data$inds ,
					unobs_data_idx = stan.input$NA.idx ,						
					chol_pop_corr = stan.input$chol.pop.corr ,
					rt_pop_drifts = stan.input$rt.pop.drifts ,
					chol_study_corr = stan.input$chol.ind.corr ,
					rt_ind_vars = stan.input$rt.ind.vars 
				)

my.model <- stan_model ( model_code = code )

blah <- sampling ( my.model , data = the.data , iter = 100 , chains = 2 , control = list ( refresh = 10 ) )
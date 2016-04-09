setwd("~/Documents/Academics/adapt-div/")
source ( "scripts/ForAnna/DataGen.R")
test.data <- GenData ( 4 , 3 , 3 , 10 , 5 )
require ( rstan )

PrepStanInput <- function ( data.matrix , pop.matrix , ind.var.matrix ) {
	
	#recover()
	NA.idx <- which ( is.na ( data.matrix ) , arr.ind = T )
	has.data.idx <- which ( !is.na ( data.matrix ) , arr.ind = T )
	
	
	flat.has.data.idx <- which ( !is.na ( c ( data.matrix ) ) )
	flat.NA.idx <- which ( is.na ( c ( data.matrix ) ) )
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
	
	return ( list ( obs.data = obs.data , has.data.idx = has.data.idx , NA.idx = NA.idx , n.obs.data.pts = n.obs.data.pts , n.unobs.data.pts = n.unobs.data.pts , trait.groupings = trait.groupings , pop.corr = pop.corr , chol.pop.corr = chol.pop.corr , rt.pop.drifts = rt.pop.drifts , chol.ind.corr = chol.ind.corr , rt.ind.vars = rt.ind.vars , flat.has.data.idx = flat.has.data.idx , flat.NA.idx = flat.NA.idx ) )
	
}

stan.input <- PrepStanInput ( test.data$obs.trait.vals , test.data$pop.matrix , test.data$theta.b.mat)

the.data <- list ( 	N_inds = nrow ( test.data$inds ) ,
					N_obs_data_pts = stan.input$n.obs.data.pts , 
					N_unobs_data_pts = stan.input$n.unobs.data.pts , 
					N_pops = test.data$n.pops , 
					N_traits = test.data$n.envs * test.data$n.trait , 
					obs_trait_vals = stan.input$obs.data , 
					obs_data_idx = stan.input$flat.has.data.idx ,
					unobs_data_idx = stan.input$flat.NA.idx ,		
					ind_record = test.data$inds ,				
					chol_pop_corr = stan.input$chol.pop.corr ,
					rt_pop_drifts = stan.input$rt.pop.drifts ,
					chol_study_corr = stan.input$chol.ind.corr ,
					rt_ind_vars = stan.input$rt.ind.vars 
				)


options(mc.cores = parallel::detectCores())
my.model <- stan_model ( file = "scripts/ForAnna/DriftSel_MV.stan" )

blah <- sampling ( my.model , data = the.data , iter = 2500 , chains = 4 , control = list ( max_treedepth = 12 , adapt_delta = 0.95 , stepsize = 0.01) , refresh = 1  , pars = c ( "A" , "S" , "chol_G_corr" , "G_vars" , "V_e" , "mu" , "lnrm_mu" , "lnrm_sig") , cores = 4 , open_progress = FALSE , diagnostic_file = "sims/diagnostics" , sample_file = "sims/samples")
save ( test.data , file = "sims/test.new.driftsel.simdata.Robj")
save ( blah , file = "sims/testing.new.driftsel.stanfit.Robj")
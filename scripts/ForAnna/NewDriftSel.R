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
	
	#obs.pheno.record <- cbind ( has.data.idx , obs.data )
	n.traits <- length ( unique ( has.data.idx [ , 2 ] ) )
	n.inds <- length ( unique ( has.data.idx [ , 1 ] ) )
	unobs.pheno.record <- list ()
	for ( i in 1 : n.traits ) {
		these <- has.data.idx [ , 2 ] == i
		these.missing <- which ( !(1:n.inds %in% has.data.idx [ these , 1 ]) )
		unobs.pheno.record [[ i ]] <- cbind ( these.missing ,  i )
	}
	unobs.pheno.record <- do.call ( rbind , unobs.pheno.record )
	unobs.pheno.record <- unobs.pheno.record [ order ( unobs.pheno.record [ , 2 ] , unobs.pheno.record [ , 1 ] ) , ] 

	# 
	pop.corr <- cov2cor(pop.matrix)
	chol.pop.corr <- t ( chol ( pop.corr ) )
	rt.pop.drifts <- sqrt ( diag ( pop.matrix ) )
	
	ind.corr <- cov2cor ( ind.var.matrix )
	chol.ind.corr <- t ( chol ( ind.corr ) )
	rt.ind.vars <- sqrt ( diag ( ind.var.matrix ) )
	
	return ( list ( obs.data = obs.data , has.data.idx = has.data.idx , NA.idx = NA.idx , n.obs.data.pts = n.obs.data.pts , n.unobs.data.pts = n.unobs.data.pts , trait.groupings = trait.groupings , pop.corr = pop.corr , chol.pop.corr = chol.pop.corr , rt.pop.drifts = rt.pop.drifts , ind.corr = ind.corr , chol.ind.corr = chol.ind.corr , rt.ind.vars = rt.ind.vars , flat.has.data.idx = flat.has.data.idx , flat.NA.idx = flat.NA.idx ,  unobs.pheno.record = unobs.pheno.record ) )
	
}

stan.input <- PrepStanInput ( test.data$obs.trait.vals , test.data$pop.matrix , test.data$theta.b.mat)

the.data <- list ( N_inds = nrow ( test.data$inds ) ,
					N_obs_data_pts = stan.input$n.obs.data.pts , 
					N_unobs_data_pts = stan.input$n.unobs.data.pts , 
					N_pops = test.data$n.pops , 
					N_traits = test.data$n.envs * test.data$n.trait , 
					obs_trait_vals = stan.input$obs.data , 		
					ind_record = test.data$inds ,
					obs_pheno_record = stan.input$has.data.idx ,
					unobs_pheno_record = stan.input$NA.idx ,			
					pop_corr = stan.input$pop.corr ,
					rt_pop_drifts = stan.input$rt.pop.drifts ,
					study_corr = stan.input$ind.corr ,
					rt_ind_vars = stan.input$rt.ind.vars 
				)


options(mc.cores = 1)
#options(mc.cores = parallel::detectCores())
my.model <- stan_model ( file = "scripts/ForAnna/DriftSel_MV.stan" )
my.matrixvar.model <- stan_model ( file = "scripts/ForAnna/DriftSel_MatVN.stan" )

blah <- sampling ( my.matrixvar.model , 
					data = the.data , 
					iter = 2500 , 
					chains = 2 , 
					control = list ( max_treedepth = 12 , 
										adapt_delta = 0.95
									) , 
					refresh = 1  , 
					pars = c ( "A" , 
								"S" , 
								"chol_G_corr" , 
								"G_vars" , 
								"V_e" , 
								"mu" , 
								"lnrm_mu" , 
								"lnrm_sig"
							) , 
					cores = 1 , 
					open_progress = FALSE 
				)
save ( test.data , file = "sims/test.new.driftsel.simdata.Robj")
save ( blah , file = "sims/testing.new.driftsel.stanfit.Robj")
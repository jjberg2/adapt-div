library ( mvtnorm )

RemoveMatches <- function ( x , matches ) 	 x [ !(x %in% matches) ]

ThetaBMatrix <- function ( ind1 , ind2 , mom1 , dad1 , mom2 , dad2 , pop.matrix , mom1.pop , dad1.pop ) {
	mom.drift <- pop.matrix [ mom1.pop , mom1.pop ]
	dad.drift <- pop.matrix [ dad1.pop , dad1.pop  ]
	theta.s.mom <- 1/2 + mom.drift/2
	theta.s.dad <- 1/2 + dad.drift/2
	if ( ind1 == ind2 ) return ( 1/2 - 1/4 * ( mom.drift + dad.drift ) )
	if ( mom1 == mom2 & dad1 == dad2 ) return ( 1/4 * ( theta.s.mom + theta.s.dad - mom.drift - dad.drift ) )
	if ( dad1 == dad2 & mom1 != mom2 ) return ( 1/4 * ( theta.s.dad - dad.drift ) )
	if ( dad1 != dad2 & mom1 == mom2 ) return ( 1/4 * ( theta.s.mom - mom.drift ) )
	
	return ( 0 )
	
}

GenData <- function ( n.pops , n.traits , n.envs , n.dad.per.pop , n.inds.per.env ) {

	pop.matrix <- rWishart ( 1 , 10 , diag ( n.pops ) )[ , , 1 ]/30
	g.matrix <- rWishart ( 1 , 30 , diag ( n.traits * n.envs ) )[,,1] / 30
	env.vars <- diag ( g.matrix ) * rnorm ( n.traits * n.envs , 0.4285 , sd = 0.03 )
	
	
	
	
	total.n.inds <- n.pops*n.dad.per.pop*n.envs*n.inds.per.env
	
	these.dads <- unlist ( sapply (1:n.pops , function ( x ) sample ( 0:99+x*1000 , size = n.dad.per.pop ), simplify = F ) )
	these.moms <- sample (  RemoveMatches ( 1000:(1999 + ( n.pops-1)*1000) , these.dads ) , total.n.inds )
	
	
	my.inds <- data.frame ( 
					ind.id = seq_len ( total.n.inds ) ,
					dads = rep ( these.dads , each = n.inds.per.env*n.envs ) ,
					moms = these.moms ,
					env = rep ( rep ( seq_len ( n.envs ) , each = n.inds.per.env ) , times = n.dad.per.pop*n.pops ) ,
					dad.pops = ( floor ( these.dads / 1000 ) ) ,
					mom.pops = ( floor ( these.moms / 1000 ) )
					)
	theta.b.mat <- matrix ( NA , nrow = nrow ( my.inds ) , ncol = nrow ( my.inds ) )
	for ( i in seq_len ( nrow ( theta.b.mat ) ) ) {
		for ( j in i:ncol ( theta.b.mat ) ) {
			
			theta.b.mat [ i , j ] <- 
			theta.b.mat [ j , i ] <- 
				ThetaBMatrix ( 
					my.inds$ind.id [ i ] , 
					my.inds$ind.id [ j ] ,
					my.inds$moms [ i ] ,
					my.inds$dads [ i ] ,
					my.inds$moms [ j ] ,
					my.inds$dads [ j ] ,
					pop.matrix ,
					my.inds$dad.pops [ i ] ,
					my.inds$mom.pops [ i ]
					)
		}
	}
	
	
	g.by.pop.matrix <- 2*g.matrix %x% pop.matrix
	g.by.ped.matrix <- 2*g.matrix %x% theta.b.mat
	
	tmp1 <- rmvnorm ( 1 , mean = rnorm ( n.traits * n.envs * n.pops ) , sigma = ( g.by.pop.matrix ) )
	pop.effects <- array ( tmp , dim = c ( n.pops , n.traits , n.envs  ) )
	tmp2 <- rmvnorm ( 1 , mean = rep ( 0 , total.n.inds*n.envs*n.traits ) , sigma = g.by.ped.matrix )
	ind.effects <- t ( matrix ( tmp2 , nrow = n.traits*n.envs ) )
	trait.ids <- paste ( rep ( sprintf ( "Trait%d" , 1:n.traits ) , n.envs ) , rep ( sprintf("_Env%d" , 1:n.envs ) , each = n.traits ) , sep = "" )
	colnames ( ind.effects ) <- paste ( trait.ids , "_IND" , sep = "" )
	
	pop.effects.by.ind <- matrix ( NA , nrow = nrow ( ind.effects ) , ncol = ncol ( ind.effects ) )
	colnames ( pop.effects.by.ind ) <- paste ( trait.ids , "_POP" , sep = "" )
	for ( i in seq_len ( total.n.inds ) ) {
		pop.effects.by.ind [ i , ] <- c ( pop.effects[ my.inds$dad.pops [ i ] , , ] + pop.effects [ my.inds$mom.pops [ i ] , , ] )
	}
	
	sum.genetic.effects <- pop.effects.by.ind + ind.effects
	ind.env.effects <- sapply ( env.vars , function ( x ) rnorm ( total.n.inds , 0 , sd = sqrt ( x ) ) )
	trait.vals <- sum.genetic.effects + ind.env.effects
	
	obs.trait.vals <- matrix ( NA , nrow = nrow ( trait.vals ) , ncol = ncol ( trait.vals ) )
	colnames ( trait.vals ) <- colnames ( obs.trait.vals ) <- trait.ids
	for ( i in seq_len ( total.n.inds ) ) {
		obs.trait.vals [ i , ( my.inds$env [ i ] - 1 ) * n.traits + seq_len(n.traits) ] <-  trait.vals [ i , ( my.inds$env [ i ] - 1 ) * n.traits + seq_len(n.traits) ]
	}
	return ( list ( inds = my.inds , ind.effects = ind.effects , pop.effects=pop.effects, pop.effects.by.ind=pop.effects.by.ind , ind.env.effects=ind.env.effects , trait.vals=trait.vals , obs.trait.vals=obs.trait.vals , pop.matrix=pop.matrix , g.matrix=g.matrix , theta.b.mat=theta.b.mat , env.vars=env.vars , n.pops = n.pops , n.traits = n.traits , n.envs = n.envs ) )
	
	
	

}



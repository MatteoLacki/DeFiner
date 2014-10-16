get.div.2 <- function( x ){
	W <- x %% 2

	return(  
		(0:((x-W)/2))*2 + W	
	)
}

get.get.lucky.configurations <- function( neutronsNo ) function(){
	X 	<- get.div.2( neutronsNo )
	XYZ	<- do.call(
		rbind,
		lapply( 
			X, 
			function( x ) data.frame( 
				x = x,
				y = get.div.2( (neutronsNo - x)/2)
			)
		)
	)

	XYZ$z<- apply(
		XYZ,1,
		function( w ) 	as.integer(
			( 
				(
					(neutronsNo - w[1])/2
				) - w[2] 
			)/2
		)
	)

	return( XYZ )
}


get.get.lucky.logProb 	<- function( intensities ){

	force(intensities)
	log.intensities 	<- log( unlist(intensities$lucky) )

	result 	<- function( configuration ) return(
			crossprod( configuration, log.intensities ) - 
			sum( lfactorial(configuration) )
		) 	

	return( result )
} 
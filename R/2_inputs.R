	# Probabilities and masses of elements

isotopeProbs <- list(
	C = c(98.93,   1.07)   /100,
	H = c(99.9885, 0.0115) /100,
	N = c(99.632,  0.368)  /100,
	O = c(99.757, 0.038, 0.205)/100,
	S = c(94.93, 0.76, 4.29, 0, 0.02)/100
)

elementsNames   <- names( isotopeProbs )

configurationNames  <- c(
	'c1', 'h1', 'n1', 'o1', 's1', 'o2', 's2', 's4'
)

isotopeMasses  <- list(
	C   = c( 12.0000000000, 13.0033548378),
	H   = c( 1.0078250321, 2.0141017780),
	N   = c( 14.00307400522, 15.0001088984),
	O   = c( 15.9949146, 16.9991312, 17.9991603),
	S   = c( 31.97207070, 32.97145843, 33.96786665, 0, 35.96708062)
)

get.get.basic.isoNo <- function( chemicalCompound ) lapply(
	chemicalCompound,
	function( element ){
		force( element )	
		function( x ) element - sum( x )
	} 
)


get.multi.is.assymptotic.nuissance <- function( chemicalCompound ){
	
	aC 	<<- as.integer(unlist(chemicalCompound))
	
	f 	<- function( configuration, ... ) 	configuration <= aC
	
	return( f )
} 

get.binom.is.assymptotic.nuissance <- function( chemicalCompound ){
	
	aC 	<<- as.integer(unlist(chemicalCompound))
	
	f 	<- function( configuration, ... ) 	configuration <= aC[4:5]
	
	return( f )
} 



get.get.mass<- function( chemicalCompound ){

	get.basic.isoNo 	<- lapply(
		chemicalCompound,
		function( element ){
			force( element )	
			function( x ) element - sum( x )
		}
	)

	names( get.basic.isoNo ) <- elementsNames

	get.mass 	<- function( configuration ){

		masses 	<- sapply( 
			elementsNames,
			function( elementName ){

				heavyPart 	<- configuration[[elementName]]	
				mass 		<- crossprod(
					isotopeMasses[[elementName]],
					c( 
						get.basic.isoNo[[elementName]](
							heavyPart
						),
						if( elementName == 'S' ) 
							c( heavyPart[1:2], 0L, heavyPart[3L] )	
						else heavyPart 
					)
				)

				return( mass )
			} 	
		)

		mass 	<- sum( masses )

		return( mass )
	}

	return( get.mass )
}

isoNames2           <- c('O', 'S')
distributionNames   <- c('multi', 'binom', 'lucky')



get.intensities <- function( chemicalCompound ){

	partial <- lapply(
		elementsNames,
		function( elementName )     isotopeProbs[[elementName]][-1]*chemicalCompound[[elementName]] 
	)

	names( partial ) <- elementsNames

	lucky   <- lapply( 
		c(1,2,4),
		function( noOfNeutrons )    sum( 
			sapply( partial, '[', noOfNeutrons ),
			na.rm   = TRUE 
		)
	)

	names( lucky )  <- c('mu', 'eta', 'gamma')

	result <- list(
		multi = sapply( partial, '[', 1)/lucky$mu,
		binom = sapply( partial[ isoNames2 ], '[', 2)/lucky$eta,
		lucky = lucky
	)

	return( result )
}


######################################################################

elementsNo  <- length( elementsNames )

Directions  <- expand.grid(1:elementsNo,1:elementsNo) 

Directions  <- do.call(
	cbind,
	apply(
		Directions,1,
		function( row ) if (row[1]!=row[2]) return(row)
	)
)

makeDirection   <- function( m ){ 
	res     <- integer( elementsNo )
	res[m]  <- c(1L,-1L)

	return( res )
}

Directions <- apply(
	Directions,2,
	makeDirection
)

rownames( Directions )  <- elementsNames

BinomialDirections 				<- matrix(c(1,-1,-1,1),2,2)
rownames( BinomialDirections ) 	<- isoNames2

Directions 	<- list(
	multi 	= Directions,
	binom 	= BinomialDirections
)
   
rm( BinomialDirections )
######################################################################
	# Input control.

makeCompound    <- function( C=0L, H=0L, N=0L, O=0L, S=0L, ... ){ 
	result          <- c(C,H,N,O,S)
	names(result)   <- elementsNames

	return( result )
}

unlistCompound  <- function( chemicalCompound ) do.call( makeCompound, chemicalCompound )



######################################################################
	# Calculation of the uncontrained probability.

get.real.log.pdf 	<- function( chemicalCompound ) {

	prob.functions  <- lapply(
		elementsNames, 
		function( elementName ){

			probs   <- isotopeProbs[[elementName]]  

			probs   <- probs[ nonNegProbs <- probs > 0 ]
			probsNo <- sum( nonNegProbs )

			maxNo   <- chemicalCompound[[elementName]] 

			result  <- if( probsNo == 2 )   
				function( configuration )
					dbinom(
						x       = configuration,
						size    = maxNo, 
						prob    = probs[2L],
						log     = TRUE
					) 
				else    
				function( configuration ) {

					basicIsoNo  <- maxNo - sum(configuration)

					if( basicIsoNo < 0 ) return( -Inf )
					else return(
						dmultinom(
							x       = c( 
								basicIsoNo, 
								configuration 
							),
							prob    = probs,
							log     = TRUE
						)
					)
				}
				
			return( result )
		}
	)

	names( prob.functions ) <- elementsNames

	real.log.pdf   <- function( configuration ) sum(
		sapply( 
			elementsNames,
			function( elementName )     return(
				prob.functions[[elementName]]( 
					configuration[[elementName]] 
				)
			)   
		)
	)
	

	return( real.log.pdf )
}   

######################################################################


conf2character <- function( configuration ){

	result <- sapply(
		configuration,
		paste,
		sep='',
		collapse='.'
	)

	result <- paste0(
		paste0( elementsNames, result),
		collapse = '_'
	)

	return( result )
}


c2char 	<- function( configuration ) paste0( 
	configuration, 
	collapse = '.' 
)

char2c 	<- function( charConfiguration )  	as.integer(
	strsplit(
		charConfiguration, 
		split="\\."
	)[[1]]
)



######################################################################


get.is.on.simplex <- function( neutronsNo )   function( configuration ){

	configuration   <- unlist( configuration)
	result          <- all( 
		configuration >= 0     & configuration <= neutronsNo 
	)   &   sum( configuration ) == neutronsNo

	return( result )
}   


######################################################################


split.up <- function( configuration ){

	strange <- configuration$S[3]
	names( strange )    <- 'S'

	list(
		multi    = sapply( 
			configuration, '[', 1
		),
		binom       = sapply( 
			configuration[c('O','S')], '[', 2 
		),
		lucky          = strange 
	)
}    


split.new.old   <- function( configuration, score ){

	splitting   <- rep.int( 'old', length( configuration ) )
	splitting[which.max( score )] <- 'new'

	return(
		split(
			configuration,
			splitting
		)
	)
}

  
merge.up  <- function(
	multi,
	binom,
	lucky
){
	result      <- as.list( multi )
	result$O[2] <- binom[1]
	result$S[2] <- binom[2]
	result$S[3] <- lucky

	return( result )
}




######################################################################
	library('hash')

storage <- hash()

is.stored   <- function( configuration )    has.key(
	conf2character( configuration ),
	storage
) 

store       <- function( configuration, prob )    storage[
	conf2character( configuration ) 
] <- prob


get.do.we.have<- function( hashtable ) function( configuration )
	has.key(  
		c2char( configuration ),
		hashtable
	)


get.store.locally	<- function( hashtable ) function( configurations )
{
	apply(
		configurations, 2,
		function( configuration ){ 
			probInd 	<- length(configuration)
			hashtable[[ c2char( configuration[-probInd ] ) ]] <- 
				 configuration[ probInd ]
			invisible()
		}
	)	

	invisible()
}



######################################################################

normalize <- function( logProb ){

	maxLogP 	<- max( logProb )

	logOdds 	<- logProb - maxLogP

	sumOdds 	<- sum(exp( logOdds ))

	return(
		exp( logOdds - log(sumOdds) )
	)
}


sum.logs <- function( logProb ){
	maxLogP  	<- max( logProb )

	logOdds 	<- logProb - maxLogP

	logSumOdds 	<- log(sum(exp( logOdds )))

	return( exp( maxLogP + logSumOdds ) )
}


######################################################################

	
makeAvergine <- function( n ){ 

    meanAtoms <- c( 
        C = 4.9384,  
        H = 7.7583, 
        O = 1.4773,
        N = 1.3577, 
        S = 0.0417 
    )

    res <- lapply(
        meanAtoms,
        function( meanAtom )    as.integer( round( n*meanAtom ))
    )

    names(res) <- elementsNames

    return( res )
}

makeAvergine(100)

	# That's ok: 6]
	# it's simply
	# C*P(^13 C) + H*P(^2 H) + N*P(^15 N) + O*P(^17 O) + 2*O*P(^18 O) +
	# S*(^33 S) + 2*S*P(^34 S) + 4*S*P(^36 S)
meanOneAtom <- sapply( 
	isotopeProbs,
	function( isotope ) crossprod( isotope, 0:(length(isotope)-1L) ) 
)


makeAvergines   <- function( averginesSize ){

	meanAtoms <- c( 
        C = 4.9384,  
        H = 7.7583, 
        O = 1.4773,
        N = 1.3577, 
        S = 0.0417 
    )

    avergines       <- lapply( 
        averginesSize, 
        function( n ) list( 
            aC      = (aC <- makeAvergine(n)),
            k       = as.integer(round( crossprod( meanOneAtom, unlist(aC) ) ) )  
        )
    )


    names(avergines) <- format( 
        averginesSize, 
        scientific = TRUE
    )


    return( avergines )
}
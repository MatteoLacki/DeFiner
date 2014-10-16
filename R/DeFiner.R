DeFiner <- function( chemicalCompound, neutronsNo ){

	chemicalCompound<- avergine[[1]]$aC 


	cores 			<- parallel::detectCores()


	neutronsNo 		<- avergine[[1]]$k

	get.mass 		<- get.get.mass( chemicalCompound )
	intensities 	<- get.intensities( chemicalCompound )

	# lucky.parameter <- .99
	lucky.parameter <- 1.99

	probabilityBuster <- get.probabilityBuster( chemicalCompound, intensities )


	log.pdf 		<- get.real.log.pdf( chemicalCompound )



		# Preparing the lucky distribution.
	get.lucky.configurations<- get.get.lucky.configurations( neutronsNo )
	luckyConfigurations 	<- get.lucky.configurations()

		# Getting rid of what is real only for Poisson: s_4 > S
	luckyConfigurations 	<- subset( luckyConfigurations, z <= chemicalCompound$S)


	get.lucky.logProb 		<- get.get.lucky.logProb( intensities )

	luckyConfigurations$logProb <- apply(
		luckyConfigurations, 1,
		get.lucky.logProb
	)

	luckyConfigurations 	<- luckyConfigurations[order(luckyConfigurations$logProb),]


	luckyConfigurations$prob<- normalize( luckyConfigurations$logProb )

	luckyProbs <- luckyConfigurations$prob
	luckyProbs <- luckyProbs[order(luckyProbs, decreasing= TRUE)]

	luckyCumulated  <- cumsum(luckyProbs)

		# Truncating to a given amount of probability.
		#
		# This might now be numerically stable.

	important 	<- luckyCumulated < lucky.parameter
	important 	<- rev(important)


	luckyConfigurations <- luckyConfigurations[important,]
	luckyConfigurations <- luckyConfigurations[nrow(luckyConfigurations):1,]

	# luckyConfiguration 	<- luckyConfigurations[1,]

	luckyConfigurations <- apply(
	 	luckyConfigurations,1,
	 	function( w ) list(w)
	)

	luckyConfigurations <- lapply(
		luckyConfigurations, 
		'[[',
		1
	)


	# luckyConfiguration <- luckyConfigurations[[1]]

	result 	<- parallel::mclapply(
		luckyConfigurations,
		function( luckyConfiguration ){


			xyz 	<- as.integer( luckyConfiguration[1:3] )

			configurations 	<- list(
				multi 	= probabilityBuster$multi( xyz[1] ),
				binom 	= probabilityBuster$binom( xyz[2] )
			)

			colnames( configurations$multi ) <- c(configurationNames[1:5], "probMulti")
			colnames( configurations$binom ) <- c(configurationNames[6:7], "probBinom")


			result 	<- apply(
				configurations$multi, 1,
				function( multi.conf ){
					result <- apply(
						configurations$binom, 1,
						function( binom.conf ){
	# multi.conf <- configurations$multi[1,]
	# binom.conf <- configurations$binom[1,]

							conf 	<- list(
								C 	= multi.conf[1],
								H 	= multi.conf[2],
								N 	= multi.conf[3],
								O 	= c( multi.conf[4], binom.conf[1] ),
								S 	= c( multi.conf[5], binom.conf[2], xyz[3] )
							)

								# Checking whether the configuration is ok.
								# feasible by the product of multinomials.


							logProb <- log.pdf( conf )

							mass 	<- get.mass( conf )

							return(c(
								unlist(conf),
								logProb = logProb,
								mass 	= mass
							))
						}
					)

					result <- t(result)

					return( list(result))
				}
			)

			result <- lapply( result, function( x ) x[[1]])
			result <- do.call( rbind, result)


			names( result ) <- configurationNames

			return( list(result) )
		},
		mc.cores = cores
	)




	result <- lapply( result, function( x ) x[[1]])
	result <- do.call( rbind, result)

	result <- as.data.frame( result )
	result <- subset( result, logProb != -Inf )

	result$prob <- normalize( result$logProb )

	ourProb 	<- sum.logs( result$logProb )

	brainProb   <- BRAIN::calculateIsotopicProbabilities(
	    chemicalCompound,
	    nrPeaks = neutronsNo + 1L
	)[neutronsNo + 1L]

	result 	<- list(
		coverage 		= exp( log(ourProb) - log( brainProb) ),
		configurations 	= result
	)

	return( result )
}

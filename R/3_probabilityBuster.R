require('hash')

get.probabilityBuster <- function( 
){

	ProbabilityBuster 	<- setRefClass(
		Class 	= "ProbabilityBuster",
		contains= "VIRTUAL",

		fields 	= list(
			percent 		= "numeric",
			storage 		= "hash",
			bfs.storage 	= "MatrixPriorityQueue",
			sum.prob 		= "numeric",
			directions 		= "matrix",
			aC 				= "integer"
		),

		methods = list(

			initialize 		= function( 
				maxNo, 
				percent 			= .99, 
				chemicalCompound 	= NULL, 
				...
			){
				if( !is.null(chemicalCompound) ){

					callSuper(...)

					percent 	<<- percent
					storage 	<<- hash()

					Mode 		<- 	get.mode( maxNo )
					logProb 	<- 	log.pdf( Mode ) 
					sum.prob 	<<-	exp( logProb )

						# Adding initial node to newly created bfs priority queue.
					bfs.storage <<- MatrixPriorityQueue$new( 
						Mode, 
						logProb
					)

						# Adding initial node to hash table.
					storage[ c2char( Mode ) ]<<- logProb

					aC 			<<- as.integer( unlist( chemicalCompound) )

					return(.self)
				}
			},


			get.mode 	= function( maxNo, ...){},	
			log.pdf 	= function( configuration, ...){},
			

			is.stored 	= function( 
				configurations, 
				...
			){
				'Checks if provided configurations were already visited.'

				result 	<- apply(
					configurations, 2, 
					function( configuration ) 	has.key(
						c2char( configuration ),
						storage
					)
				)

				return( result )
			}, 


			store 		= function( configurations, ... ){	
				'Stores configurations on hash table and priority queue.'

				apply(
					configurations, 2,
					function( configuration ){ 
						probInd 	<- length(configuration)

						storage[[ c2char( configuration[-probInd ] ) ]] <<- 
							configuration[ probInd ]

						invisible()
					}
				)	

				bfs.storage$enqueue( configurations )

				invisible()
			},



			is.isotope.no.ok = function( ... ){},


			get.more.configurations = function(...){
				'Generates neighbouring configurations, taking into account ones outside the simplex and the already visited ones.'

				usingMethods( 'is.isotope.no.ok' )
				usingMethods( 'is.stored' )
				usingMethods( 'log.pdf' )


					# Take the highest priority diophantine solution.
				currentConf <- bfs.storage$poll()

					# Find all neighbouring solutions.
				candidates 	<- currentConf + directions			

					# Check if they are still really solutions.
				onSimplex 	<- apply( 
					candidates, 2,
					function( col ) all( col >= 0L )		 
				)
					# ... there will always be at least one: the node's father.

				candidates 	<- candidates[, onSimplex, drop = FALSE ]

				visited 	<- is.stored( candidates )


				if( !all( visited ) ){
					candidates 	<- candidates[, !visited, drop = FALSE ]

					isotope.no.ok <- apply(
						candidates, 2,
						is.isotope.no.ok
					)

					if( any(isotope.no.ok)){

						candidates 	<- candidates[, isotope.no.ok, drop = FALSE ]


							# Count proxy log-probability.
						candidates 	<- rbind(
							candidates,
							logProb = logProb <- apply( 
								candidates, 2,
								log.pdf
							)		
						)


						store( candidates )

							# This might be quite tricky numerically.
						sum.prob 	<<- sum.prob + 	sum( exp( logProb ) )
					}
				} 

				invisible()
			},


			get.most.prob 	= function(...){
				'This function performs prioritised bfs to find at least the prespecified amount of probability.'

					# This is the most strange thing I ever ever seen.


				timing <- system.time({
					while( sum.prob < percent && bfs.storage$size() > 0 ) 	
						get.more.configurations()				
				})

				# print( timing )

				localStorage 	<- as.list( storage )
				localStorage 	<- cbind(
					t(
						sapply( 
							names( localStorage ),
							char2c,
							USE.NAMES 	= FALSE
						)
					),
					unlist( localStorage, use.names = FALSE )
				)

				return( localStorage )
			}
		) 	
	)


	multi 	<- setRefClass(
		Class 	= "Multi",
		contains= "ProbabilityBuster",

		methods = list(

			initialize 	= function( 
				maxNo, 
				percent 	= .99,
				chemicalCompound = NULL,
				... 
			){
				if( !is.null(chemicalCompound) ){
					callSuper( maxNo, percent, chemicalCompound, ... )

					directions 	<<- Directions$multi
				}	
			},


			is.isotope.no.ok = function( 
				configuration, 
				... 
			){
				res 	<- all( configuration <= aC ) 
				return( res )
			},


			get.mode 	= function( maxNo, ...){
				'Gets a node in the vicinity of the mode of the multinomial distribution.'

				dimens  <- length( 		intensities$multi )
				minInd  <- which.min( 	intensities$multi )
				others  <- setdiff( 	1:dimens, 	minInd )
				result  <- integer( 	dimens )

				result[ others ]    <- floor( 
					(maxNo + 1) * intensities$multi[ others ] 
				)

				othersSum           <- sum( result[ others ] )

				if( othersSum <= maxNo ){ 
					result[minInd]  <- maxNo - othersSum

				} else {

					otherRes 		<- result[ others ]
					cond 			<- otherRes > 0

					if( sum( cond) > 0 ){
						otherRes 		<- otherRes[ cond ]
						
						minInd2         <- which.min( otherRes )
						
						result[minInd2] <- result[minInd2] - 1L 
					} else stop('Something is wrong with get.mode function for multinomial distribution.')
				}

				return( result )
			}
		)
	)



	get.multinomial.configurations 	<- function( 
		maxNo,
		percent 			= .99,
		chemicalCompound,
		... 
	){	
		multinomial.buster <- multi$new( maxNo, percent, chemicalCompound ) 	

		return( multinomial.buster$get.most.prob() )
	}


	binom 	<- setRefClass(
		Class 	= "Binom",
		contains= "ProbabilityBuster",

		methods = list(

			initialize 	= function( 
				maxNo,
				percent = .99,
				chemicalCompound = NULL,
				... 
			){
				if( !is.null(chemicalCompound) ){
					callSuper( maxNo, percent, chemicalCompound, ... )

					directions 	<<- Directions$binom
				}
			},



			is.isotope.no.ok = function( 
				configuration, 
				... 
			){
				res 	<- all( configuration <= aC[4:5] )
				return( res )
			},


			get.mode 	= function( 
				maxNo, 
				...
			){
				'Gets a node in the vicinity of the mode of the binomial distribution.'

				res <- c(
				 	o2 	<- max(
				 		0L,
				 		floor( 
					 		(maxNo+1)*intensities$binom[1] 
					 	) - 1L 
				 	), 
				 	maxNo - o2	
				)

				names( res ) <- c('O', 'S') 

				return( res )
			}
		)
	)


	get.binomial.configurations 	<- function( 
		maxNo,
		percent 	= .99,
		chemicalCompound,
		... 
	){	
		binomial.buster <- binom$new( maxNo, percent, chemicalCompound ) 	

		return( binomial.buster$get.most.prob() )
	}

}


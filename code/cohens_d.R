cohens_d <- function(x, y, DIM=1, SIGN=TRUE, na.rm=TRUE) {
#
# Function will compute cohen's d effect size.
# Generalized to work on either matrices, data frames, vectors
#
# INPUT
#	x <- matrix or numeric vector
#	y <- matrix or numeric vector
#	DIM <- specify the dimension which samples are along
#
# Example usage
#
# x <- cbind(rnorm(100,0,1), rnorm(100,0,1), rnorm(100,0,1))
# y <- cbind(rnorm(100,1,1), rnorm(100,2,1), rnorm(100,3,1))
# d <- cohens_d(x, y, 1)
#
# written by mvlombardo - 28.08.2015
#

	library(matrixStats)
	
	# if x and y are vectors, coerce them into matrices
	if (class(x)=="numeric" | class(x)=="integer") {
		x <- as.matrix(x)
	} # if

	if (class(y)=="numeric" | class(y)=="integer") {
		y <- as.matrix(y)
	}# if
  
  if (na.rm==TRUE){
    missingValDecision = TRUE
  } else {
    missingValDecision = FALSE
  }
  
	# n-1 for x and y
	lx <- dim(x)[DIM]-1
	ly <- dim(y)[DIM]-1

	# if samples are along the rows
	if (DIM==1){
		if (SIGN){
			# mean difference (numerator)
			md <- colMeans(x, na.rm = missingValDecision) - colMeans(y, na.rm = missingValDecision)
		} else{
			# mean difference (numerator)
			md <- abs(colMeans(x, na.rm = missingValDecision) - colMeans(y, na.rm = missingValDecision))
		}# if (SIGN)
		# pooled variance (denominator), but before any sqrt is done
		csd <- (lx * rowVars(t(x),na.rm = missingValDecision)) + (ly * rowVars(t(y), na.rm = missingValDecision))

	# else if samples are along the columns
	} else if (DIM==2){
		if (SIGN){
			# mean difference (numerator)
			md <- rowMeans(x, na.rm = missingValDecision) - rowMeans(y, na.rm = missingValDecision)
		} else{
			# mean difference (numerator)
			md <- abs(rowMeans(x, na.rm = missingValDecision) - rowMeans(y, na.rm = missingValDecision))
		}# if (SIGN)
		# pooled variance (denominator), but before any sqrt is done
		csd <- lx * rowVars(x, na.rm = missingValDecision) + ly * rowVars(y, na.rm = missingValDecision)
	}# end if

	# divide pooled variance by sum of n-1 for x and y and then square root it
	csd <- sqrt(csd/(lx + ly))
	# compute cohen's d
	cd  <- md/csd
}# end cohens_d <- function(x, y, DIM)

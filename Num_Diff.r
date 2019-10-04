#############################
# A simple function for 
# numerical differentiation
#############################


# Inputs 
# f: function to differentiate
# a: lower bound of study interval 
# b: upper bound of study interval
# order: order of differentiation (default=1)
# n: grid size (default = 100)

# Output
# A list with components
#   x:     evaluation grid
#   deriv: numerical derivative
#   order: differentiation order
 
numDiff1 <- function(f,a,b,order=1,n=100)
{
	stopifnot(a < b)
	k <- as.integer(order)	
	x <- seq(a,b,len=n+k)
	h <- x[2] - x[1]
	fx <- sapply(x,f)
	Dkfx <- diff(fx, differences=k) / h^k
	return(list(x=tail(x,n), deriv=Dkfx, order=k))	
}


# >@@@@ How would you modify the function   @@@@< #
# >@@@@ to calculate backward differences?  @@@@< #



########
# Tests
########



## Test 1: first derivative of sine function
result1 <- numDiff1(sin,0,pi/2,order=1,n=100)
plot(result1$x, result1$deriv, type="l", xlab="x", ylab="f'(x)")
curve(cos, 0, pi/2, add=TRUE, lty=2)
legend("topright", c("True derivative","Numerical approximation"), 
	lty=2:1, cex=.9)


## Test 2: second derivative of logarithmic function
f <- log
D2f <- function(x) -1/x^2
result2 <- numDiff1(f,1,2,order=2,n=100)
plot(result2$x, result2$deriv, type="l", xlab="x", ylab="f''(x)")
curve(D2f, 1, 2, add=TRUE, lty=2)
legend("topleft", c("True derivative","Numerical approximation"), 
	lty=2:1, cex=.9)



###############################
# A more accurate function for 
# numerical differentiation
###############################


numDiff2 <- function(f,a,b,order=1,n=100)
{
	stopifnot(a < b)
	k <- as.integer(order)	
	if (k == 0) { # case: no derivative
		x <- seq(a,b,len=n)
		fx <- tryCatch(f(x), error = function(e) NULL)
		if (is.null(fx) || length(fx) != n)
			fx <- sapply(x,f)	
		Dkfx <- fx		
	} else { 
		odd <- (k %% 2 == 1)
		m <- ifelse(odd, k+2, k+1) # number of support points
		x <- seq(a,b,len=n+m-1) # evaluation grid
		h <- x[2] - x[1] # mesh size
		halfm <- (m-1)/2
		pos <- (-halfm):halfm
		lhs <- matrix(0,m,m)
		for (i in 1:m) lhs[i,] <- pos^(i-1)
		rhs <- numeric(m)
		rhs[k+1] <- factorial(k)/h^k
		coefs <- solve(lhs,rhs) # coefficients in derivative formula
		fx <- tryCatch(f(x), error = function(e) NULL)
		if (is.null(fx) || length(fx) != n)
			fx <- sapply(x,f)	
		x <- x[(halfm+1):(n+halfm)]
		Dkfx <- filter(fx,rev(coefs))
		Dkfx <- na.omit(Dkfx)		
	}

	return(list(x=x, deriv=Dkfx, order=k))	
	
}


# >@@@@ How would you modify the function to        @@@@< #
# >@@@@ calculate forward and backward differences  @@@@< #
# >@@@@ near the boundaries? That is, use endpoint  @@@@< #
# >@@@@ formula instead of midpoint at boundaries   @@@@< #





######################
# Redo previous tests 
# with new function
######################

## Test 1: first derivative of sine function
result1.new <- numDiff2(sin, 0, pi/2, order=1,n=100)
plot(result1.new$x, result1.new$deriv, type="l", xlab="x", ylab="Df(x)")
curve(cos, 0, pi/2, add=TRUE, lty=2, lwd=3)
legend("topright", c("True derivative","Numerical approximation"), 
	lty=2:1, lwd=c(3,1), cex=.9)


## Test 2: second derivative of logarithmic function
f <- log
D2f <- function(x) -1/x^2
result2.new <- numDiff2(f, 1, 2, order=2, n=100)
plot(result2.new$x, result2.new$deriv, type="l", xlab="x", ylab="D2f(x)")
curve(D2f, 1, 2, add=TRUE, lty=2, lwd=3)
legend("topleft", c("True derivative","Numerical approximation"), 
	lty=c(2,1), lwd=c(3,1), cex=.9)


# Compare accuracy of approximations 
# in terms of mean absolute error
tab <- matrix(,2,2)
dimnames(tab) <- list(method=c("forward","higher"), 
	test=c("sin","log"))
tab["forward","sin"] <- mean(abs(result1$deriv-cos(result1$x)))
tab["forward","log"] <- mean(abs(result2$deriv +1/(result2$x)^2))
tab["higher","sin"] <- mean(abs(result1.new$deriv-cos(result1.new$x)))
tab["higher","log"] <- mean(abs(result2.new$deriv +1/(result2.new$x)^2))
tab



# >@@@@ Try running numDiff1 and numDiff2 with             <@@@@ #
# >@@@@ different values of n and different test functions <@@@@ #
# >@@@@ Can you find an example for which numerical        <@@@@ #
# >@@@@ differentiation is very inaccurate?                <@@@@ #



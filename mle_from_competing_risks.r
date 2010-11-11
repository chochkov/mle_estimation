
#=============================================================================================
# MLE of three parameters from a Competing Risks Model :: nikola chochkov c 2010                                          
#---------------------------------------------------------------------------------------------

RNGkind("Mersenne-Twister")                              # to set the kind of RNG in use (not necessary)
set.seed(13)                                             # to specify seeds
library(splines)

# Parameter true values
assign('lambda0', 2E-05)
assign('mu0', 7.0)
assign('sigma0', 0.5)

# Sample size
assign('n', 1000)
assign('repetitions', 100)

repetitions.choice<-as.integer(readline("Number of repetitions, default 100: \n"))
if(repetitions.choice > 0 && !is.na(repetitions.choice)){repetitions <- repetitions.choice}

# results matrix
res.numerical <- res <- matrix(0, nrow = repetitions, ncol = 3)

# Hessian matrix
H <- matrix(0, nrow=3, ncol=3)

#=============================================================================================
# 1. Define the functions to compute Log-Likelihood, Gradient vector and Hessian matrix                                          
#---------------------------------------------------------------------------------------------

log_lik_fn = function(par){
	lambda <- par[1]
	mu <- par[2]
	sigma <- par[3]
	
	loglik <- Nf*log(lambda) - Nf*log(1 - 1/sqrt(2*pi)) - lambda*sum(X) - 1/(2*sigma^2)*sum((log(X) - mu)^2) + 1/(2*sqrt(2*pi)*sigma^2)*sum((log(Xf) - mu)^2) - Nsurv*log(sigma*sqrt(2*pi)) - sum(log(Xsurv))
	
	# we return the negative of the two parts because nlminb() is a minimizer.
	return(-loglik)
	}

log_lik_grad = function(par){
	lambda <- par[1]
	mu <- par[2]
	sigma <- par[3]
	
	d_lam <- Nf/lambda - sum(X)
	d_mu  <- 1/sigma^2 * sum(log(X) - mu) - 1/(sqrt(2*pi)*sigma^2) * sum(log(Xf) - mu)
	d_sig <- 1/sigma^3 * sum((log(X) - mu)^2) - 1/(sqrt(2*pi)*sigma^3) * sum((log(Xf) - mu)^2) - Nsurv/sigma
	
	return(c(-d_lam, -d_mu, -d_sig)) 
	}

log_lik_hassian = function(par){
	lambda <- par[1]
	mu <- par[2]
	sigma <- par[3]
	
	d_lam_lam <- -Nf/lambda	
	d_sig_sig <- 3/(sqrt(2*pi)*sigma^4) * sum((log(Xf) - mu)^2) -3/sigma^4 * sum((log(X) - mu)^2) + Nsurv/sigma^2
	d_mu_mu   <- Nf/(sqrt(2*pi)*sigma^2) - n/sigma^2
	d_mu_sig  <- 2/(sqrt(2*pi)*sigma^3) * sum(log(Xf) - mu) - 2/sigma^3 * sum(log(X) - mu)

	H[1,1] 	 		 <- -(d_lam_lam)
	H[2,2]	  		 <- -(d_mu_mu)
	H[3,3]   		 <- -(d_sig_sig)
	H[2,3] <- H[3,2] <- -(d_mu_sig)
	
    return(H)
	}

#=============================================================================================
# 2. Simulate the data-generating process for a Competing Risks Model
#---------------------------------------------------------------------------------------------

for(rep in 1:repetitions){
	
	# Data holder - first column will store the data, second column : 0 at non-failure, 1 at failure
	y <- matrix(0, nrow = n, ncol = 2)

	Lnorm <- rlnorm(n, meanlog = mu0, sdlog = sigma0 )       # Generate the Lognormally distributed data    
	Expon <- rexp(n, rate = lambda0)                         # and the Exponential random deviates

	# Apply the competing risks model:
	for(i in 1:n){

		# If the Exponential i-th value is smaller than the Lnorm, then we record a failure 
		if(Lnorm[i] > Expon[i]){
			y[i, 1] = Expon[i]
			y[i, 2] = 1
		}
		else{
			y[i, 1] = Lnorm[i]	
		}
	}

	# Number of failed and survived cases for this repetition.
	Nf = sum(y[,2])
	Nsurv = n - Nf

	# Data vectors for failed and survived cases in this repetition
	Xf = matrix(0, nrow=Nf, ncol=1)
	Xsurv = matrix(0, nrow=Nsurv, ncol=1)
	X = y[,1]
	
	cnt_f = cnt_s = 1
	for(i in 1:n){
		
		# the i-th item failed
		if(y[i,2] == 1){
			Xf[cnt_f] = y[i,1]
			cnt_f = cnt_f + 1
		}
		else{
			Xsurv[cnt_s] = y[i,1]
			cnt_s = cnt_s + 1
		}
	}
		
	mle <- nlminb(c(lambda0*.8, mu0*.8, sigma0*.8), obj=log_lik_fn, grad=log_lik_grad, hes=log_lik_hassian, lower=c(lambda0/2,mu0/2,sigma0/2), upper=c(lambda0*2,mu0*2,sigma0*2), control=list(eval.max = 300, iter.max = 300))

	lambda_hat = Nf/sum(X)
	mu_hat = (sum(log(Xsurv)) - 1/(2*sqrt(2*pi)))/Nsurv
	sigma_square_hat = sqrt((sum((log(X) - mu_hat)^2) - 1/(sqrt(2*pi)) * sum((log(Xf) - mu_hat)^2))/Nsurv)
	
	res.numerical[rep,] <- mle$par
	res[rep,] <- c(lambda_hat, mu_hat, sigma_square_hat) # mle$par

#	cat(c(rep, res[rep,], Nf, Nsurv, mle$message, log_lik_fn(res[rep,])))
#	cat("\n")
}

res.average = c(mean(res[,1]), mean(res[,2]), mean(res[,3]))
res.numerical.average = c(mean(res.numerical[,1]), mean(res.numerical[,2]), mean(res.numerical[,3]))

res.precision = (res.average - c(lambda0, mu0, sigma0))/c(lambda0, mu0, sigma0) * 100
res.numerical.precision = (res.numerical.average - c(lambda0, mu0, sigma0))/c(lambda0, mu0, sigma0) * 100

cat('precision of numerical optimization of Loglikelihood:\n')
cat(res.numerical.precision)
cat('\n')
cat('precision of the directly derived formulas:\n')
cat(res.precision)

#plot( ecdf(Expon), main="", xlab="" , fg = gray(0.7))
#lines( ecdf(Lnorm), main="", xlab="",  col=2)

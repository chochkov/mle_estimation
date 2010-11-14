
#=============================================================================================
# MLE of three parameters from a Competing Risks Model :: nikola chochkov c 2010                                          
#---------------------------------------------------------------------------------------------

# set the kind of RNG in use (not necessary)
RNGkind("Mersenne-Twister")                              

# specify seeds
set.seed(17)
# library(splines)

# histograms display parameters
graphics.off()
par(mfrow=c(4,3))

# Parameters true value matrix - we study 4 parameter sets
params <- matrix(c(6E-05,7.3,0.7, 4E-05,7.5,0.75, 2E-05,7.75,0.5, 1E-06,9,0.75), nrow = 4, ncol=3, byrow=TRUE,
               dimnames = list(c("Set 1", "Set 2", "Set 3", "Set 4"),
                               c("Lambda", "Mu", "Sigma")))

# Precision of the estimation for each parameter set
precision <- matrix(0, nrow = 4, ncol=3, byrow=TRUE,
               dimnames = list(c("Set 1", "Set 2", "Set 3", "Set 4"),
                               c("Bias Lambda %", "Bias Mu %", "Bias Sigma %")))

# Sample size /e.g. total units under test/
assign('n', 1000)
assign('repetitions', 500)

#repetitions.choice<-as.integer(readline("Number of repetitions, press 'enter' for default 100: \n"))
#if(repetitions.choice > 0 && !is.na(repetitions.choice)){repetitions <- repetitions.choice}

# results matrix
res.numerical <- res <- matrix(0, nrow = repetitions, ncol = 3)

cat('Matrix of true parameter values: \n')
print(params)

# Loop over each set of true values in params matrix
for (params.set in 1:4){

	assign('lambda0', params[params.set, 1])
	assign('mu0', params[params.set, 2])
	assign('sigma0', params[params.set, 3])

	#=============================================================================================
	# 1. Simulate the data-generating process for a Competing Risks Model
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
		
		#=============================================================================================
		# 2. Estimator formulas:
		#---------------------------------------------------------------------------------------------
		
		lambda_hat = Nf/sum(X)
		mu_hat = (sum(log(Xsurv)) - 1/(2*sqrt(2*pi)))/Nsurv
		sigma_hat = sqrt((sum((log(X) - mu_hat)^2) - 1/(sqrt(2*pi)) * sum((log(Xf) - mu_hat)^2))/Nsurv)
	
		res[rep,] <- c(lambda_hat, mu_hat, sigma_hat)

	}

	#=============================================================================================
	# 3. Save the results for this iteration over the true values sets matrix - param
	#---------------------------------------------------------------------------------------------

	hist(res[,1],freq=FALSE,main=c(paste("Set ", params.set, ": Lognorm:", Nsurv, ", Expon:", Nf, sep='')), xlab=paste("Lambda, true val.", lambda0), ylab="")
	abline(v=lambda0, col='blue')
	abline(v=mean(res[,1]), col='red', lty=2)
	xfit<-seq(min(res[,1]),max(res[,1]),length=40)
	yfit<-dnorm(xfit,mean=mean(res[,1]),sd=sd(res[,1]))
	lines(xfit, yfit, lwd=2, lty=3, col='green') 

	h<-hist(res[,2],freq=FALSE,main='',xlab=paste("Mu, true val.", mu0), ylab="")
	abline(v=mu0, col='blue')
	abline(v=mean(res[,2]), col='red', lty=2)
	xfit<-seq(min(res[,2]),max(res[,2]),length=40)
	yfit<-dnorm(xfit,mean=mean(res[,2]),sd=sd(res[,2]))
	lines(xfit, yfit, lwd=2, lty=3, col='green') 

	h<-hist(res[,3],freq=FALSE,main='',xlab=paste("Sigma, true val.", sigma0), ylab="")
	abline(v=sigma0, col='blue')
	abline(v=mean(res[,3]), col='red', lty=2)
	xfit<-seq(min(res[,3]),max(res[,3]),length=40)
	yfit<-dnorm(xfit,mean=mean(res[,3]),sd=sd(res[,3]))
	lines(xfit, yfit, lwd=2, lty=3, col='green') 

	#=============================================================================================
	# 4. Print the results table for this parameters set:
	#---------------------------------------------------------------------------------------------
	res.average = c(mean(res[,1]), mean(res[,2]), mean(res[,3]))
	precision[params.set,] = (res.average - c(lambda0, mu0, sigma0))/c(lambda0, mu0, sigma0) * 100

	res.list<- list(True.Value =params[params.set,], Avg.MLE =res.average, Bias.percentage =precision[params.set,], Lognormal.items =Nsurv, Exponential.items =Nf)
	res.frame<- data.frame(res.list)
	rownames(res.frame)<- c(paste("Set:", params.set, 'lam'), paste("Set:", params.set, 'mu'), paste("Set:", params.set, 'sigma'))

	print(res.frame)

}

par(mfrow=c(1,1))
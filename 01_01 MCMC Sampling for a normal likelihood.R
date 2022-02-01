#' ---
#' title: "MCMC Sampling for a normal likelihood"
#' author: "Tomas Gonzalez <structuredguy@gmail.com>"
#' date: "`r Sys.Date()`"
#' fig.align: "center"
#' bibliography: Statistics.bib
#' ---
#'
#' Given a sample assumed to be extracted from normal distribution, find the mean and variance by Gibbs sampling.
#' See @{congdon2014applied}

#+ unknown distribution
set.seed(1234L)
x <- rnorm (100, 10, 3); n <- 100; mn.x <- mean (x) ; var.x <- var(x)
cat ("Observed mean ", mn.x, "\n ")
cat ("Observed variance ", var.x, "\n ")


#' Gibbs Sampling in conjunction with Normal/Gamma Priors.
#' Parameters for Normal priors on mu, and gamma prior on tau.
#+ priors
m <- 0; V <- 100; h <- 1/V; alph <- 1; beta <- 1

#' MCMC sample settings, arrays for holding samples, initial
#' values (0 for mu, 1 for tau)
T <- 10000; B <- 1000; TB <- T-B
mu <- tau <- sigma2 <- numeric(T); mu[1] <- 0; tau[1] <- 1

# start Gibbs sampler loop
for (t in 2:T) {
  # full conditional for mu
  m1 <- (h* m + n* tau[t-1] * mn.x) / (h+n* tau[t-1])
  V1 <- 1/(h + n*tau[t-1])
  mu[t] <- rnorm(1, m1, sqrt(V1))
  # full conditional for tau
  alph1 <- alph + (n/2); beta1 <- (sum((x-mu[t])^2)/2) + beta
  tau[t] <- rgamma(1, alph1, beta1)
  sigma2[t] <- 1/tau [t]
}
# end loop

#' Note that for the particular sample of y-values considered here, we have mean(y)=
{{ mean(x) }}
#'  and var(y) =
{{ var(x) }}
#'.
#'  Numerical and graphical summaries may be obtained
#' (for iterations t > B) using the commands

# Retain samples after Burn-in
mu <- mu[B+1 : T] ; length(mu) <- TB;
sigma2 <- sigma2[B+1: T]; length(sigma2) <- TB
# parameter summaries
summary (mu); quantile (mu [1: TB], c (0.025, 0.05, 0.90, 0.975))
summary (sigma2); quantile (sigma2 [1: TB], c (0.025, 0.05, 0.90, 0.975))

#+ subplots
par (mfrow=c(2, 2))
# Trace Plots
plot (mu, type="l"); plot (sigma2, type="l")
# Marginal Posterior Densities
plot (density (mu), col="red",
      main=expression (paste ("Posterior Density of ", mu)))
plot (density (sigma2), col="blue",
      main=expression (paste ("Posterior Density of ", sigma^2)))

#+ Bibliography
#'

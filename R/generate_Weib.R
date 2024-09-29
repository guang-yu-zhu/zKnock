#' Function that simulates response from Cox model with Weibull baseline hazard.
#'
#' This function generates survival data using a Weibull baseline hazard. The event times are simulated
#' based on a linear predictor and mild censoring is applied to the data.
#'
#' @param N Sample size (number of subjects).
#' @param lambda0 Baseline hazard scale parameter (λ).
#' @param rho Baseline hazard shape parameter (ρ).
#' @param lp Linear predictor, typically a linear combination of covariates and their coefficients.

#' @return A survival object representing the simulated event times with mild censoring.
#' 
#' @examples
#' # Simulate 10 Gaussian covariate predictors:
#' X <- generate_X(n = 100, p = 10, p_b = 0, cov_type = "cov_equi", rho = 0.2)
#'
#' # Create linear predictor with first 5 beta-coefficients = 1 (all other zero)
#' lp <- generate_lp(X, p_nn = 5, a = 1)
#'
#' # Simulate from Weibull hazard with baseline hazard h0(t) = λ * ρ * t^(ρ-1)
#' # and linear predictor, whose first 3 coefficients are non-zero:
#' y <- generate_Weib(N = nrow(X), lambda0 = 0.01, rho = 1, lp = lp)
#' @family generate
#' @importFrom survival Surv
#' @export
generate_Weib <- function(N, lambda0, rho, lp) {

  # Censoring times ~ Exponential(lambdaC)
  lambdaC = 0.0005 # very mild censoring
  
  # Simulate Weibull latent event times
  v <- runif(n = N)
  Tlat <- (-log(v) / (lambda0 * exp(lp)))^(1 / rho)
  
  # Simulate censoring times
  C <- rexp(n = N, rate = lambdaC)

  # Calculate follow-up times and event indicators
  time <- pmin(Tlat, C) # follow-up times
  status <- as.numeric(Tlat <= C) # event indicators

  # Return survival object
  survival::Surv(time = time, event = status)
}
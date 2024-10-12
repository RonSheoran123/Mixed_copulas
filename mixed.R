library(quantmod)
library(copula)
library(PerformanceAnalytics)


start_time <- proc.time()


# Fetch stock data
getSymbols(c("GOOGL", "MSFT"), src = "yahoo", from = "2015-01-01", to = Sys.Date())
data <- merge(Cl(GOOGL), Cl(MSFT))
colnames(data) <- c("GOOGL", "MSFT")

# Calculate log returns
returns <- diff(log(data))[-1, ]
returns <- as.data.frame(returns)  # Convert to data frame for compatibility

# Transform returns to uniform margins
u <- pobs(returns)

# Function for Gaussian Copula
gaussian_copula <- function(u) {
  gaussian_cop <- normalCopula(dim = 2)
  fit_gaussian <- fitCopula(gaussian_cop, u, method = "ml")
  
  # Extract goodness of fit metrics
  return(list(copula = fit_gaussian@copula,
              logLik = logLik(fit_gaussian),
              AIC = AIC(fit_gaussian),
              BIC = BIC(fit_gaussian)))
}

# Function for Gumbel Copula
gumbel_copula <- function(u) {
  gumbel_cop <- gumbelCopula()
  fit_gumbel <- fitCopula(gumbel_cop, u, method = "ml")
  
  # Extract goodness of fit metrics
  return(list(copula = fit_gumbel@copula,
              logLik = logLik(fit_gumbel),
              AIC = AIC(fit_gumbel),
              BIC = BIC(fit_gumbel)))
}

clayton_copula <- function(u) {
  clayton_cop <- claytonCopula()
  fit_gumbel <- fitCopula(clayton_cop, u, method = "ml")
  
  # Extract goodness of fit metrics
  return(list(copula = fit_gumbel@copula,
              logLik = logLik(fit_gumbel),
              AIC = AIC(fit_gumbel),
              BIC = BIC(fit_gumbel)))
}

# Function for Mixed Copula
mixed_copula <- function(u) {
  # Fit Gaussian Copula
  gaussian_fit <- gaussian_copula(u)
  
  # Fit Gumbel Copula
  gumbel_fit <- gumbel_copula(u)
  
  # Fit Clayton Copula
  clayton_fit <- clayton_copula(u)
  
  # Combine using a weighted average of the two copulas
  #weights <- c(0.0, 0.0, 0.0) # Example weights; adjust as necessary.
  
  max_logLik <- -1e9
  min_AIC <- 1e9
  min_BIC <- 1e9
  
  w1 <- 0.00000001
  w2 <- 0.00000001
  #count <- 0
  while (w1 < 1) {
    while (w2 < 1) {
      if (w1 + w2 < 1) {
        
        w3 <- 1 - w1 - w2
        
        #count <- count + 1
        
        # Create a mixture using the fitted copulas
        combined_density <- function(x, y) {
          density_gaussian <- dCopula(cbind(x, y), gaussian_fit$copula)
          density_gumbel <- dCopula(cbind(x, y), gumbel_fit$copula)
          density_clayton <- dCopula(cbind(x, y), clayton_fit$copula)
          return(w1 * density_gaussian + w2 * density_gumbel + w3 * density_clayton)
        }
        
        # Log-likelihood and AIC/BIC are calculated based on individual fits
        total_logLik <- gaussian_fit$logLik + gumbel_fit$logLik + clayton_fit$logLik
        total_AIC <- gaussian_fit$AIC + gumbel_fit$AIC + clayton_fit$AIC
        total_BIC <- gaussian_fit$BIC + gumbel_fit$BIC + clayton_fit$BIC
        
        if(total_AIC < min_AIC){
          w1_f <- w1
          w2_f <- w2
          w3_f <- w3
          min_AIC <- total_AIC
          min_BIC <- total_BIC
          max_logLik <- total_logLik
        }
      }
      w2 <- w2 + 0.00000001
    }
    w1 <- w1 + 0.00000001
  }
  
  return(list(density = combined_density,
              logLik = max_logLik,
              AIC = min_AIC,
              BIC = min_BIC, weights = c(w1_f, w2_f, w3_f)))
}

# Run the copulas on the transformed returns data
gaussian_fit <- gaussian_copula(u)
gumbel_fit <- gumbel_copula(u)
clayton_fit <- clayton_copula(u)
mixed_fit <- mixed_copula(u)

# Combine results into a data frame for comparison
fit_results <- data.frame(
  Model = c("Gaussian Copula", "Gumbel Copula", "Clayton Copula", "Mixed Copula"),
  LogLikelihood = c(gaussian_fit$logLik, gumbel_fit$logLik, clayton_fit$logLik, mixed_fit$logLik),
  AIC = c(gaussian_fit$AIC, gumbel_fit$AIC, clayton_fit$AIC, mixed_fit$AIC),
  BIC = c(gaussian_fit$BIC, gumbel_fit$BIC, clayton_fit$BIC, mixed_fit$BIC),
  weights = c(mixed_fit$weights, 0)
)



end_time <- proc.time()

time_taken <- end_time - start_time


print(fit_results)
print("Time taken:")
print(time_taken)

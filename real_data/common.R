find_best_survreg_model <- function(time, status, covariates, distributions = c("weibull", "exponential", "gaussian", "logistic", "lognormal", "loglogistic")) {
  # Charger le package nécessaire
  if (!requireNamespace("survival", quietly = TRUE)) stop("Le package 'survival' est requis.")
  
  # # Créer la formule
  # formula <- as.formula(paste("Surv(", time, ",", status, ") ~ ", paste(covariates, collapse = " + ")))
  # 
  # Initialisation
  aic_values <- numeric(length(distributions))
  valid_distributions <- rep(FALSE, length(distributions))
  
  # Boucle sur les distributions
  for (i in seq_along(distributions)) {
    # Ajustement du modèle pour la distribution courante
    
        survival::survreg(Surv(time, status)~covariates, dist = distributions[i])
        
    
    # Vérification et collecte de l'AIC
    if (!is.null(fit)) {
      valid_distributions[i] <- TRUE
      aic_values[i] <- -2 * as.numeric(logLik(fit)) + 2 * length(coef(fit))
    } else {
      aic_values[i] <- Inf
    }
  }
  
  # Vérification des distributions valides
  if (!any(valid_distributions)) {
    stop("Aucune distribution valide n'a pu être ajustée.")
  }
  
  # Identifier la meilleure distribution
  names(aic_values) <- distributions
  best_distribution <- names(which.min(aic_values[valid_distributions]))
  
  # Message sur la meilleure distribution
  cat("La distribution avec le meilleur AIC est :", best_distribution, 
      "avec un AIC de", min(aic_values[valid_distributions]), "\n")
  
  # Ajuster le modèle avec la meilleure distribution
  best_model <- survival::survreg(Surv(time, status), dist = best_distribution)
  
  # Retourner le résultat
  list(
    best_model = best_model,
    best_distribution = best_distribution,
    aic_values = aic_values
  )
}

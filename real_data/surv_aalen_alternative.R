##' @title Fitting Latent Factor Mixed Models (Least squares algorithm)
##' @description Latent Factor Mixed Models (LFMMs) are factor regression models. 
##' The lfmm2 function estimates latent factors based on an exact least-squares
##'  approach.
##' @param input Continuous intermediary variables matrix  encompassing potential mediators with n rows and p columns.
##' @param env An explanatory variable matrix with n rows and d columns.
##' @param K latent factor number
##' @param lambda ridge penalization parameter
##' @param effect.sizes true or false to obtain effect sizes
##' @return an object with the following attributes
##'  U
##'  V
##' @export 
##' @author Florence Pittion, Magali Richard, Olivier Francois
##' @examples
##' data(simu_data)
##' K = 5
##' mod.lfmm1 = lfmm2_med(input = simu_data$M1, 
##' env = simu_data$X_binary, 
##' K = K,
##' effect.sizes = FALSE)

lfmm2_med = function(input,
                     env, 
                     K, 
                     lambda = 1e-5,
                     effect.sizes = FALSE){
  
  ## Check input response matrix 
  ## LEA  
  if (is.character(input)){
    Y <- read.lfmm(input)
    lst.unique <- unique(as.numeric(Y))
    if (9 %in% lst.unique){
      stop("'input' file contains missing data (9's). Use the 'impute()' function to impute them.")
    }
    if (-9 %in% lst.unique){
      stop("'input' file contains missing data (-9's). Use the 'impute()' function to impute them.")
    }
  } else {
    ## Y is an R object       
    if (is.null(input)){
      stop("NULL value for argument 'input'.")
    }
    Y <- as.matrix(input)
    Y[Y == 9] <- NA
    Y[Y == -9] <- NA
    if (anyNA(Y)) {
      stop("The input matrix contains missing values: NA, 9 or -9 not allowed.")
    }
  }
  
  ## Check independent/covariate env matrix  
  ## LEA 
  if (is.character(env)){
    X <- read.env(env)
    if (anyNA(X)){
      stop("'env' file contains missing data (NA).")
    }
  } else {
    if (is.null(env)){
      stop("NULL value for argument 'env'.")
    }
    X <- as.matrix(env)
    if (anyNA(X)) {
      stop("The environmental matrix contains NA.")
    }
  }
  
  if (length(K) > 1){
    stop("Multiple values of K not allowed.")
  }
  if (lambda <= 0){
    stop("The ridge regularization parameter must be positive.")
  }
  
  d <-  ncol(X) #number of environmental variables
  n <-  nrow(X) #number of individuals
  
  if (nrow(Y) != n){
    stop("Number of rows in the input matrix not equal to the number of rows in the 'env' matrix")    
  }
  
  if (n < d) {
    stop("The environmental covariate matrix X contains more columns (d) than rows (n).")
  }
  
  # centering  
  Xs <- scale(X, scale = FALSE)
  Ys <- scale(Y, scale = FALSE)
  
  # run SVD of X: X = Q Sigma R
  svx <- svd(x = Xs, nu = n)
  Q <- svx$u
  
  d_lambda <- c(sqrt(lambda/(lambda + svx$d)), rep(1, n-d))
  d_lambda_inv <- c(sqrt((lambda + svx$d)/lambda), rep(1, n-d))
  D_inv <- diag(d_lambda_inv)
  D  <- diag(d_lambda)
  
  # run SVD of modified Y    
  svk <- svd(D %*% t(Q) %*% Ys, nu = K)
  
  if (K > 1) {
    Sigma_k <- diag(svk$d[1:K])
  } else {
    Sigma_k <- as.matrix(svk$d[1])
  }
  
  # compute the latent matrix W
  W <- Q %*% D_inv %*% tcrossprod(svk$u %*% Sigma_k, svk$v[,1:K])
  
  # compute LFMM factors U and loadings V
  # Non orthogonal factors
  U <- crossprod(t(Q %*% D_inv), svk$u %*% Sigma_k)
  V <- svk$v[,1:K]
  
  # compute environmental effect sizes 
  if (effect.sizes){
    B <- (t(Ys - W) %*% Xs) %*% solve(t(Xs) %*% Xs + diag(lambda, nrow = d, ncol = d))
    B <- as.matrix(B)
  } else
  {B <-  matrix(NA)}
  
  obj= list()
  obj$U <- as.matrix(U)
  obj$V <- as.matrix(V)
  
  ## LEA 
  class(obj) = "lfmm2"
  return(obj)
}

#---------------------------------------------
##' @title P-values adjusted for latent factors computed by lfmm2
##' @description The function returns a vector of p-values for association 
##' between potential mediators and exposure/outcome variables adjusted for latent factors 
##' computed by \code{lfmm2}. As input, it takes an object of class
##' \code{lfmm2Class} with the data that were used to adjust the LFMM. 
##' If \code{full} is set to \code{FALSE}, the function computes significance 
##' values (p-values) for each exposure variable, otherwise it returns 
##' p-values for the full set of exposure variables.
##' @param object lfmm2Class object
##' @param input a response variable matrix with n rows and p columns
##' @param env An explanatory variable matrix with n rows and d columns.
##' @param covar covariables
##' @param genomic.control correct pvalue with genomic inflation factor
##' @param linear true or false (else is logistic)
##' @param family of logistic reg
##' @param full compute partial regression FALSE/TRUE
##' @return an object with the following attributes 
##' P-values, fscores, zscores, adj.squared, gif
##' @importFrom stats binomial glm lm median pchisq pf prcomp qchisq qf
##' @importFrom utils read.table
##' @export
##' @author Florence Pittion, Magali Richard, Olivier Francois
##' @examples 
##' attach(simu_data)
##' K = 5
##' mod.lfmm1 = lfmm2_med(input = simu_data$M1, 
##' env = simu_data$X_binary, 
##' K = K,
##' effect.sizes = FALSE)
##' 
##' res_reg1 = lfmm2_med_test(mod.lfmm1, 
##' input = simu_data$M1, 
##' env = simu_data$X_binary,
##' covar = cbind(simu_data$age, simu_data$gender),
##' genomic.control = TRUE)
##' 

lfmm2_med_test= function(object, 
                         input,
                         env,
                         covar,
                         full=FALSE,
                         genomic.control=TRUE, 
                         linear=TRUE,
                         family=binomial(link = "logit"))
{
  ## check object
  if (class(object)!="lfmm2"){stop("the object is not lfmm2 type")}
  ## Check input matrix   
  ## LEA  
  if (is.character(input)){
    warning("Reading large input files with 'read.lfmm()' may be slow. See 'data.table::fread()' for fast import.")
    Y <- read.lfmm(input)
    lst.unique <- unique(as.numeric(Y))
    if (9 %in% lst.unique){
      stop("'input' file contains missing data (9's). Use the 'impute()' function to impute them.")
    }
    if (-9 %in% lst.unique){
      stop("'input' file contains missing data (-9's). Use the 'impute()' function to impute them.")
    }
  } else {
    ## Y is an R object       
    if (is.null(input)){
      stop("NULL value for argument 'input'.")
    }
    Y <- as.matrix(input)
    Y[Y == 9] <- NA
    Y[Y == -9] <- NA
    if (anyNA(Y)) {
      stop("The input matrix contains missing values (NA or 9).")
    }
  }
  
  ## Check independent/covariate matrix  
  ## LEA 
  if (is.character(env)){
    X <- read.env(env)
    if (anyNA(X)){
      stop("'env' file contains missing data (NA).")
    }
  } else {
    if (is.null(env)){
      stop("NULL value for argument 'env'.")
    }
    X <- as.matrix(env)
    if (anyNA(X)) {
      stop("The environmental matrix contains NA.")
    }
  }
  
  d <-  ncol(X) #number of environmental variables
  n <-  nrow(X) #number of individuals
  
  if (nrow(Y) != n){
    stop("Number of rows in the input matrix not equal to the number of rows in the 'env' matrix")    
  }
  
  if (n < d) {
    stop("The environmental covariate matrix X contains more columns (d) than rows (n).")
  }
  
  p <- ncol(Y)
  gif <-  NULL
  p_value <- NULL
  z_score <- NULL
  f_score <- NULL           
  r_squared <- NULL
  
  if (full){
    ## a single p-value is returned (f-test)
    ## Check linear models  
    if (linear == FALSE){
      stop("Option full == TRUE is available only for linear models.")
    }
    
    if (is.null(covar)){
      ## partial regression
      mod_Y = lm(Y ~ ., data = data.frame(object$U)) 
      res_Y = mod_Y$residuals
      mod_X = lm(X ~ ., data = data.frame(object$U))
      res_X = mod_X$residuals
      
      mod_lm =  lm(res_Y ~ res_X)
      sm = summary(mod_lm)
      r_squared <- sapply(sm, FUN = function(x) x$adj.r.squared)
      f_score <- sapply(sm, FUN = function(x) x$fstat[1])
      p_value <- sapply(sm, FUN = function(x) pf(x$fstat[1], x$fstat[2], x$fstat[3], lower.tail = F))
    } else {
      mod_Y = lm(Y ~ ., data = data.frame(covar, object$U)) 
      res_Y = mod_Y$residuals
      mod_X = lm(X ~ ., data = data.frame(covar, object$U))
      res_X = mod_X$residuals
      
      mod_lm =  lm(res_Y ~ res_X)
      sm = summary(mod_lm)
      r_squared <- sapply(sm, FUN = function(x) x$adj.r.squared)
      f_score <- sapply(sm, FUN = function(x) x$fstat[1])
      p_value <- sapply(sm, FUN = function(x) pf(x$fstat[1], x$fstat[2], x$fstat[3], lower.tail = F))
    }
  } else {
    ## All p-values returned    
    if (linear){
      if(is.null(covar)){
        mod_lm <- lm(Y ~ ., data = data.frame(X, object$U)) 
        sm <- summary(mod_lm)
        p_value <- sapply(sm, FUN = function(x) x$coeff[2:(d + 1), 4])
        z_score <- as.matrix(sapply(sm, FUN = function(x) x$coeff[2:(d + 1), 3]))
      } else {
        mod_lm <- lm(Y ~ ., data = data.frame(X, covar, object$U)) 
        sm <- summary(mod_lm)
        p_value <- sapply(sm, FUN = function(x) x$coeff[2:(d + 1), 4])
        z_score <- as.matrix(sapply(sm, FUN = function(x) x$coeff[2:(d + 1), 3]))
      }
    } else {
      if(is.null(covar)){
        for (j in 1:p) {
          mod_glm <- glm(Y[, j] ~ ., data = data.frame(X, object$U), family = family)
          sm <- summary(mod_glm)
          p_value <- rbind(p_value, sm$coeff[2:(d + 1), 4])
          z_score <- rbind(z_score, sm$coeff[2:(d + 1), 3])
        }
      } else {
        for (j in 1:p) {
          mod_glm <- glm(Y[, j] ~ ., data = data.frame(X, covar, object$U), family = family)
          sm <- summary(mod_glm)
          p_value <- rbind(p_value, sm$coeff[2:(d + 1), 4])
          z_score <- rbind(z_score, sm$coeff[2:(d + 1), 3])
        }
      }
    }
  }
  
  if (genomic.control){
    if (!full){
      if (d == 1){
        gif <- median(z_score^2)/qchisq(0.5, df = 1, lower.tail = FALSE)
      } else {
        gif <- apply(z_score^2, 1, median)/qchisq(0.5, df = 1, lower.tail = FALSE)
      }
      p_value <- pchisq(z_score^2/gif, df = 1, lower.tail = FALSE)
    } else {
      gif <- median(f_score)/qf(0.5, d, n - d - 1, lower.tail = FALSE)
      p_value <- pf(f_score/gif, d, n - d - 1, lower.tail = FALSE)
    }
  }
  
  if (!full & anyNA(z_score)) {
    warning("NA in significance values and z-scores. Check the input matrix for non-variable genomic sites.")
  }
  if (full & anyNA(r_squared)) {
    warning("NA in significance values and R-squared. Check the input matrix for non-variable genomic sites.")
  }
  
  res <- list(pvalues = p_value, zscores = z_score, fscores = f_score, adj.r.squared = r_squared,  gif = gif)
  return(res)
}

read.env <- function(input.file) {
  # test arguments
  if(missing(input.file))
    stop("'input.file' argument is missing.")
  else if (!is.character(input.file))
    stop("'input.file' argument has to be of type character.")
  # check extension 
  test_extension(input.file, "env")
  return(as.matrix(read.table(input.file)));
}

read.lfmm <- function(input.file) {
  # test arguments
  if(missing(input.file))
    stop("'input.file' argument is missing.")
  else if (!is.character(input.file))
    stop("'input.file' argument has to be of type character.")
  # check extension 
  test_extension(input.file, "lfmm")
  return(as.matrix(read.table(input.file)))
}

test_extension <- function(name, extension)
{
  # obtain the extension of name
  ext = getExtension(basename(name))
  
  # if not the correct extension, stop
  if (ext != extension) {
    p = paste("'input_file' format and extension have to be \".", 
              extension, "\" (not \".",ext,"\").", sep="")
    stop(p)
  } 
  return(ext);
}

getExtension <- function(file) {
  l = strsplit(file, "\\.")[[1]]
  return(l[length(l)])
}

#########################################################################
######### TODO roxygen comments
surv_cox_med_test = function(object,
                             M,
                             survival_time,
                             censoring_status,
                             exposure,
                             covar=NULL) {
  Y = survival::Surv(survival_time, censoring_status)
  p = ncol(M)
  pvalues = zscores = c()
  for (j in 1:p) {
    if(is.null(covar)){
      cox_model = survival::coxph(Y ~ M[, j] + exposure + object$U)
      pvalues = c(pvalues, summary(cox_model)$coefficients[1, 5])
      zscores = c(zscores, summary(cox_model)$coefficients[1, 4])
    } else {
      cox_model = survival::coxph(Y ~ M[, j] + exposure + object$U + covar)
      pvalues = c(pvalues, summary(cox_model)$coefficients[1, 5])
      zscores = c(zscores, summary(cox_model)$coefficients[1, 4])
    }
  }
  res = list(pvalues = pvalues, zscores = zscores)
}

#########################################################################
######### TODO roxygen comments
surv_aalen_med_test = function(object,
                               M,
                               survival_time,
                               censoring_status,
                               exposure,
                               covar=NULL) {
  base = survival::Surv(survival_time, censoring_status)
  p = ncol(M)
  pvalues = c()
  
  for (j in 1:p) { 
    if (j%%100==0) {print(paste0("adjusting feature ",j,"/",p))}
    if(is.null(covar)){
      
      param_model = timereg::aalen(base ~ M[,j] + exposure + object$U)
      pvalues = c(pvalues, param_model$pval.testBeq0[2])
    } else {
      param_model = timereg::aalen(base ~ M[,j] + exposure + object$U + covar)
      pvalues = c(pvalues, param_model$pval.testBeq0[2])
    }
  }
  res = list(pvalues = pvalues)
}

#########################################################################
run_AS_surv_aalen = function(exposure,
                             survival_time,
                             censoring_status,
                             M, 
                             K,
                             covar = NULL,
                             covar_U = NULL,
                             suppl_covar = NULL,
                             each_var_pval = FALSE) {
  
  ## Check exposure is a data.frame
  check_argument_exposure(exposure) 
  
  # ## Check outcome is vector or single column data.frame
  # check_argument_outcome(outcome)
  
  ## Check survival_time
  
  ## Check censoring_status
  
  ## Check Mediator matrix is a numeric matrix
  check_argument_mediators_matrix(M)
  
  ## Check K provided and is integer
  check_K(K)
  
  if(!is.null(covar)){
    check_covar(covar)
  }
  
  if(!is.null(suppl_covar)){
    check_covar(suppl_covar)
  }
  
  # Exposure and outcome before pretreatment
  exposure_input = exposure
  survival_time_input = survival_time
  censoring_status_input = censoring_status
  
  ## Exposure data frame pretreatment
  # numeric are needed
  if (is.vector(exposure)){
    expo_var_n = 1
    expo_var_types =  typeof(exposure)
    expo_var_ids = "univariate"
    if (length(unique(exposure))<=1){
      stop("Categorial exposome must have at least two levels")
    }
    if (expo_var_types == "character"){
      #message("The input exposome is categorial")
      # model matrix transformation (-1 column to avoid colinearity)
      exposure = as.factor(exposure)
      exposure = model.matrix(~exposure)
      exposure = exposure[,-1]
      new_expo_var_type = typeof(exposure)
    } 
    else if (expo_var_types== "integer"||expo_var_types== "logical"||expo_var_types== "double"){
      #message("The input exposome is continuous or binary" )
      exposure = as.numeric(exposure)
      new_expo_var_type = typeof(exposure)
    } 
  } else if (is.factor(exposure)){
    expo_var_n = 1
    expo_var_types =  typeof(exposure)
    expo_var_ids = "univariate"
    #message("The input exposome is categorial")
    # model matrix transformation (-1 column to avoid colinearity)
    exposure = model.matrix(~exposure)
    exposure = exposure[,-1]
    new_expo_var_type = typeof(exposure)
  } else if(is.data.frame(exposure)){
    expo_var_n = dim(exposure)[2]
    expo_var_ids = colnames(exposure)
    expo_var_types = sapply(exposure, typeof)
    new_expo_var_types = list()
    exposures = c()
    for(expo_var in 1:expo_var_n) {
      if (expo_var_types[expo_var] == "character"){
        #message(paste("The input exposome no ", expo_var," is categorial"))
        # model matrix transformation
        new_exposure = as.factor(exposure[,expo_var])
        new_exposure = model.matrix(~new_exposure)
        new_exposure = new_exposure[,-1]
        new_expo_var_type =  typeof(new_exposure)
      } else if (is.factor(exposure)){ 
        #message(paste("The input exposome no ", expo_var," is categorial"))
        # model matrix transformation
        new_exposure = model.matrix(~new_exposure)
        new_exposure = new_exposure[,-1]
        new_expo_var_type =  typeof(new_exposure)
      } else if (expo_var_types[expo_var]== "integer"||expo_var_types[expo_var]== "logical"|| expo_var_types[expo_var]== "double"){
        #message(paste("The input exposome no ",expo_var, "is continuous or binary" ))
        new_exposure = exposure[,expo_var]
        new_expo_var_type = typeof(new_exposure)
      } 
      col_name = paste("Var", expo_var, sep="_")
      exposures= cbind(exposures, stats::setNames(new_exposure,col_name))
      new_expo_var_types[expo_var] = new_expo_var_type
    }
    expo_var_ids = colnames(exposure)
  } else {
    stop("Unsupported exposure variable type")
  }
  
  ## outcome pretreatment
  # exposure and outcome after pretreatment
  # if(expo_var_n == 1){
  #   if(is.data.frame(X)){
  #     exposure_output = exposures
  #   } else {
  #     exposure_output = exposure
  #   }
  # }
  # if(expo_var_n > 1){
  #   exposure_output = exposures
  # }
  # outcome_output = outcome
  
  res = list()
  ##################################
  # Run first regression : M ~ X ###
  ##################################
  # In univariate situation
  if(expo_var_n == 1){
    message("Running first regression with univariate exposure variable.")
    env_U = as.matrix(cbind(exposure, covar_U))
    if (expo_var_types == "character"|| is.factor(exposure_input)){
      mod.lfmm1 = lfmm2_med(input = M, 
                            env = env_U, 
                            K = K,
                            effect.sizes = FALSE)
      res_reg1 = lfmm2_med_test(mod.lfmm1, 
                                input = M, 
                                env = exposure,
                                covar = covar,
                                full = TRUE, #parameter to compute a single p-value for the global categorial design matrix using partial regressions
                                genomic.control = TRUE)
      if(each_var_pval == TRUE){
        #message("Generating detailed pvalues for each explanatory variable.")
        mod.lfmm1 = lfmm2_med(input = M, 
                              env = env_U, 
                              K = K,
                              effect.sizes = FALSE)
        res_reg1 = lfmm2_med_test(mod.lfmm1, 
                                  input = M, 
                                  env = exposure,
                                  full = FALSE,
                                  covar = covar,
                                  genomic.control = TRUE)
        pvals_1 = as.matrix(res_reg1$pvalues)
        names(pvals_1) = colnames(M)
      } else {
        pvals_1 = NA
      }
    } else if (is.vector(exposure) && (expo_var_types== "integer"||expo_var_types== "logical"|| expo_var_types== "double")){
      mod.lfmm1 = lfmm2_med(input = M, 
                            env = env_U, 
                            K = K,
                            effect.sizes = FALSE)
      res_reg1 = lfmm2_med_test(mod.lfmm1, 
                                input = M, 
                                env = exposure,
                                covar = covar,
                                genomic.control = TRUE)
    } else if (is.data.frame(exposure) && (expo_var_types== "integer"||expo_var_types== "logical"|| expo_var_types== "double")){
      env_U = as.matrix(cbind(exposures, covar_U))
      mod.lfmm1 = lfmm2_med(input = M, 
                            env = env_U, 
                            K = K,
                            effect.sizes = FALSE)
      res_reg1 = lfmm2_med_test(mod.lfmm1, 
                                input = M, 
                                env = exposures,
                                covar = covar,
                                genomic.control = TRUE)
    }
    pval1 = as.double(res_reg1$pvalues)
    names(pval1) = colnames(M)
    U1 = mod.lfmm1$U
    V1 = mod.lfmm1$V
    zscores1 = res_reg1$zscores
    fscores1 = res_reg1$fscores
    adj_rsquared1 = res_reg1$adj.r.squared
    gif1 = res_reg1$gif
    reg1 = list(pval1,
                U1, 
                zscores1,
                fscores1,
                adj_rsquared1,
                gif1)
    names(reg1) = c("pval","U","zscores","fscores", "adj_rsquared", "gif")
  }
  # In multivariate situation
  if(expo_var_n > 1){
    exposure = exposures
    env_U = as.matrix(cbind(exposure, covar_U))
    message("Running first regression with multivariate exposure variables.")
    # Computes a global pvalue for regression 1
    mod.lfmm1 = lfmm2_med(input = M, 
                          env = env_U, 
                          K = K,
                          effect.sizes = FALSE)
    res_reg1 = lfmm2_med_test(mod.lfmm1, 
                              input = M, 
                              env = exposure, 
                              full = TRUE, #parameter to compute a single p-value for the global multivariate model using partial regressions 
                              covar = covar,
                              genomic.control = TRUE)
    pval1 = res_reg1$pvalues
    names(pval1) = colnames(M)
    U1 = mod.lfmm1$U
    V1 = mod.lfmm1$V
    zscores1 = res_reg1$zscores
    fscores1 = res_reg1$fscores
    adj_rsquared1 = res_reg1$adj.r.squared
    gif1 = res_reg1$gif
    
    # Computed a single p-value for each explanatory variable
    if (each_var_pval == TRUE & expo_var_n == 1) {
      stop("Cannot perform detailed analysis for univariate exposome. Detailed analysis is only applicable for multivariate exposomes.")
    }
    if(each_var_pval == TRUE){
      message("Generating detailed pvalues for each explanatory variable.")
      mod.lfmm1 = lfmm2_med(input = M, 
                            env = env_U, 
                            K = K,
                            effect.sizes = FALSE)
      res_reg1 = lfmm2_med_test(mod.lfmm1, 
                                input = M, 
                                env = exposure,
                                full = FALSE,
                                covar = covar,
                                genomic.control = TRUE)
      pvals_1 = as.matrix(res_reg1$pvalues)
      names(pvals_1) = colnames(M)
    } else {
      pvals_1 = NA
    }
    reg1 = list(pval1,
                U1, 
                V1, 
                zscores1,
                fscores1,
                adj_rsquared1, 
                gif1,
                pvals_1 )
    names(reg1) = c("pval","U","V","zscores","fscores", "adj_rsquared","gif","each_var_pval")
  }
  res[[1]] = reg1  
  
  #########################################
  ### Run second regression : Y ~ X + M ###
  #########################################
  message("Running second regression.")
  if(!is.null(suppl_covar)){
    covars = cbind(covar, suppl_covar)
  } else {
    covars = covar
  }
  if(expo_var_n == 1){
    if(is.vector(exposure)||is.factor(exposure)||is.matrix(exposure)){
      res_reg2 = surv_aalen_med_test(mod.lfmm1,
                                     M = M,
                                     survival_time = survival_time,
                                     censoring_status = censoring_status,
                                     exposure = exposure,
                                     covar = covars)
    } else if(is.data.frame(exposure)){
      res_reg2 = surv_aalen_med_test(mod.lfmm1,
                                     M = M,
                                     survival_time = survival_time,
                                     censoring_status = censoring_status,
                                     exposure = exposures,
                                     covar = covars)
    }
  } else if(expo_var_n > 1){
    res_reg2 = surv_aalen_med_test(mod.lfmm1,
                                   M = M,
                                   survival_time = survival_time,
                                   censoring_status = censoring_status,
                                   exposure = exposures,
                                   covar = covars)
  }
  pval2 = as.double(res_reg2$pvalues)
  names(pval2) = colnames(M)
  #zscores2 = res_reg2$zscores
  reg2 = list(pval2 = pval2)
  # ,
  #             zscores2 = zscores2
  # )
  names(reg2) = c("pval")#, "zscores")
  res[[2]] = reg2
  
  ########################
  ### max-squared test ###
  ########################
  message("Running max-squared test.")
  # max2 test for global model
  max2_pval <- apply(cbind(pval1, pval2), 1, max)^2
  max2 = max2_pval
  res[[3]] = max2
  if (each_var_pval == TRUE){
    #message("Generating max2 pvalues for each explanatory variable.")
    max2_each_var_pval = list()
    if(expo_var_n == 1){
      for (x in 1:dim(exposure)[2]){
        max2_pval <- apply(cbind(pvals_1[x,], pval2), 1, max)^2
        names(max2_pval) = colnames(M)
        max2_each_var_pval[[colnames(exposures)[x]]] = max2_pval 
      }
      res[[4]] = max2_each_var_pval
    } else {
      for (x in 1:dim(exposures)[2]){
        max2_pval <- apply(cbind(pvals_1[x,], pval2), 1, max)^2
        names(max2_pval) = colnames(M)
        max2_each_var_pval[[colnames(exposures)[x]]] = max2_pval 
      }
      res[[4]] = max2_each_var_pval
    } 
  } else {
    #message("Not generating max2 pvalues for each explanatory variable.")
    res[[4]] = NA
  }
  input = list(
    exposure_input,
    expo_var_types,
    expo_var_ids,
    covar, 
    suppl_covar,
    survival_time_input, 
    censoring_status_input#,
    #survival_distribution
  )
  names(input) = c("exposure_input", "expo_var_types", "expo_var_ids" , "covar", "suppl_covar", "survival_time_input", "censoring_status_input")#, "survival_distribution")
  res[[5]] = input
  names(res) <- c("AS_1", "AS_2", "max2_pvalues",  "max2_each_var_pvalues", "input")
  class(res) = "hdmax2_step1"
  return(res)
}

check_argument_exposure = function(argument){
  if(is.data.frame(argument)) {
    #message("The exposure argument is a data frame")
    if (ncol(argument) == 1) {
      #message("The exposure argument is a data frame with a single column.")
    } else if (ncol(argument) > 1) {
      #message("The exposure argument is a data frame with more than one column.")
    }
  } else if (is.vector(argument)) {
    #message("The exposure argument is a vector.")
  } else if (is.factor(argument)) {
    #message("The exposure argument is a factor.")
  } else {
    stop("The exposure  is not a data frame,  nor a vector ")
  }
}

check_argument_outcome = function(argument) {
  if (is.vector(argument)) {
    #message("The outcome argument is a vector.")
  } else if (is.data.frame(argument)) {
    if (ncol(argument) == 1) {
      #message("The outcome argument is a data frame with a single column.")
    } else {
      stop("The outcome data frame must have a single column.")
    }
  } else if (is.matrix(argument)) {
    if (ncol(argument) == 1) {
      #message("The outcome matrix has a single column.")
    } else {
      stop("The outcome matrix must have a single column.")
    }
  } else {
    stop("The outcome argument is neither a vector, nor a data frame, nor a matrix with a single column.")
  }
  if (is.numeric(argument)) {
    #message("The outcome argument is numeric")
  } else if (is.integer(argument)) {
    #message("The outcome argument is integer")
  } else if (is.logical(argument)) {
    #message("The outcome argument is logical")
  } else {
    stop("The outcome argument is neither numeric, nor integer, nor logical")
  }
}

check_argument_mediators_matrix = function(argument){
  if(is.matrix(argument)){
    #message("Potential mediators matrix is actually a matrix")
  } else {
    stop("Potential mediators matrix must be a matrix")
  }
}

check_K = function(argument){
  if (!is.null(argument)) {
    #message(paste("provided K =",argument))
    if(is.integer(argument)) {
      #message("K value is integer")
    } else {
      K= as.integer(argument)
      #message("K value has been transformed as integer")
    }
  }
  else {
    stop("K is not provided")
  }
}

check_covar = function(argument){
  if (is.data.frame(argument)){
    #message("Adjutment factors is data frame")
  } else if(is.matrix(argument)){
    #message("Adjutment factors is matrix")
  } else {
    stop("Adjutment factors must be a data frame or a matrix")
  }
  for (i in 1:dim(argument)[2]){
    if(!is.numeric(argument[,1])){
      stop("adjusment factors must be numeric")
    }
  }
}

# estimate_effect_surv_aalen <- function(object , m, boots = 1000, ...) {
#
#   if (class(object)!="hdmax2_step1_surv_aalen"){
#     stop("The object is not of class hdmax2_step1. This function only compute hdmax2 objects generated by hdmax2::run_AS")
#   }
#
#
#   exposure_mat = object$input$exposure_input
#   survival_time = object$input$survival_time
#   censoring_status = object$input$censoring_status_input
#   #survival_distribution = object$input$survival_distribution
#   base = survival::Surv(survival_time, censoring_status)
#   M = m
#
#
# M = top10
# exposure = alcohol
#
#
# # ########## multi mediation using multimediate
# # lmodel.m = vector("list", length = ncol(M))
# #   for (i in 1:ncol(M)) {
# # modM = lm(M[,i] ~ exposure)
# # lmodel.m[[i]] = modM
# # }
# # lmodel.m
# #
# # mediator_vars <- paste(colnames(M), collapse = " + ")
# # formula <- as.formula(paste("Surv(survival_time, censoring_status) ~ const(exposure) +", mediator_vars))
# #
# #
# # library(timereg)
# # model.y = timereg::aalen(formula, data = as.data.frame(M))
# #
# # multi.media = multimediate::multimediate(lmodel.m,correlated=TRUE,model.y)#,treat='exposure')#,treat.value=1,control.value=0,J=1000,conf.level=0.95,data=data5)
#
# ###########################
#
# ############### univariate mediation
# for (i in 1:ncol(M)) {
#   i = 1
#   exposure = alcohol
#   survival_time = time
#   censoring_status = cens
#   data = data.frame(mediator = M[,i], exposure = exposure, time = survival_time, status = censoring_status)
#
#   mod1 = lm(mediator ~ exposure, data = data)
#   mod2 = timereg::aalen(Surv(survival_time, censoring_status) ~ timereg::const(exposure) + timereg::const(mediator), data= data)
#   lmodel.m = list(mod1)
#   model.y = mod2
#   med = multimediate::multimediate(lmodel.m, correlated=FALSE, model.y, treat = "exposure", data = data)
#
#   med = multimediate_survival(lmodel.m, correlated=FALSE, model.y, treat = "exposure", data = data)
# #
#
#   ACME[i, ] <- c(med$d0, med$d0.ci[1], med$d0.ci[2], med$d0.p)
#   ADE[i, ] <- c(med$z0, med$z0.ci[1], med$z0.ci[2], med$z0.p)
#   PM[i, ] <- c(med$n0, med$n0.ci[1], med$n0.ci[2], med$n0.p)
#   TE[i, ] <- c(med$tau.coef, med$tau.ci[1], med$tau.ci[2], med$tau.p)
# }
#
#

##' @title The function hdmax2::estimate_effect() takes as input an object
##' hdmax2_step1 and a list of potential mediators MS to be analyzed in subsequent steps.
##' @description  For each univariate exposure variable and the subset of mediators MS,
##' the hdmax2::estimate_effect() function computes several estimates
##' to evaluate the indirect effects in the path between exposure variables
##' and the survival variable. Initially, this function assesses each mediator
##' variable MS_j individually and computes causal measures of interest
##' such as (i) the Average Causal Mediated Effect (ACME, corresponding to
##' the indirect effect) and (ii) the Proportion Mediated (PM). The ACME
##' differs from the Average Direct Effect (ADE), which represents the
##' unmediated effect, and from the Total Effect (TE) which is equal to the sum
##' of direct and indirect effect. PM corresponds to the proportion of
##' the total effect that is mediated by the mediator (ratio of the indirect
##' effect to the total effect). ACME and PM are computed by the
##'  mediation::mediate() function of the package mediation, that
##' automatically detects the type of statistical model used in the mediation
##' analysis (Tingley et al. 2014). The function mediation::mediate()
##' calculates uncertainty estimates by a quasi-Bayesian Monte Carlo approach
##' described (Imai et al. 2010). In addition, it estimates the intermediary
##' effect sizes  a_j and b_j and their standard deviations. Eventually,
##' hdmax2 calculates an Overall Indirect Effect (OIE) from a single model
##' that includes all mediators MS simultaneously. The OIE corresponds
##' to the sum of the indirect effect associated with all mediators.
##' The confidence interval (CI) of the OIE is estimated by a bootstrap
##' approach. Along with the OIE, hdmax2 estimates the Overall Total
##' Effect (OTE) corresponding to the effect of exposure variables on
##' the outcome variable, as well as the Overall Direct Effect (ODE)
##' corresponding to the effect of exposure variables on the outcome variable
##' when the mediators MS are included in the model.
##'
##' @param object results from hdmax2 step 2
##' @param m a response variable matrix with n rows and p columns corresponding to mediators selected at step1.
##' Response variables must be encoded as numeric. No NAs allowed.
##' @param boots number of bootstrap
##' @param ... arguments for inner functions
##' @return hdmax2_step2_surv_param object
##'
##'  - ACME, estimation of the average causal mediation effect (the indirect effect)
##'  - ADE, estimation average direct effect
##'  - PM, estimation of the proportion mediated
##'  - TE, estimation of the total effect
##'
##' Regressions:
##'  - xm, regression X on M
##'  - my, regression M on Y
##'
##'
##' Overall effect
##'  - oie, overall indirect effect
##'  - oie_med , oie median
##'  - oie_sd , oie standard deviation
##'  - ote, overall total effect
##'  - ode, overall direct effect
##'
##' @details
##'
##' We use the mediate function of the mediation package on the set of markers having Q-value lower
##' than the FDR threshold. It estimates their indirect effects and
##' tests their significance.
##'
##' @export
##' @author Florence Pittion, Magali Richard, Olivier Francois, Basile Jumentier
##' @examples
##' # Load example dataset
##' attach(simu_data)
##' K = 5
##' # Run {hdmax2} step 1
##' hdmax2_step1 = hdmax2::run_AS(
##'   exposure = simu_data$X_binary,
##'   outcome =  simu_data$Y_continuous,
##'   M =  simu_data$M1,
##'   K = K
##' )
##' # Select mediators
##' mediators_subset = names(sort(hdmax2_step1$max2_pvalues)[1:10])
##' mediators_top10 = simu_data$M1[, mediators_subset]
##' # Run {hdmax2} step 2
##' hdmax2_step2 = hdmax2::estimate_effect(object = hdmax2_step1,
##'                                        m = mediators_top10)

estimate_effect_surv_aalen <- function(object , m, boots = 1000, ...) {
  # if (class(object)!="hdmax2_step1"){
  #   stop("The object is not of class hdmax2_step1. This function only compute hdmax2 objects generated by hdmax2::run_AS")
  # }
  library(timereg)
  exposure_mat = object$input$exposure_input
  survival_time = object$input$survival_time
  #survival_time = as.numeric(survival_time)
  censoring_status = object$input$censoring_status_input
  #censoring_status = as.numeric(censoring_status)
  #survival_distribution = object$input$survival_distribution
  base = survival::Surv(survival_time, censoring_status)
  M = m
  if (!is.matrix(M)){
    stop("m must be a matrix")
  }
  if (is.null(colnames(M))) {
    colnames(M) <- 1:ncol(M)
  }
  expo_var_ids =object$input$expo_var_ids
  ncol_var = length(expo_var_ids)
  effects = list()
  # TODO x = as.dataframe de exposure
  
  for(expo_var_id in expo_var_ids){
    if (is.vector(exposure_mat)){
      exposure = exposure_mat
    } else if (is.data.frame(exposure_mat)|| is.matrix(exposure_mat)){
      exposure = exposure_mat[, expo_var_id]
    }
    if (ncol_var == 1){
      message("Estimating indirect effect for univariate exposome.")
      if (is.null(object$input$covar)) {
        covars = data.frame(latent_factors = object$AS_1$U)
        if (is.null(object$input$covar_sup_reg2)) {
          covars_2 = data.frame(latent_factors = object$AS_1$U)
        } else  {
          covars_2 = data.frame(obs_covar_2 = object$input$covar_sup_reg2, latent_factors = object$AS_1$U)
        }
      } else  {
        covars = data.frame(obs_covar = object$input$covar, latent_factors = object$AS_1$U)
        if (is.null(object$input$covar_sup_reg2)) {
          covars_2 = data.frame(latent_factors = object$AS_1$U, obs_covar = object$input$covar)
        } else  {
          covars_2 = data.frame(obs_covar = object$input$covar, obs_covar_2 = object$input$covar_sup_reg2, latent_factors = object$AS_1$U)
        }
      }
    } else if (ncol_var > 1) {
      message("Estimating indirect effect for multivariate exposome.")
      extra_expo_vars = expo_var_ids[-which(expo_var_ids %in% expo_var_id)]
      df_extra = exposure_mat[,which(expo_var_ids %in% extra_expo_vars)]
      if(is.vector(df_extra)){
        df_extra = t(t(df_extra))
      }
      for(col in 1:dim(df_extra)[2]){
        if(typeof(df_extra[,col])=="character"){
          df_extra[,col] = as.factor(df_extra[,col])
          message(paste("Categorial column" , col , "transformed in factors in covariable data frame"))
        } else if (is.factor(df_extra[,col])) {
          message(paste("Categorial column" , col , "is factors in covariable data frame"))
        } else if (typeof(df_extra[,col])=="integer"||typeof(df_extra[,col])== "logical"||typeof(df_extra[,col])== "double"){
          df_extra[,col] = as.numeric(df_extra[,col])
          message(paste("Column" , col , "is Continuous or Binary in covariable data frame"))
        }
      }
      if (is.null(object$input$covar)) {
        covars = data.frame(latent_factors = object$AS_1$U, df_extra = df_extra)
        if (is.null(object$input$covar_sup_reg2)) {
          covars_2 = data.frame(latent_factors = object$AS_1$U, df_extra = df_extra)
        } else  {
          covars_2 = data.frame(latent_factors = object$AS_1$U, obs_covar_2 = object$input$covar_sup_reg2,  df_extra = df_extra)
        }
      } else  {
        covars = data.frame( latent_factors = object$AS_1$U, obs_covar = object$input$covar,  df_extra = df_extra)
        if (is.null(object$input$covar_sup_reg2)) {
          covars_2 = data.frame(latent_factors = object$AS_1$U, obs_covar = object$input$covar, df_extra = df_extra)
        } else  {
          covars_2 = data.frame(latent_factors = object$AS_1$U, obs_covar = object$input$covar, obs_covar_2 = object$input$covar_sup_reg2,  df_extra = df_extra)
        }
      }
    }
    
    expo_var_type =  typeof(exposure)
    if (expo_var_type == "character"||is.factor(exposure)){
      message("The input exposome is categorial")
      if(is.factor(exposure)==FALSE){
        exposure_fact = as.factor(exposure)
      } else {
        exposure_fact = exposure
      }
      exposure_dm = stats::model.matrix(~exposure_fact)
      message("categorial exposome design matrix transformation")
      exposure = exposure_dm[,-1]
      cn = colnames(exposure)

      k_effects = list()
      for (k in 1:length(cn)) {
        # from package mediation
        ACME <- matrix(ncol = 4, nrow = ncol(M))
        ADE <- matrix(ncol = 4, nrow = ncol(M))
        PM <- matrix(ncol = 4, nrow = ncol(M))
        TE <- matrix(ncol = 4, nrow = ncol(M))
        # from linear models
        xm <- matrix(ncol = 4, nrow = ncol(M))
        my <- matrix(ncol = 4, nrow = ncol(M))
        for (i in 1:ncol(M)) {#numeric
          dat.x <- data.frame(exposure_k = exposure[,k], exposure_minus_k = exposure[,-k], Mi = M[, i], covars = covars)
          dat.y <- data.frame(exposure_k = exposure[,k], exposure_minus_k = exposure[,-k], Mi = M[, i], survival_time = survival_time, censoring_status = censoring_status, covars = covars_2)
          mod1 = stats::lm(Mi ~ exposure_k + ., data = dat.x)
          message(paste0("Generate regression 1 for categorial exposure and mediator ", i))
          mod2 = timereg::aalen(Surv(survival_time, censoring_status) ~ const(exposure) + const(Mi) + const(as.vector(covars$latent_factors)), data = dat.y)
          med = multimediate::multimediate(lmodel.m = mod1, correlated=FALSE, model.y = mod2, treat = "exposure")
          
          ACME[i, ] <- c(med$d0, med$d0.ci[1], med$d0.ci[2], med$d0.p)
          ADE[i, ] <- c(med$z0, med$z0.ci[1], med$z0.ci[2], med$z0.p)
          PM[i, ] <- c(med$n0, med$n0.ci[1], med$n0.ci[2], med$n0.p)
          TE[i, ] <- c(med$tau.coef, med$tau.ci[1], med$tau.ci[2], med$tau.p)
          xm[i, ] <- summary(mod1)$coeff[2, ] # effect of X
          my[i, ] <- mod2$coefficients[3] # effect of M
          
          #### on est supposer preciser l'outcome selon doc de mediation car objet survreg ne le fait , à verifier
          med = mediation::mediate(mod1, mod2, treat = "exposure_k", mediator = "Mi")
          ACME[i, ] <- c(med$d0, med$d0.ci[1], med$d0.ci[2], med$d0.p)
          ADE[i, ] <- c(med$z0, med$z0.ci[1], med$z0.ci[2], med$z0.p)
          PM[i, ] <- c(med$n0, med$n0.ci[1], med$n0.ci[2], med$n0.p)
          TE[i, ] <- c(med$tau.coef, med$tau.ci[1], med$tau.ci[2], med$tau.p)
        }
        ACME <- as.data.frame(ACME)
        ADE <- as.data.frame(ADE)
        PM <- as.data.frame(PM)
        TE <- as.data.frame(TE)
        xm <- as.data.frame(xm)
        my <- as.data.frame(my)
        
        colnames(ACME) <- c("est", "CI_2.5", "CI_97.5", "pval")
        colnames(ADE) <- c("est", "CI_2.5", "CI_97.5", "pval")
        colnames(PM) <- c("est", "CI_2.5", "CI_97.5", "pval")
        colnames(TE) <- c("est", "CI_2.5", "CI_97.5", "pval")
        colnames(xm) <- c("Estimate", "Std.Error", "t.Value", "pValue")
        colnames(my) <- c("Estimate", "Std.Error", "z.Value", "pValue")
        
        ACME$feat <- colnames(M)
        ADE$feat <- colnames(M)
        PM$feat <- colnames(M)
        TE$feat <- colnames(M)
        xm$feat <- colnames(M)
        my$feat <- colnames(M)
        
        tmp = list(ACME=ACME, ADE=ADE, PM=PM, TE=TE, xm=xm, my=my)#, oie=oie, oie_med=oie_med, oie_sd=oie_sd, ote=ote, ode=ode)
        k_effects[[paste0("cat_", k)]] = tmp
      }
      effects[[expo_var_id]] = k_effects
    } else if (expo_var_type== "integer"||expo_var_type== "logical"||expo_var_type== "double"){
      message("The input exposome is continuous or binary" )
      # boolean transformed as numeric
      exposure = as.numeric(exposure)
      
      # To collect data from package mediation
      ACME <- matrix(ncol = 4, nrow = ncol(M))
      ADE <- matrix(ncol = 4, nrow = ncol(M))
      PM <- matrix(ncol = 4, nrow = ncol(M))
      TE <- matrix(ncol = 4, nrow = ncol(M))
      
      # To collect data from linear models
      xm <- matrix(ncol = 4, nrow = ncol(M))
      my <- matrix(ncol = 4, nrow = ncol(M))
      
      for (i in 1:ncol(M)) {#numeric
        dat.x <- data.frame(exposure = exposure, Mi = M[, i], covars = covars)
        dat.y <- data.frame(exposure = exposure, Mi = M[, i], survival_time = survival_time, censoring_status = censoring_status, latent_factor_1 = covars$latent_factors.1, latent_factor_2 = covars$latent_factors.2, latent_factor_3 = covars$latent_factors.3, latent_factor_4 = covars$latent_factors.4, latent_factor_5 = covars$latent_factors.5)#, covars = covars_2)
        
        mod1 = stats::lm(Mi ~ exposure + . , data = dat.x)
        message(paste0("Generate regression 1 for continuous or binary exposure and mediator ", i))
        lmodel.m = list(mod1)
        
        library(timereg)
        mod2 = timereg::aalen(Surv(survival_time, censoring_status) ~ const(exposure) + const(Mi) +  const(latent_factor_1) +  const(latent_factor_2) +  const(latent_factor_3) +  const(latent_factor_4) +  const(latent_factor_5) , data = dat.y)#+ ., data = dat.y)
        model.y = mod2
        
        # mod2 = timereg::aalen(as.formula(paste("Surv(survival_time, censoring_status) ~ const(exposure) + const(Mi) +",paste(colnames(covars_2), collapse = "+"))) , data = dat.y)#+ ., data = dat.y)
        # model.y = mod2
        # message(paste0("aalen done ",i))
        # med = multimediate::multimediate(lmodel.m, correlated=FALSE, model.y, treat = "exposure", treat.value = 1, control.value = 0, J = 1000,  data = dat.y)
        
        med = multimediate::multimediate(lmodel.m, correlated=FALSE, model.y, treat = "exposure", treat.value = 1, control.value = 0, J = 1000,  data = dat.y)
        
        ACME[i, ] <- c(med$d0, med$d0.ci[1], med$d0.ci[2], med$d0.p)
        ADE[i, ] <- c(med$z0, med$z0.ci[1], med$z0.ci[2], med$z0.p)
        PM[i, ] <- c(med$n0, med$n0.ci[1], med$n0.ci[2], med$n0.p)
        TE[i, ] <- c(med$tau.coef, med$tau.ci[1], med$tau.ci[2], med$tau.p)
        
        xm[i, ] <- summary(mod1)$coeff[2, ] # effect of X
        my[i, ] <- mod2$gamma[2] # effect of M
      }
      ACME <- as.data.frame(ACME)
      ADE <- as.data.frame(ADE)
      PM <- as.data.frame(PM)
      TE <- as.data.frame(TE)
      xm <- as.data.frame(xm)
      my <- as.data.frame(my)
      
      colnames(ACME) <- c("est", "CI_2.5", "CI_97.5", "pval")
      colnames(ADE) <- c("est", "CI_2.5", "CI_97.5", "pval")
      colnames(PM) <- c("est", "CI_2.5", "CI_97.5", "pval")
      colnames(TE) <- c("est", "CI_2.5", "CI_97.5", "pval")
      colnames(xm) <- c("Estimate", "Std.Error", "t.Value", "pValue")
      colnames(my) <- c("Estimate", "Std.Error", "z.Value", "pValue")
      
      ACME$feat <- colnames(M)
      ADE$feat <- colnames(M)
      PM$feat <- colnames(M)
      TE$feat <- colnames(M)
      xm$feat <- colnames(M)
      my$feat <- colnames(M)
      
      tmp = list(ACME=ACME, ADE=ADE, PM=PM, TE=TE, xm=xm, my=my)#, oie=oie, oie_med=oie_med, oie_sd=oie_sd, ote=ote, ode=ode)
      effects[[expo_var_id]] = tmp
    }
  }
  obj = list(effects = effects,
             input = object$input)
  class(obj) = "hdmax2_step2"
  return(obj)
}

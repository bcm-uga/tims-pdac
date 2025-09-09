library(MASS)
survival_data_simulation = function (n = n,
                                     p = p,
                                     #prop.causal.probes = prop.causal.probes,
                                     nb.causal.probes = nb.causal.probes,
                                     rho = rho,
                                     sd.A = sd.A,
                                     mean.A = mean.A,
                                     sd.B = sd.B,
                                     mean.B = mean.B,
                                     sigma = 0.5,## modif du 05/02
                                     #sigma = 0.1,## modif du 03/02
                                     probes_corr_value = 0.5,
                                     K=5,
                                     causal_probes_overlap = causal_probes_overlap,
                                     lambda_time=lambda_time){



  K = K
  probes_corr_value = probes_corr_value
  sigma = sigma
  #freq = NULL

  prop.causal.probes = nb.causal.probes/p
  prop.causal.ylx = causal_probes_overlap ### MAJOR MODIFICATION 05/02
  prop.causal.x = prop.causal.probes * (1 / prop.causal.ylx)
  prop.causal.y = prop.causal.x
  #prop.variance.y = 0.5
  prop.variance.x = 0.2
  sd.U = 0.1
  sd.V = 0.1
  corr_value = probes_corr_value

  idx_causal_x = idx_causal_y = idx_causal_ylx = c()

  nb_causal_bloc = nb.causal.probes/10

    # choose among 16 blocks, 10 causal probes per block (4x10, 8x10, 16x10(allblocks) => according to the number of probes) take account causal probes is intersection betxeen xm causal probes and my causal probes

  blocks=sample(1:16,nb_causal_bloc)
  A = matrix(0, p, 1)
  B = matrix(0, p, 1)

  for (i in blocks){
    nb_causal_probes_in_bloc = 10
    idx_bloc_start = (i - 1)*round(p/16) +1
    idx_bloc_end = idx_bloc_start + round(p/16)

    idx_x_start = sample(idx_bloc_start:(idx_bloc_end-nb_causal_probes_in_bloc*(2/prop.causal.ylx -1)),1)
    idx_x_end = idx_x_start + nb_causal_probes_in_bloc/prop.causal.ylx -1
    idx_y_end = idx_x_start + nb_causal_probes_in_bloc*(2/prop.causal.ylx -1)-1
    idx_y_start = idx_y_end - nb_causal_probes_in_bloc/prop.causal.ylx +1

    idx_causal_x_tmp = idx_x_start:idx_x_end # causal probes first regression
    idx_causal_y_tmp = idx_y_start:idx_y_end# causal probes second regression
    idx_causal_ylx_tmp = idx_y_start:idx_x_end# causal probes mediation (intersection of 2 above)

    # sign effects all A+/B- or all A-/B+ variable depending on the block
    signs = sample(c(-1, 1), 1)
    # if(signs==1){
    #   A[idx_causal_x_tmp, 1] = rnorm(length(idx_causal_x_tmp), mean.A, sd.A) * 1
    #   B[idx_causal_y_tmp, 1] = rnorm(length(idx_causal_y_tmp), mean.B, sd.B) * -1
    # }else{
    #   A[idx_causal_x_tmp, 1] = rnorm(length(idx_causal_x_tmp), mean.A, sd.A) * -1
    #   B[idx_causal_y_tmp, 1] = rnorm(length(idx_causal_y_tmp), mean.B, sd.B) * 1
    # }
    if(signs==1){
      A[idx_causal_x_tmp, 1] = rnorm(length(idx_causal_x_tmp), mean.A, sd.A) * 1
      B[idx_causal_y_tmp, 1] = rnorm(length(idx_causal_y_tmp), mean.B, sd.B) * 1
    }else{
      A[idx_causal_x_tmp, 1] = rnorm(length(idx_causal_x_tmp), mean.A, sd.A) * -1
      B[idx_causal_y_tmp, 1] = rnorm(length(idx_causal_y_tmp), mean.B, sd.B) * -1
    }
      idx_causal_x = c(idx_causal_x,idx_causal_x_tmp) # causal probes first regression
      idx_causal_y = c(idx_causal_y, idx_causal_y_tmp)# causal probes second regression
      idx_causal_ylx = c(idx_causal_ylx, idx_causal_ylx_tmp)

  }


  x.nb = length(idx_causal_x)
  y.nb = length(idx_causal_y)

  # if (is.null(freq))
  #   freq =  runif(n = p, min =  0.2, max =  0.8) # mean of methylation for each site

  # if (prop.variance.y + rho ^ 2 > 1)
  #   stop("prop.variance.y + rho^2 > 1")
  if (prop.variance.x + rho ^ 2 > 1)
    stop("prop.variance.x + rho^2 > 1")

  # constructing the covariance matrix
  cs.x = runif(K, min = -1, max = 1)
  #cs.x = runif(K, min = -0.5, max = 0.5)  ## modif du 03/02
  theta.x = sqrt(prop.variance.x / sum((cs.x / sd.U) ^ 2))
  Sigma = diag(x = sd.U ^ 2, nrow = K, ncol = K)
  Sigma = rbind(Sigma, matrix(cs.x * theta.x, nrow = 1))
  Sigma = cbind(Sigma, matrix(c(cs.x * theta.x, 1), ncol = 1))
  UX = MASS::mvrnorm(n, mu = rep(0, K + 1), Sigma = Sigma)
  U = UX[, 1:K, drop = FALSE]   # confounders
  X = UX[, K + 1, drop = FALSE] # exposure
  # To obtain a binary exposure, transform the continuous column X into a binary X with 0 if <0.5 and 1 if >0.5.
  X_bin = rep(NA, dim(UX)[1])
  X_bin[which(X<0)]=0
  X_bin[which(X>=0)]=1
  #latent factors loadings
  V = MASS::mvrnorm(p, mu = rep(0, K), Sigma = sd.V ^ 2 * diag(K))
  #V = rnorm(p, mu = rep(0, K), Sigma = sd.V ^ 2 * diag(K)) #rnorm simple

  # # methylation matrix generation
  Epsilon =
    apply(matrix(rep(0, p), nrow = 1), 2, function(x)
      rnorm(n, x, sigma))
  # replicate à la place de apply(suggestion Hugo)
  #hist(Epsilon)

  # Z = U %*% t(V) + X %*% t(A) + Epsilon
   Z = tcrossprod(U,V) + tcrossprod(X,A) + Epsilon # continuous exposition
  #Z = tcrossprod(U,V) + tcrossprod(X_bin,A) + Epsilon # binary exposition
  #hist(Z)

  ###############
  generate_corr_matrix = function(block_sizes, corr_value) {
    total_size = sum(block_sizes)
    corr_matrix = matrix(0, nrow = total_size, ncol = total_size)
    start_idx = 1
    for (block_size in block_sizes) {
      end_idx = start_idx + block_size - 1
      # Générer un vecteur aléatoire pour les corrélations dans le bloc
      corr_value = corr_value #runif(1, min = corr_range[1], max = corr_range[2])
      # Créer une matrice de corrélation pour le bloc
      block = matrix(corr_value, nrow = block_size, ncol = block_size)
      diag(block) = 1
      # Insérer le bloc dans la matrice principale
      corr_matrix[start_idx:end_idx, start_idx:end_idx] = block
      # Avancer à l'index suivant
      start_idx = end_idx + 1
    }
    return(corr_matrix)
  }
  ##############
  corr_value = corr_value
  nb_block = 16
  block_sizes = c(rep(p/nb_block, nb_block))
  corr_matrix = generate_corr_matrix(block_sizes, corr_value)

  # Vérifier si elle est définie positive
  is.positive.definite = function(mat) all(eigen(mat)$values > 0)
  print(is.positive.definite(corr_matrix))

  #### code d'Hugo
  data_tmp = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = corr_matrix)

  data_tmp_unif = pnorm(data_tmp) # correlated freq in [0;1]
  hist(data_tmp_unif)

  data_tmp_unif_02_08 = data_tmp_unif*0.6 + 0.2 # correlated freq in [0.2;0.8]
  # hist(data_tmp_unif_02_08)
  # # sum(data_tmp_unif_02_08 >=0.8) + sum(data_tmp_unif_02_08 <=0.2)
  # hist(data_tmp_unif[,1])
  # hist(data_tmp_unif_02_08[,1])

  M = qnorm(data_tmp_unif_02_08) + Z
  #hist(M)
  M = pnorm(M)
  #hist(M)
  dim(M)

####################
  # survival outcome generation
  u = runif(n, 0, 1)
  Time = matrix(0, n, 1)
  MX = cbind(M, X) # continuous exposure
  #MX = cbind(M, X_bin) # binary exposure
  Bx = rbind(B, matrix(rho, 1, 1))

  effects = MX %*% Bx

  ## TESTER AVEC OU SANS
  # Normalisation (optionnelle) si les valeurs sont trop extrêmes
  #effects = scale(effects)

  lambda_time = lambda_time  # Taux de base pour la survie #### A VOIR POUR FAIRE VARIER
  Time = -log(u) / (lambda_time * exp(effects))
  Time = as.numeric(Time)
  # Replacement of extreme values
  Time[Time > 150] = 150
  # 150 corresponding to 5 five year , unit is 1 month = 2.5
  #plot(density(Time))

  # Génération des temps de censure avec une distribution gamma
  shape = 2  # Ajuster pour contrôler la dispersion
  scale = quantile(Time, probs = 0.8) / shape # increase  probs to decrease censoring rate
  C = rgamma(n, shape = shape, scale = scale)

  # Calcul des temps observés et du statut de censure
  status = Time <= C
  OT = pmin(Time, C)  # Temps observés (min entre censure et survie)
  censoring_rate = 1 - mean(as.numeric(status))

  colnames(M) = paste0("M", 1:p)
  rownames(M) = paste0("S", 1:n)
  names(idx_causal_ylx) = colnames(M[,idx_causal_ylx])

  file = paste0("simu_p", p, "_n", n, "_mediators", length(idx_causal_ylx),
                "_mean.A", mean.A, "_mean.B", mean.B, "_sd.A",  sd.A, "_sd.B",
                sd.B, "_rho", rho)

  return(
    list(
      M = M,
      time  = Time,
      OT = OT,
      status = as.numeric(status),
      X = X,
      X_bin = X_bin,
      B = B,
      A = A,
      mediators = idx_causal_ylx,
      causal.x = idx_causal_x,
      causal.y = idx_causal_y,
      U = U,
      V = V,
      K = K,
      #freq = freq,
      #Sigma = Sigma,
      censoring_rate = censoring_rate,
      file = file
    )
  )

}

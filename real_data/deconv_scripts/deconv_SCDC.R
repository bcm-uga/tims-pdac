# https://github.com/meichendong/SCDC/blob/master/R/Deconvolution.R

getESET_sc <- function(ref_sc){
  exprs = ref_sc@assays$RNA$counts
  fdata = rownames(exprs)
  pdata = cbind(cell_type = ref_sc@meta.data$cell_type,
                cellname = colnames(exprs),
                sample = ref_sc@meta.data$sample,
                patient = ref_sc@meta.data$patient)
  pdata <- as.data.frame(pdata)
  fdata <- as.data.frame(fdata)
  exprs <- as.matrix(exprs)
  rownames(pdata) <- colnames(exprs)
  rownames(fdata) <- rownames(exprs)
  eset <- ExpressionSet(exprs,
                        AnnotatedDataFrame(pdata),
                        AnnotatedDataFrame(fdata))
  return(eset)
}

getESET_bulk <- function(counts){
  exprs = counts
  fdata = rownames(counts)
  pdata = cbind(cellname = colnames(counts))
  pdata <- as.data.frame(pdata)
  fdata <- as.data.frame(fdata)
  exprs <- as.matrix(exprs)
  rownames(pdata) <- colnames(exprs)
  rownames(fdata) <- rownames(exprs)
  eset <- ExpressionSet(exprs,
                        AnnotatedDataFrame(pdata),
                        AnnotatedDataFrame(fdata))
  return(eset)
}

SCDC_basis <- function(x, ct.varname, sample){
  # select only the subset of cell types of interest
  ct.sub <- unique(x@phenoData@data[,ct.varname])
  ct.sub <- ct.sub[!is.na(ct.sub)]
  x.sub <- x[,x@phenoData@data[,ct.varname] %in% ct.sub]
  # qc: remove non-zero genes
  x.sub <- x.sub[rowSums(exprs(x.sub)) > 0,]
  # calculate sample mean & sample variance matrix: genes by cell types
  countmat <- exprs(x.sub)
  ct.id <- droplevels(as.factor(x.sub@phenoData@data[,ct.varname]))
  sample.id <- as.character(x.sub@phenoData@data[,sample])
  ct_sample.id <- paste(ct.id,sample.id, sep = '%')
  mean.mat <- pbapply::pbsapply(unique(ct_sample.id), function(id){
    y = as.matrix(countmat[, ct_sample.id %in% id])
    apply(y,1,sum, na.rm = TRUE)/sum(y)
  })
  mean.id <- do.call('rbind',strsplit(unique(ct_sample.id), split = '%'))
  sigma <- pbapply::pbsapply(unique(mean.id[, 1]), function(id) {
    y = mean.mat[, mean.id[, 1] %in% id]
    if (is.null(dim(y))){
      res = rep(0, length(y))
      message("Warning: the cell type [", id,"] is only available in at most 1 subject!")
    } else {res = apply(y, 1, var, na.rm = TRUE)}
    return(res)})
  sum.mat2 <- pbapply::pbsapply(unique(sample.id), function(sid) {
    sapply(unique(ct.id), function(id) {
      y = as.matrix(countmat[, ct.id %in% id & sample.id %in% sid])
      if (ncol(y)>0){
        out = sum(y)/ncol(y)
      } else {out = 0}
      return(out)})})
  rownames(sum.mat2) <- unique(ct.id)
  colnames(sum.mat2) <- unique(sample.id)
  # library size factor calculated from the samples:
  sum.mat <- rowMeans(sum.mat2, na.rm = T)
  basis <- pbapply::pbsapply(unique(mean.id[,1]), function(id){
    z <- sum.mat[mean.id[,1]]
    mean.mat.z <- t(t(mean.mat)*z)
    y = as.matrix(mean.mat.z[,mean.id[,1] %in% id])
    apply(y,1,mean, na.rm = TRUE)})
  
  # weighted basis matrix
  my.max <- function(x, ...) {
    y <- apply(x, 1, max, na.rm = TRUE)
    if (median(y, na.rm = T) == 0){
      outx = y
    } else {outx = y/median(y, na.rm = T)}
    return(outx)
  }
  
  # MATCH DONOR, CELLTYPE, GENES!!!!!!!!!!!!!!!!
  print("Computing var.adj...")
  var.adj <- pbapply::pbsapply(unique(sample.id), function(sid) {
    my.max(sapply(unique(ct.id), function(id) {
      y = countmat[, ct.id %in% id & sample.id %in% sid,
                   drop = FALSE]
      if (ncol(y)>0){
        out = apply(y, 1, var, na.rm = T)
      } else {out = rep(0, nrow(y))}
      return(out)}), na.rm = T)})
  colnames(var.adj) <- unique(sample.id)
  q15 <- pbapply::pbapply(var.adj, 2, function(zz) {
    z1 = min(zz[zz > 0])
    z2 = quantile(zz, 0.15, na.rm = T)
    return(max(z1, z2))})
  q85 <- apply(var.adj, 2, quantile, probs = 0.85, na.rm = T)
  var.adj.q <- t(apply(var.adj, 1, function(y) {
    y[y < q15] <- q15[y < q15]
    y[y > q85] <- q85[y > q85]
    return(y)}))
  var.adj.q <- t(apply(var.adj, 1, function(y){
    y[y<q15] <- q15[y<q15]
    y[y>q85] <- q85[y>q85]
    return(y)})) #+ 1e-4
  message("Creating Basis Matrix adjusted for maximal variance weight")
  mean.mat.mvw <- pbapply::pbsapply(unique(ct_sample.id), function(id){
    sid = unlist(strsplit(id,'%'))[2]
    y = as.matrix(countmat[, ct_sample.id %in% id])
    yy = sweep(y, 1, sqrt(var.adj.q[,sid]), '/')
    apply(yy,1,sum, na.rm = TRUE)/sum(yy)})
  basis.mvw <- pbapply::pbsapply(unique(mean.id[,1]), function(id){
    z <- sum.mat[mean.id[,1]]
    mean.mat.z <- t(t(mean.mat.mvw)*z)
    y = as.matrix(mean.mat.z[,mean.id[,1] %in% id])
    apply(y,1,mean, na.rm = TRUE)})
  
  # reorder columns
  basis.mvw <- basis.mvw[,ct.sub]
  sigma <- sigma[, ct.sub]
  basis <- basis[, ct.sub]
  sum.mat <- sum.mat[ct.sub]
  return(list(basis = basis, sum.mat = sum.mat,
              sigma = sigma, basis.mvw = basis.mvw, var.adj.q = var.adj.q))
}

getCPM0 <- function(x){
  if (is.null(dim(x))){
    vec = as.matrix(x/sum(x))
    vec
  } else {
    cpm <- t(t(x)/apply(x,2,sum))
    cpm
  }
}

SCDC_qc <- function(sc.eset, ct.varname, sample, scsetname = "Single Cell",
                    ct.sub, qcthreshold = 0.7){
  iter.max = 1000
  nu = 1e-04
  epsilon = 0.001
  sc.basis = SCDC_basis(x = sc.eset, ct.varname = ct.varname, sample = sample)
  M.S <- sc.basis$sum.mat[ct.sub]
  xsc <- getCPM0(exprs(sc.eset)[rownames(sc.basis$basis.mvw),])
  N.sc <- ncol(xsc)
  m.basis <- sc.basis$basis.mvw[, ct.sub]
  sigma <- sc.basis$sigma[, ct.sub]
  valid.ct <- (colSums(is.na(sigma)) == 0) & (colSums(is.na(m.basis)) == 0) & (!is.na(M.S))
  if (sum(valid.ct) <= 1) {
    stop("Not enough valid cell type!")
  }
  message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
  m.basis <- m.basis[, valid.ct]
  M.S <- M.S[valid.ct]
  sigma <- sigma[, valid.ct]
  prop.qc <- NULL
  for (i in 1:N.sc) {
    basis.temp <- m.basis
    xsc.temp <- xsc[, i]
    sigma.temp <- sigma
    ### weighting scheme:
    lm.qc <- nnls::nnls(A=basis.temp,b=xsc.temp)
    delta <- lm.qc$residuals
    wt.gene <- 1/(nu + delta^2 + colSums((lm.qc$x)^2*t(sigma.temp)))
    x.wt <- xsc.temp*sqrt(wt.gene)
    b.wt <- sweep(basis.temp,1,sqrt(wt.gene),"*")
    lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
    prop.wt <- lm.wt$x/sum(lm.wt$x)
    delta <- lm.wt$residuals
    for (iter in 1:iter.max){
      wt.gene <- 1/(nu + delta^2 + colSums((lm.wt$x)^2*t(sigma.temp)))
      x.wt <- xsc.temp*sqrt(wt.gene)
      b.wt <- sweep(basis.temp,1,sqrt(wt.gene),"*")
      lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
      delta.new <- lm.wt$residuals
      prop.wt.new <- lm.wt$x/sum(lm.wt$x)
      if (sum(abs(prop.wt - prop.wt.new) < epsilon )){
        prop.wt <- prop.wt.new
        delta <- delta.new
        if (i%%100==0) {message(paste0("Step ",i,"/",N.sc,", converged at iteration ", iter))}
        break
      }
      prop.wt <- prop.wt.new
      delta <- delta.new
    }
    prop.qc <- rbind(prop.qc, prop.wt)
  }
  # name col and row
  colnames(prop.qc) <- colnames(m.basis)
  rownames(prop.qc) <- colnames(xsc)
  prop.qc.keep <- rowSums(prop.qc > qcthreshold) == 1 # truncated values -> F or T
  sc.eset.qc <- sc.eset[,prop.qc.keep]
  return(list(prop.qc = prop.qc, sc.eset.qc = sc.eset.qc))
}

SCDC_prop <- function (bulk.eset, sc.eset, ct.varname, sample, ct.sub, iter.max = 1000){
  nu = 1e-04
  epsilon = 0.001
  bulk.eset <- bulk.eset[rowSums(exprs(bulk.eset)) > 0, , drop = FALSE]
  ct.sub <- intersect(ct.sub, unique(sc.eset@phenoData@data[, ct.varname]))
  sc.basis <- SCDC_basis(x = sc.eset, ct.varname = ct.varname, sample = sample)
  commongenes <- intersect(rownames(sc.basis$basis.mvw), rownames(bulk.eset))
  if (length(commongenes) < 0.2 * min(dim(sc.eset)[1], dim(bulk.eset)[1])) {
    stop("Too few common genes!")
  }
  message(paste("Used", length(commongenes), "common genes..."))
  basis.mvw <- sc.basis$basis.mvw[commongenes, ct.sub]
  xbulk <- getCPM0(exprs(bulk.eset)[commongenes, , drop = F])
  sigma <- sc.basis$sigma[commongenes, ct.sub]
  ALS.S <- sc.basis$sum.mat[ct.sub]
  N.bulk <- ncol(bulk.eset)
  valid.ct <- (colSums(is.na(sigma)) == 0) & (colSums(is.na(basis.mvw)) == 0) & (!is.na(ALS.S))
  if (sum(valid.ct) <= 1) {
    stop("Not enough valid cell type!")
  }
  message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
  basis.mvw <- basis.mvw[, valid.ct]
  ALS.S <- ALS.S[valid.ct]
  sigma <- sigma[, valid.ct]
  prop.est.mvw <- NULL
  yhat <- NULL
  yhatgene.temp <- rownames(basis.mvw)
  for (i in 1:N.bulk) {
    basis.mvw.temp <- basis.mvw
    xbulk.temp <- xbulk[, i]*100
    sigma.temp <- sigma
    #message(paste(colnames(xbulk)[i], "has common genes",
    #              sum(xbulk[, i] != 0), "..."))
    lm <- nnls::nnls(A = basis.mvw.temp, b = xbulk.temp)
    delta <- lm$residuals
    wt.gene <- 1/(nu + delta^2 + colSums((lm$x * ALS.S)^2 *
                                           t(sigma.temp)))
    x.wt <- xbulk.temp * sqrt(wt.gene)
    b.wt <- sweep(basis.mvw.temp, 1, sqrt(wt.gene), "*")
    lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
    prop.wt <- lm.wt$x/sum(lm.wt$x)
    delta <- lm.wt$residuals
    for (iter in 1:iter.max) {
      wt.gene <- 1/(nu + delta^2 + colSums((lm.wt$x * ALS.S)^2 *
                                             t(sigma.temp)))
      x.wt <- xbulk.temp * sqrt(wt.gene)
      b.wt <- sweep(basis.mvw.temp, 1, sqrt(wt.gene), "*")
      lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
      delta.new <- lm.wt$residuals
      prop.wt.new <- lm.wt$x/sum(lm.wt$x)
      if (sum(abs(prop.wt.new - prop.wt)) < epsilon) {
        prop.wt <- prop.wt.new
        delta <- delta.new
        R2 <- 1 - var(xbulk.temp - basis.mvw.temp %*%
                        as.matrix(lm.wt$x))/var(xbulk.temp)
        if (i%%10==0) {message(paste0("Step ",i,"/",N.bulk,", WNNLS Converged at iteration ",iter))}
        break
      }
      prop.wt <- prop.wt.new
      delta <- delta.new
    }
    R2 <- 1 - var(xbulk.temp - basis.mvw.temp %*% as.matrix(lm.wt$x))/var(xbulk.temp)
    prop.est.mvw <- rbind(prop.est.mvw, prop.wt)
    yhat.temp <- basis.mvw.temp %*% as.matrix(lm.wt$x)
    yhatgene.temp <- intersect(rownames(yhat.temp), yhatgene.temp)
    yhat <- cbind(yhat[yhatgene.temp, ], yhat.temp[yhatgene.temp,
    ])
  }
  colnames(prop.est.mvw) <- colnames(basis.mvw)
  rownames(prop.est.mvw) <- colnames(xbulk)
  colnames(yhat) <- colnames(xbulk)
  yobs <- exprs(bulk.eset)
  return(list(prop.est.mvw = prop.est.mvw, basis.mvw = basis.mvw,
              yhat = yhat))
}

SCDC_ENSEMBLE <- function(bulk.eset, sc.eset.list = NULL, ct.varname, sample, ct.sub){
  iter.max = 2000
  # STEP 1: CALCULATE PROPORTION USING REF SEPARATELY.
  prop.list <- lapply(sc.eset.list, function(zz){
    SCDC_prop(bulk.eset = bulk.eset, sc.eset = zz, ct.varname = ct.varname, sample = sample,
                ct.sub = ct.sub, iter.max = iter.max)})
  row.list <- sapply(1:length(prop.list), function(x){
    rownames(prop.list[[x]]$yhat)})
  gene.prop <- Reduce("intersect", row.list)
  gene.prop2 <- intersect(gene.prop, rownames(bulk.eset))
  subj.order <- colnames(bulk.eset)
  ycpm <- getCPM0(exprs(bulk.eset)[gene.prop2,subj.order])
  g.filter <- rowSums(ycpm) < quantile(rowSums(ycpm), 0.95) & rowSums(ycpm) > quantile(rowSums(ycpm), 0.15)
  gene.use <- gene.prop2[g.filter]
  length(gene.use)
  yv <- c(getCPM0(exprs(bulk.eset)[gene.use,subj.order]))*1e5 #vectorize y to make computing faster. scale to 100,000
  y.list <- do.call(cbind, lapply(prop.list, function(x){
    c(getCPM0(x$yhat[gene.use,subj.order]))*1e5}))
  sse <- function(x,y) {
    sum((x-y)^2, na.rm = T)
  }
  sae <- function(x,y) {
    sum(abs(x-y), na.rm = T)
  }
  rmsd <- function(x,y) {
    sqrt(mean((x-y)^2, na.rm = T))
  }
  
  # STEP2: WEIGHT CALCULATION
  # -------------------------------
  # sum of squared errors
  message("Searching ENSEMBLE weight by Sum of Squared Errors or Sum of Abs Errors ......")
  sses <- apply(y.list, 2, function(x){
    sse(yv, x)
  })
  sse.wt <- 1/sses / sum(1/sses)
  # sum of absolute errors
  saes <- apply(y.list, 2, function(x){
    sae(yv, x)
  })
  sae.wt <- 1/saes / sum(1/saes)
  # RMSD
  rmsds <- apply(y.list, 2, function(x){
    rmsd(yv, x)
  })
  rmsd.wt <- 1/rmsds / sum(1/rmsds)
  #-------------------------
  # combo of proportions according to selected weights
  combo <- lapply(prop.list, function(x){
    x$prop.est.mvw
  })
  # ----------------------------------------------
  # summarize all weights and performances
  weight.mat <- rbind(
    sse.wt,
    sae.wt,
    rmsd.wt)
  rownames(weight.mat) <- c("inverse SSE", "inverse SAE","inverse RMSD")
  wt_y <- function(wt, y.list = y.list){
    wt <- as.numeric(wt)
    combo.list <- list()
    for (i in 1:ncol(y.list)){
      combo.list[[i]] <- y.list[,i]*wt[i]
    }
    combo.y <- Reduce("+", combo.list)
    return(combo.y)
  }
  out <- round(weight.mat,2)
  return(list(w_table = out, prop.only = combo))
}

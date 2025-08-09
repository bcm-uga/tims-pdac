##' @title Identification of differentially methylated regions (AMR)
##' @description  Function adapt from the combp function() of the ENmix package
##' @param data A data frame from bed format file with colname name
##' "V1","V2", "V3","V4","V5",V1 indicate chromosome (1,2,3,...,X,Y),
##' V2 is chromosome position, V4 is for P value and V5 for name of CpGs.
##' @param dist.cutoff Maximum distance in base pair to combine adjacent AMRs.
##' @param bin.size bin size for autocorrelation calculation.
##' @param seed FDR significance threshold for initial selection of AMR region.
##' @param nCores Number of computer cores used in calculation
##' @return
##' Results of the AMRs analysis.
##' - result.fdr, table of selected AMRs. For each AMR include chromosomic position, P-value, and FDR
##' @export
##' @details
##' The input should be a data frame with column name V1-V5, indicating chromosome, start position,end position,
##' pValues and probe names. The function will use a modified comb-p method to identify
##' differentially methylated regions.
##' @author Basile Jumentier
##' @example 
##' data = hdmax2::helper_ex
##' chr = data$annotation$chr
##' start = data$annotation$start
##' end = data$annotation$end
##' pval = hdmax2_step1$max2_pvalues
##' cpg = data$annotation$cpg
##' data_combp <- data.frame(chr, start, end, pval, cpg)
##' colnames(data_combp) <- paste0("V", 1:5)
##' res <- combp2(data_combp, seed = 0.6,  nCores = 2, ...)
##'
combp2 <- function (data, dist.cutoff = 1000, bin.size = 310, seed = 0.01, nCores = 10) {
  
  ##### a function to get a table of p-values for estimating acf
  #####loc should be increasing;
  acf.table<-function(x,loc,dist.cutoff){
    flag=TRUE; lag=1; result=NULL
    while(flag){
      x1=utils::head(x,-lag); x2=utils::tail(x,-lag); dist=diff(loc,lag=lag)
      index=(dist<dist.cutoff)
      if(all(!index)){flag=FALSE}else{
        result=rbind(result,data.frame(x1=x1[index],x2=x2[index],dist=dist[index]))
        lag=lag+1
      }
    }
    return(result)
  }
  
  ##### a function to estimate acf
  get.acf <- function(data, dist.cutoff, bin.size) {
    temp <- NULL
    for (chr in unique(data$V1)) {
      y <- data[data$V1 == chr, ]
      y <- y[order(y$V3), ]
      temp <- rbind(temp, acf.table(y$V4, y$V3, dist.cutoff))
    }
    
    bin.label <- findInterval(temp$dist, seq(bin.size, dist.cutoff, bin.size))
    # new line add try catch
    temp.stouffer <- by(temp, bin.label, FUN = function(x) {
      tryCatch(
        stats::cor.test(stats::qnorm(x$x1), stats::qnorm(x$x2), alternative = "greater"),
        error = function(e) return(NA)
      )
    }, simplify = FALSE)
    # add line when cor test failed cor = 0 and pvalue = 1
    cor.stouffer <- sapply(temp.stouffer, function(x) {
      if (inherits(x, "htest") && !is.null(x$estimate)) {
        unname(x$estimate)
      } else {
        0
      }
    })
    
    p.stouffer <- sapply(temp.stouffer, function(x) {
      if (inherits(x, "htest") && !is.null(x$p.value)) {
        x$p.value
      } else {
        1
      }
    })
    
    if (any(p.stouffer > 0.05, na.rm = TRUE)) {
      index = min(which(p.stouffer > 0.05))
      cor.stouffer[index:length(cor.stouffer)] = 0
    }
    
    return(cor.stouffer)
  }
  
  if (nCores > parallel::detectCores()) {
    nCores = parallel::detectCores()
  }
  data = as.data.frame(data)
  acf <- get.acf(data, dist.cutoff, bin.size)
  result <- parallel::mclapply(unique(data$V1), function(chr) {
    y = data[data$V1 == chr, ]
    y = y[order(y$V3), ]
    pos = y$V3
    p = stats::qnorm(y$V4)
    temp = sapply(pos, function(i) {
      index.i = (abs(pos - i) < bin.size)
      if (sum(index.i) > 1) {
        int <- findInterval(c(stats::dist(pos[index.i])), c(bin.size,
                                                            2 * bin.size))
        sd <- sqrt(sum(acf[int + 1]) * 2 + sum(index.i))
        return(stats::pnorm(sum(p[index.i]), mean = 0, sd = sd))
      }
      else {
        return(y$V4[index.i])
      }
    })
    return(data.frame(chr, start = pos, end = pos, s.p = temp))
  }, mc.cores = nCores)
  result <- do.call("rbind", result)
  names(result) = c("chr", "start", "end", "s.p")
  result = result[stats::p.adjust(result$s.p, method = "fdr") < seed,]
  result = result[!is.na(result$start) & !is.na(result$end), ]  # ligne ajoutée
  
  result.fdr = NULL
  if (nrow(result) > 0) {
    for (chr in unique(result$chr)) {
      y = data[data$V1 == chr, ]
      y = y[order(y$V3), ]
      pos = y$V3
      p = stats::qnorm(y$V4)
      result.chr = result[result$chr == chr, ]
      result.chr = result.chr[!is.na(result.chr$start) & !is.na(result.chr$end), ]  # ligne ajoutée
      if (nrow(result.chr) == 0) next  # passe si plus rien  # ligne ajoutée
      a = IRanges::IRanges(start = result.chr$start, end = result.chr$end)
      b = IRanges::reduce(a, min.gapwidth = dist.cutoff)
      start = IRanges::start(b)
      end = IRanges::end(b)
      region.max <- max(Biostrings::width(b))
      temp = sapply(1:length(b), function(i) {
        index.i = (pos >= start[i] & pos <= end[i])
        
        # print(sum(index.i))
        
        if (sum(index.i) > 1) {
          int <- findInterval(c(stats::dist(pos[index.i])),
                              seq(bin.size, region.max + bin.size, bin.size))
          sd <- sqrt(sum(ifelse(int < length(acf), acf[int +
                                                         1], 0)) * 2 + sum(index.i))
          return(stats::pnorm(sum(p[index.i]), mean = 0, sd = sd))
        }
        else {
          return(y$V4[index.i])
        }
      })
      result.fdr = rbind(result.fdr, data.frame(chr, start,
                                                end, p = temp))
    }
    result.fdr$fdr = stats::p.adjust(result.fdr$p, method = "fdr")
    result.fdr <- result.fdr[order(result.fdr$p), ]
    result.fdr$start = (result.fdr$start - 1)
  }
  
  return(result.fdr)
}

#------------------------------------------------------------------

##' @title Identifying aggregated mediator regions (AMR)
##'
##' @description Identify aggregated methylated regions (AMR) from the P-values from function max2 using a modified comb-p method. 
##' Compute the P-value and the FDR for each AMR detected.
##' @param chr chromosomes
##' @param start chromosomal position of markers (start)
##' @param end chromosomal position of markers (end)
##' @param pval P-values for each markers, from the max2 function
##' @param cpg name of each markers
##' @param ... see help of combp of ENmix package
##' @return
##' - res, table of selected AMRs. For each AMR include chromosomic position, P-value, and FDR
##' - data, matrix of all cpg, with annotation and provided P-values
##' @details
##' The function uses a modified comb-p method to identify
##' aggregated methylated regions (AMRs).
##' @export
##' @author Basile Jumentier
##' @examples
##' data = hdmax2::helper_ex
##' K=5
##' ## run hdmax2 step1
##' hdmax2_step1 = hdmax2::run_AS(exposure = data$exposure,
##'                              outcome = data$phenotype,
##'                              M = data$methylation,
##'                              K = K)
##'
##' ##Detecting AMR
##' chr = data$annotation$chr
##' start = data$annotation$start
##' end = data$annotation$end
##' pval = hdmax2_step1$max2_pvalues
##' cpg = data$annotation$cpg
##'
##' res.amr_search = hdmax2::AMR_search(
##' chr = data$annotation$chr,
##' start = data$annotation$start,
##' end = data$annotation$end,
##' pval = hdmax2_step1$max2_pvalues,
##' cpg = data$annotation$cpg,
##' seed = 0.7, #Careful to change this parameter when working with real data
##' nCores = 2)
##' res.amr_search$res
##'

AMR_search <- function(chr, start, end, pval, cpg, ...) {
  
  tmp <- data.frame(chr, start, end, pval, cpg)
  colnames(tmp) <- paste0("V", 1:5)
  
  tmp <- combp2(tmp, ...)
  
  return(list(res = tmp,
              data = data.frame(chr, start, end, pval, cpg)))
}

#------------------------------------------------------------------

##' @title Build AMR vector
##'
##' @description Build AMR from the result of function AMR_search
##' @param res result object of function AMR_search
##' @param methylation a matrix of methylation profile.
##' @param nb_cpg threshold of minimal number of CpG in the AMR
##' @return
##' A set of build AMRs.
##'  - res, selected AMR
##'  - CpG_for_each_AMR, list of markers present on each AMR.
##'  - AMR_acp, first components of PCA for each AMR
##'  - AMR_mean, mean value of CpG on the AMR 
##' @details
##' We use the series of pValues (one pValue per CpGs) obtained with the mEWAS
##' regression method and the combination of pValue max2.
##' To determine the potential AMRs used the combp method present in the ENmix package (Xu et al. 2016).
##' This method uses the Fisher method to combine the pValues and also the base pair distance (bP)
##' between CpGs (1000 bP maximum between nb_cpg CpGs on the same AMR).
##' The information for each AMR is summarized by doing the mean (by row) of each CpG.
##' @export
##' @author Basile Jumentier
##' @examples
##' data = hdmax2::helper_ex
##' K=5
##' ## run hdmax2 step1
##' hdmax2_step1 = hdmax2::run_AS(exposure = data$exposure,
##'                              outcome = data$phenotype,
##'                              M = data$methylation,
##'                              K = K)
##'
##'##Detecting AMR
##' chr = data$annotation$chr
##' start = data$annotation$start
##' end = data$annotation$end
##' pval = hdmax2_step1$max2_pvalues
##' cpg = data$annotation$cpg
##'
##' res.amr_search = hdmax2::AMR_search(
##' chr = data$annotation$chr,
##' start = data$annotation$start,
##' end = data$annotation$end,
##' pval = hdmax2_step1$max2_pvalues,
##' cpg = data$annotation$cpg,
##' seed = 0.7, #Careful to change this parameter when working with real data
##' nCores = 2)
##'
##' res.amr_search$res
##'
##' res.arm_build = hdmax2::AMR_build(res.amr_search, 
##' methylation = data$methylation, nb_cpg = 2)
##' #List of AMR selected
##' head(res.arm_build$res)
##' ## CpG in the AMR
##' res.arm_build$CpG_for_each_AMR
##'
AMR_build <- function(res, methylation, nb_cpg = 2) {
  
  data <- res$data
  res <- res$res
  
  # Number of CpG per AMR
  
  nb <- NULL
  
  for (i in 1:nrow(res)) {
    chri <- as.character(res$chr[i])
    tmp <- dplyr::filter(data, chr == chri)
    nb <- c(nb, sum((res$start[i]:res$end[i]) %in% tmp$start))
  }
  
  # Select AMRs with nb_cpg CpGs at minimum
  
  res <- cbind(res, nb)
  res <- dplyr::filter(res, nb >= nb_cpg)
  AMR.select <- list()
  
  for (i in 1:nrow(res)) {
    chri <- as.character(res$chr[i])
    tmp <- dplyr::filter(data, chr == chri)
    AMR.select[[i]] <- as.character(tmp$cpg[(tmp$start %in% (res$start[i]:res$end[i]))])
  }
  
  # Select CpGs values in the methylation matrix
  
  AMR.meth <- list()
  
  for (i in 1:length(AMR.select)) {
    AMR.meth[[i]] <- methylation[, AMR.select[[i]]]
  }
  
  # Built a vector for each AMR with the first component of PCA or with the rowmeans
  
  AMR.acp <- as.data.frame(matrix(ncol = length(AMR.meth), nrow = nrow(methylation)))
  colnames(AMR.acp) <- paste0("AMR", 1:length(AMR.meth))
  
  AMR.mean <- as.data.frame(matrix(ncol = length(AMR.meth), nrow = nrow(methylation)))
  colnames(AMR.mean) <- paste0("AMR", 1:length(AMR.meth))
  
  for (i in 1:length(AMR.meth)) {
    AMR.acp[, i] <- prcomp(AMR.meth[[i]])$x[, 1]
    AMR.mean[, i] <- rowMeans(AMR.meth[[i]])
  }
  
  # data
  
  res <- cbind(AMR = colnames(AMR.acp), res)
  names(AMR.select) <- colnames(AMR.acp)
  
  return(list(AMR_acp = AMR.acp,
              AMR_mean = AMR.mean,
              res = res,
              CpG_for_each_AMR = AMR.select))
}



##' @title Compute q-values from p-values
##' @description Compute q-values from provided p-values Estimate the proportion of H0
##' relies on the defined interquartile range
##' @param pvalues result of max squared test
##' @param theta constant defines interquartile range
##' @return qvalues, pi0 =  H0_proportion
##' @export
##' @author Olivier Francois
##' @examples
##' data = hdmax2::simu_data
##' K = 5
##' hdmax2_step1 = hdmax2::run_AS(exposure = simu_data$X_binary,
##'                               outcome = simu_data$Y_continuous,
##'                               M = simu_data$M1, 
##'                               K = K)
##' 
##'  #Select candidate mediator  
##' qv = hdmax2::hdmax2_qvalue(hdmax2_step1$max2_pvalues)
##' 
##' 
##' 
hdmax2_qvalue <- function(pvalues, theta = 0.25){
  
  ##Estimate the proportion of H0
  ##relies on the interquartile range (25%-75%) of the default histogram
  H0_proportion = 2*mean( pvalues > theta & pvalues < (theta + 0.5) )
  
  ## compute tests number
  test_number = length(pvalues)
  
  ## qvalues definition
  qvalues =  H0_proportion * test_number * pvalues / rank(pvalues)
  
  return(list(qvalues =  qvalues , pi0 =  H0_proportion))
}
#' Estimate the confidence interval of a heritability estimate.
#'
#' @param data matrix with the phenotype in the first column and the covariate model matrix
#'       in the remaining columns. Sample IDs should be in rownames(data).
#' @param K list of per chromosome kinship matrices. 
#' @param nsim integer indicating the number of simulations to perform. Default = 1000.
#'
#' @return numeric matrix containing heritabilities estimates for the observed data (in row 1) and the simulated data.
#' 
#' @importFrom regress regress
#' @importFrom mvtnorm rmvnorm
#' @export
herit_per_chr_ci = function(data, K, nsim = 1000) {

  K = qtl2::double_kinship(K)

  h2_chr = h2_per_chr(data, K)
  
  n = nrow(data)

  # Multiply each kinship matrix by the heritability and sum it.
  tmp = mapply(function(w, k) { w * k }, h2_chr, K, SIMPLIFY = FALSE)
  tmp = apply(vapply(tmp, "+", matrix(0, n, n)), 1:2, sum)
  # Add the error variance component as well.
  err_cov = tmp + (1 - sum(h2_chr)) * diag(n)

  # Simulate nsim phenotypes and estimate their heritability.
  pheno_sim = t(rmvnorm(n = nsim, mean = rep(0, n), sigma = err_cov))
  rownames(pheno_sim) = rownames(data)
  sim_data = vector("list", ncol(pheno_sim))

  for(i in 1:length(sim_data)) {

    sim_data[[i]] = cbind(pheno_sim[,i], data[,-1])

  } # for(i)

  # This returns a list of vectors. Each vector contains per-chromosome
  # heritabilities.
  sim_h2 = lapply(sim_data, h2_per_chr)
  sim_h2 = matrix(unlist(sim_h2), nrow = length(sim_h2), ncol = length(sim_h2[[1]]),
            dimnames = list(paste0("sim", 1:length(sim_h2)), names(sim_h2[[1]])), 
            byrow = TRUE)
  sim_h2 = rbind(h2_chr, sim_h2)
  rownames(sim_h2)[1] = "obs"
  colnames(sim_h2) = c(1:19, "X")

  return(sim_h2)

} # herit_per_chr_ci()
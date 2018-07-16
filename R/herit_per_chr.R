#' Estimate the heritability of a phenotype on each chromosome.
#' 
#' @param data numeric matrix containing samples in rows and phenotypes and covariates in columns. The phenotype should be in column 1 and the covariates should be generated from model.matrix() in the remaining columns.
#' @param K list containing kinship matrices for each chromosome. NOTE: these should be made using qtl2::calc_kinship(genoprobs, type = "chr")
#' 
#' @return numeric matrix containing phenotypes in rows and kinship estimates per chromosome in columns.
#' 
#' @importFrom regress regress
#' @export
herit_per_chr = function(data, K) {

  K = qtl2::double_kinship(K)

  mod = regress(data[,1] ~ data[,-1], 
        ~K[[1]]  + K[[2]]  + K[[3]]  + K[[4]]  + K[[5]] + 
         K[[6]]  + K[[7]]  + K[[8]]  + K[[9]]  + K[[10]] + 
         K[[11]] + K[[12]] + K[[13]] + K[[14]] + K[[15]] +
         K[[16]] + K[[17]] + K[[18]] + K[[19]] + K[[20]],
         pos = rep(TRUE, length(K) + 1))

  return(mod$sigma[1:length(K)] / sum(mod$sigma))

} # h2_per_chr()



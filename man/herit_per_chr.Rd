% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/herit_per_chr.R
\name{herit_per_chr}
\alias{herit_per_chr}
\title{Estimate the heritability of a phenotype on each chromosome.}
\usage{
herit_per_chr(data, K)
}
\arguments{
\item{data}{numeric matrix containing samples in rows and phenotypes and covariates in columns. The phenotype should be in column 1 and the covariates should be generated from model.matrix() in the remaining columns.}

\item{K}{list containing kinship matrices for each chromosome. NOTE: these should be made using qtl2::calc_kinship(genoprobs, type = "chr")}
}
\value{
numeric matrix containing phenotypes in rows and kinship estimates per chromosome in columns.
}
\description{
Estimate the heritability of a phenotype on each chromosome.
}

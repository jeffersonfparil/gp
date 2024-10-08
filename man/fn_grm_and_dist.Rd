% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distances.R
\name{fn_grm_and_dist}
\alias{fn_grm_and_dist}
\title{Genetic relationship matrices and their inverses}
\usage{
fn_grm_and_dist(
  G,
  maf = 0.01,
  ploidy = 2,
  diagonal_load = 0.001,
  verbose = FALSE
)
}
\arguments{
\item{G}{n samples x p loci matrix of allele frequencies (numeric ranging from 0 to 1) with non-null row and column names.
Row names can be any string of characters which identify the sample or entry or pool names.
Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name,
the second should be numeric which refers to the position in the chromosome/scaffold, and
subsequent elements are optional which may refer to the allele identifier and other identifiers.}

\item{maf}{minimum allele frequency (Default=0.01)}

\item{ploidy}{(for the GRM of Bell et al, 2017 & VanRaden et al, 2008) ploidy level which can refer to the number of haploid genomes to simulate pools (Default=2)}

\item{diagonal_load}{a small positive value to add to the diagonals of the genetic relationship matrix to ensure invertibility (Default=0.001)}

\item{verbose}{show messages? (Default=FALSE)}
}
\value{
\itemize{
\item Ok (returns several genetic relationship and distances matrices):
\itemize{
\item $grm: simple X\%*\%t(X)/ncol(X)
\item $grm_VanRaden: ploidy-aware
\item $dist_binary: Jaccard's distance, where non-zero values are converted to ones.
\item $dist_euclidean: Euclidean's distance (L2 norm)
\item $dist_taxicab: Taxicab or Manhattan distance (L1 norm)
\item $inverse_grm: inverse of $grm
\item $inverse_grm_VanRaden: inverse of $grm_VanRaden
\item $inverse_one_minus_dist_binary: inverse of 1.00 - $dist_binary
\item $inverse_one_minus_dist_euclidean: inverse of 1.00 - $dist_euclidean
\item $inverse_one_minus_dist_taxicab: inverse of 1.00 - $dist_taxicab
}
\item Err: gpError
}
}
\description{
Genetic relationship matrices and their inverses
}
\examples{
list_sim = fn_simulate_data(verbose=TRUE)
G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
list_grm_and_dist = fn_grm_and_dist(G=G, verbose=TRUE)
}

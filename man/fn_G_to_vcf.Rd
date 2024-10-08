% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{fn_G_to_vcf}
\alias{fn_G_to_vcf}
\title{Convert numeric allele frequency matrix into a vcfR::vcf object with randomly sampled depths (for fn_simulate_data(...) below)}
\usage{
fn_G_to_vcf(G, min_depth = 100, max_depth = 1000, verbose = FALSE)
}
\arguments{
\item{G}{numeric n samples x p biallelic loci matrix of allele frequencies with non-null row and column names.
Row names can be any string of characters which identify the sample or entry or pool names.
Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name,
the second should be numeric which refers to the position in the chromosome/scaffold, and
subsequent elements are optional which may refer to the allele identifier and other identifiers.
Note that this function only accepts biallelic loci.}

\item{min_depth}{minimum depth per locus for sampling simulated depths (Default=100)}

\item{max_depth}{maximum depth per locus for sampling simulated depths (Default=1000)}

\item{verbose}{show allele frequency genotype matrix to vcf conversion messages? (Default=FALSE)}
}
\value{
\itemize{
\item Ok: simulated genotype data as a vcfR object with GT, AD and DP fields
\item Err: grError
}
}
\description{
Convert numeric allele frequency matrix into a vcfR::vcf object with randomly sampled depths (for fn_simulate_data(...) below)
}
\examples{
G = simquantgen::fn_simulate_genotypes(verbose=TRUE)
vcf = fn_G_to_vcf(G=G, verbose=TRUE)
}

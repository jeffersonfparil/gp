% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load.R
\name{fn_vcf_to_G}
\alias{fn_vcf_to_G}
\title{Convert biallelic vcf data into allele frequencies where loci beyond minimum and maximum depths are set to missing}
\usage{
fn_vcf_to_G(
  vcf,
  min_depth = 0,
  max_depth = Inf,
  retain_minus_one_alleles_per_locus = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{vcf}{biallelic genotype data as a vcfR object with GT and/or AD and DP fields (where AD takes precedence over GT)}

\item{min_depth}{minimum depth per locus beyond which will be set to missing data \link{Default=0}}

\item{max_depth}{maximum depth per locus beyond which will be set to missing data \link{Default=Inf}}

\item{retain_minus_one_alleles_per_locus}{omit the alternative or trailing allele per locus? \link{Default=TRUE}}

\item{verbose}{show vcf to allele frequency genotype matrix conversion messages? \link{Default=FALSE}}
}
\value{
Ok: named n samples x p biallelic loci matrix.
Row names can be any string of characters.
Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name,
the second should be numeric which refers to the position in the chromosome/scaffold, and
subsequent elements are optional which may refer to the allele identifier and other identifiers.
Err: gpError
}
\description{
Convert biallelic vcf data into allele frequencies where loci beyond minimum and maximum depths are set to missing
}
\examples{
G = simquantgen::fn_simulate_genotypes(verbose=TRUE)
vcf = fn_G_to_vcf(G=G, verbose=TRUE)
G_back = fn_vcf_to_G(vcf=vcf, verbose=TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{fn_filter_phenotype}
\alias{fn_filter_phenotype}
\title{Filter phenotype data by converting outliers into missing data and optionally excluding missing data}
\usage{
fn_filter_phenotype(
  list_pheno,
  remove_outliers = TRUE,
  remove_NA = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{list_pheno}{list with 3 elements:
(1) $y: a named vector of numeric phenotype data,
(2) $pop: population or groupings corresponding to each element of y, and
(3) $trait_name: name of the trait.}

\item{remove_outliers}{remove missing data in the y vector?
If true, this removes the corresponding element/s in the pop vector. (Default=TRUE)}

\item{remove_NA}{remove missing data in the y vector?
If true, this removes the corresponding element/s in the pop vector. (Default=FALSE)}

\item{verbose}{show phenotype filtering messages? (Default=FALSE)}
}
\value{
\itemize{
\item Ok:
\itemize{
\item $y: named vector of numeric phenotype data
\item $pop: population or groupings corresponding to each element of y
\item $trait_name: name of the trait
}
\item Err: gpError
}
}
\description{
Filter phenotype data by converting outliers into missing data and optionally excluding missing data
}
\examples{
list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
list_pheno$y[1] = Inf
list_pheno$y[2] = NA
list_pheno_filtered = fn_filter_phenotype(list_pheno, remove_NA=TRUE, verbose=TRUE)
}

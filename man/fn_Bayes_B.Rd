% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{fn_Bayes_B}
\alias{fn_Bayes_B}
\title{Bayes B regression model
(scaled t-distributed effects with probability \eqn{\pi}; and zero effects with probability \eqn{1-\pi},
where \eqn{\pi \sim \beta(\theta_1, \theta_2)})}
\usage{
fn_Bayes_B(
  list_merged,
  vec_idx_training,
  vec_idx_validation,
  other_params = list(nIter = 12000, burnIn = 2000, out_prefix = "bglr_bayesB-"),
  verbose = FALSE
)
}
\arguments{
\item{list_merged}{list of merged genotype matrix, and phenotype vector, as well as an optional covariate matrix
\itemize{
\item $G: numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
Row names can be any string of characters which identify the sample or entry or pool names.
Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name,
the second should be numeric which refers to the position in the chromosome/scaffold, and
subsequent elements are optional which may refer to the allele identifier and other identifiers.
\item $list_pheno:
\itemize{
\item $y: named vector of numeric phenotype data
\item $pop: population or groupings corresponding to each element of y
\item $trait_name: name of the trait
}
\item $COVAR: numeric n samples x k covariates matrix with non-null row and column names.
}}

\item{vec_idx_training}{vector of numeric indexes referring to the training set}

\item{vec_idx_validation}{vector of numeric indexes referring to the validation set}

\item{other_params}{list of additional parameters,
\itemize{
\item other_params$nIter: total number of iterations (Default=12e3)
\item other_params$burnIn: number of burn-in iterations (Default=2e3)
\item other_params$out_prefix: prefix of temporary output files (Default="bglr_bayesB-")
}}

\item{verbose}{show Bayes B regression messages? (Default=FALSE)}
}
\value{
\itemize{
\item Ok:
\itemize{
\item $list_perf:
\itemize{
\item $mbe: mean bias error
\item $mae: mean absolute error
\item $rmse: root mean squared error
\item $r2: coefficient of determination
\item $corr: Pearson's product moment correlation
\item $power_t10: fraction of observed top 10 phenotype values correctly predicted
\item $power_b10: fraction of observed bottom 10 phenotype values correctly predicted
\item $var_additive: variance of predicted phenotype values (estimator of additive genetic variance)
\item $var_residual: variance of differnce between observed and predicted phenotype values (estimator of residual variance)
\item $h2: narrow-sense heritability estimate
}
\item $df_y_validation:
\itemize{
\item $id: names of the samples/entries/pools,
\item $pop: population from which the sample/entry/pool belongs to
\item $y_true: observed phenotype values
\item $y_pred: predicted phenotype values
}
\item $vec_effects: named numeric vector of estimated effects, where the names correspond to the
SNP/allele identity including chromosome/scaffold, position and optionally allele.
\item $n_non_zero: number of non-zero estimated effects (effects greater than machine epsilon ~2.2e-16)
}
\item Err: gpError
}
}
\description{
Bayes B regression model
(scaled t-distributed effects with probability \eqn{\pi}; and zero effects with probability \eqn{1-\pi},
where \eqn{\pi \sim \beta(\theta_1, \theta_2)})
}
\examples{
list_sim = fn_simulate_data(verbose=TRUE)
G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, verbose=TRUE)
n = nrow(list_merged$G)
vec_idx_training = sample(c(1:n), floor(n/2))
vec_idx_validation = c(1:n)[!(c(1:n) \%in\% vec_idx_training)]
list_ols = fn_ols(list_merged, vec_idx_training, vec_idx_validation, verbose=TRUE)
}

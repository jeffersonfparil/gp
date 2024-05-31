% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{fn_prediction_performance_metrics}
\alias{fn_prediction_performance_metrics}
\title{Genomic prediction metrics
NOTES:
\itemize{
\item Add more metrics which we may be interested in. while keeping the function signatures, i.e. input and output consistent
\item The number of non-zero estimated marker effects are extracted within each model, i.e. in \code{models.R}
\item Run time of each genomic prediction and validation on a single set is measured during each cross-validation replicate, i.e in \code{cross-validation.R}
}}
\usage{
fn_prediction_performance_metrics(y_true, y_pred, verbose = FALSE)
}
\arguments{
\item{y_true}{numeric vector observed phenotype values}

\item{y_pred}{numeric vector predicted phenotype values}
}
\value{
$mbe: mean bias error
$mae: mean absolute error
$rmse: root mean squared error
$r2: coefficient of determination
$corr: Pearson's product moment correlation
$power_t10: fraction of observed top 10 phenotype values correctly predicted
$power_b10: fraction of observed bottom 10 phenotype values correctly predicted
$h2: estimated narrow-sense heritability weighted by Pearson's correlation between observed and predicted phenotype values
Formula:
- h2_predicted = corr(observed, predicted) * h2_true
- Caveat: yields NA if the estimated heritability lies beyond 0 and 1
which happens when correlation between observed and predicted phenotypes is negative, and
the predicted phenotypes have higher variance than observed phenotypes.
This unexpected higher variance in the predictions implies that the genomic prediction model performs very poorly,
and therefore heritability estimates from them should not be considered.
Assumptions:
- correlation between predicted and true trait values relates the variance of the predicted trait values with the true heritability,
- predicted trait values are due to purely additive effects,
- variance in the predicted trait values is proportional to the additive genetic variance (~Va)
- variance in the true trait values represent the total phenotypic variance (Vp)
}
\description{
Genomic prediction metrics
NOTES:
\itemize{
\item Add more metrics which we may be interested in. while keeping the function signatures, i.e. input and output consistent
\item The number of non-zero estimated marker effects are extracted within each model, i.e. in \code{models.R}
\item Run time of each genomic prediction and validation on a single set is measured during each cross-validation replicate, i.e in \code{cross-validation.R}
}
}
\examples{
y_pred = stats::rnorm(100)
y_true = y_pred + stats::rnorm(100)
list_metrics = fn_prediction_performance_metrics(y_true=y_true, y_pred=y_pred, verbose=TRUE)
}
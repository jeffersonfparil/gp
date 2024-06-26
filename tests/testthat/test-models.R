# library(testthat)
# source("R/io.R")
# source("R/metrics.R")
# source("R/models.R")

test_that("fn_ols", {
    set.seed(123)
    list_sim = fn_simulate_data(verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, verbose=TRUE)
    n = nrow(list_merged$G)
    vec_idx_training = sample(c(1:n), floor(n/2))
    vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
    list_ols = fn_ols(list_merged, vec_idx_training, vec_idx_validation, verbose=TRUE)
    expect_equal(list_ols$list_perf$corr < 0.9, TRUE)
    expect_equal(list_ols$list_perf$corr, cor(list_ols$df_y_validation$y_true, list_ols$df_y_validation$y_pred))
    expect_equal(nrow(list_ols$df_y_validation), length(vec_idx_validation))
    expect_equal(list_ols$n_non_zero, ncol(G))
})

test_that("fn_ridge", {
    set.seed(123)
    list_sim = fn_simulate_data(verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, verbose=TRUE)
    n = nrow(list_merged$G)
    vec_idx_training = sample(c(1:n), floor(n/2))
    vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
    list_ridge = fn_ridge(list_merged, vec_idx_training, vec_idx_validation, verbose=TRUE)
    expect_equal(list_ridge$list_perf$corr < 0.9, TRUE)
    expect_equal(list_ridge$list_perf$corr, cor(list_ridge$df_y_validation$y_true, list_ridge$df_y_validation$y_pred))
    expect_equal(nrow(list_ridge$df_y_validation), length(vec_idx_validation))
    expect_equal(sum(list_ridge$vec_effects != 0.0), 1+ncol(G))
    expect_equal(list_ridge$n_non_zero, 1+ncol(G))
})

test_that("fn_lasso", {
    set.seed(123)
    list_sim = fn_simulate_data(verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, verbose=TRUE)
    n = nrow(list_merged$G)
    vec_idx_training = sample(c(1:n), floor(n/2))
    vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
    list_lasso = fn_lasso(list_merged, vec_idx_training, vec_idx_validation, verbose=TRUE)
    expect_equal(list_lasso$list_perf$corr < 0.9, TRUE)
    expect_equal(list_lasso$list_perf$corr, cor(list_lasso$df_y_validation$y_true, list_lasso$df_y_validation$y_pred))
    expect_equal(nrow(list_lasso$df_y_validation), length(vec_idx_validation))
    expect_equal(sum(list_lasso$vec_effects != 0.0) < 1+ncol(G), TRUE)
    expect_equal(list_lasso$n_non_zero < ncol(G), TRUE)
})

test_that("fn_elastic_net", {
    set.seed(123)
    list_sim = fn_simulate_data(verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, verbose=TRUE)
    n = nrow(list_merged$G)
    vec_idx_training = sample(c(1:n), floor(n/2))
    vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
    list_elastic_net = fn_elastic_net(list_merged, vec_idx_training, vec_idx_validation, verbose=TRUE)
    expect_equal(list_elastic_net$list_perf$corr < 0.9, TRUE)
    expect_equal(list_elastic_net$list_perf$corr, cor(list_elastic_net$df_y_validation$y_true, list_elastic_net$df_y_validation$y_pred))
    expect_equal(nrow(list_elastic_net$df_y_validation), length(vec_idx_validation))
    expect_equal(sum(list_elastic_net$vec_effects != 0.0) < 1+ncol(G), TRUE)
    expect_equal(list_elastic_net$n_non_zero < ncol(G), TRUE)
})

test_that("fn_Bayes_A", {
    set.seed(123)
    list_sim = fn_simulate_data(verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, verbose=TRUE)
    n = nrow(list_merged$G)
    vec_idx_training = sample(c(1:n), floor(n/2))
    vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
    list_Bayes_A = fn_Bayes_A(list_merged, vec_idx_training, vec_idx_validation, verbose=TRUE)
    expect_equal(list_Bayes_A$list_perf$corr < 0.9, TRUE)
    expect_equal(list_Bayes_A$list_perf$corr, cor(list_Bayes_A$df_y_validation$y_true, list_Bayes_A$df_y_validation$y_pred))
    expect_equal(nrow(list_Bayes_A$df_y_validation), length(vec_idx_validation))
    expect_equal(sum(list_Bayes_A$vec_effects != 0.0), 1+ncol(G))
    expect_equal(list_Bayes_A$n_non_zero, 1+ncol(G))
})

test_that("fn_Bayes_B", {
    set.seed(123)
    list_sim = fn_simulate_data(verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, verbose=TRUE)
    n = nrow(list_merged$G)
    vec_idx_training = sample(c(1:n), floor(n/2))
    vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
    list_Bayes_B = fn_Bayes_B(list_merged, vec_idx_training, vec_idx_validation, verbose=TRUE)
    expect_equal(list_Bayes_B$list_perf$corr < 0.9, TRUE)
    expect_equal(list_Bayes_B$list_perf$corr, cor(list_Bayes_B$df_y_validation$y_true, list_Bayes_B$df_y_validation$y_pred))
    expect_equal(nrow(list_Bayes_B$df_y_validation), length(vec_idx_validation))
    expect_equal(sum(list_Bayes_B$vec_effects != 0.0), 1+ncol(G))
    expect_equal(list_Bayes_B$n_non_zero, 1+ncol(G))
})

test_that("fn_Bayes_C", {
    set.seed(123)
    list_sim = fn_simulate_data(verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, verbose=TRUE)
    n = nrow(list_merged$G)
    vec_idx_training = sample(c(1:n), floor(n/2))
    vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
    list_Bayes_C = fn_Bayes_C(list_merged, vec_idx_training, vec_idx_validation, verbose=TRUE)
    expect_equal(list_Bayes_C$list_perf$corr < 0.9, TRUE)
    expect_equal(list_Bayes_C$list_perf$corr, cor(list_Bayes_C$df_y_validation$y_true, list_Bayes_C$df_y_validation$y_pred))
    expect_equal(nrow(list_Bayes_C$df_y_validation), length(vec_idx_validation))
    expect_equal(sum(list_Bayes_C$vec_effects != 0.0), 1+ncol(G))
    expect_equal(list_Bayes_C$n_non_zero, 1+ncol(G))
})

test_that("fn_gBLUP", {
    set.seed(123)
    list_sim = fn_simulate_data(verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, verbose=TRUE)
    n = nrow(list_merged$G)
    vec_idx_training = sample(c(1:n), floor(n/2))
    vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
    list_gBLUP = fn_gBLUP(list_merged, vec_idx_training, vec_idx_validation, verbose=TRUE)
    expect_equal(list_gBLUP$list_perf$corr < 0.9, TRUE)
    expect_equal(list_gBLUP$list_perf$corr, cor(list_gBLUP$df_y_validation$y_true, list_gBLUP$df_y_validation$y_pred))
    expect_equal(nrow(list_gBLUP$df_y_validation), length(vec_idx_validation))
    expect_equal(sum(list_gBLUP$vec_effects != 0.0), 1+nrow(G))
    expect_equal(list_gBLUP$n_non_zero, 1+nrow(G))
})

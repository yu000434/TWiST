library(TWiST)
load("example_data/example_data_stage1training.rda")

# Recover the individual-level genotypes from the example cell data.
ids <- unique(rownames(geno.cell))
donor <- match(rownames(geno.cell), ids)
geno <- geno.cell[match(ids, rownames(geno.cell)), , drop=FALSE]

set.seed(1)
twist <- twist_train_model(y=gene_exp_i, geno_cell=geno.cell, pt=pt,
    libsize=libsize, covar=covar, nlambda=10)

set.seed(1)
fast <- fast_twist_train_model(y=gene_exp_i, geno=geno, donor=donor, pt=pt,
    libsize=libsize, covar=covar, nlambda=10)

set.seed(1)
fpca <- fpca_train_model(y=gene_exp_i, geno=geno, donor=donor, pt=pt,
    libsize=libsize, covar=covar)

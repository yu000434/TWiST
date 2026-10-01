# TWiST

TWiST (**TW**AS **i**n p**S**eudo**T**ime) is an R package for single-cell TWAS analysis of heterogeneous cell types, where gene expression and eQTL effects can vary along a continuous cell state within the cell type.
Cell state is defined by pseudotime. This package implements two main analyses:

**Stage 1:** Train models to predict genetically regulated gene expression from cis-SNPs using single-cell eQTL data. The function `twist_train_model()` uses fast-TWiST by default, with TWiST and FPCA available through the `method` argument.

**Stage 2:** Conduct association analysis between gene expression and a trait using GWAS summary statistics.

<img src="example_data/overview.png" alt="overview" width="800"/>

**Pre-trained TWiST models are provided for CD4+ T cells, CD8+ T cells, and B cells. Users interested in these cell types may use the pre-trained models and skip Stage 1 (see example below).**

## 1. Installation

Install `TWiST` with fast-TWiST and FPCA from GitHub:

```r
devtools::install_github("yu000434/TWiST", ref="fast-twist-fpca")
```

In addition, install the `plink2R` package to read genotype data (PLINK files) into R:

```r
devtools::install_github("gabraham/plink2R/plink2R")
```

## 2. Pre-trained models

Download pre-trained models for three immune cell types from [pretrained_models](pretrained_models): `twist_weights_T_CD4.rda` (CD4+ T cells), `twist_weights_T_CD8.rda` (CD8+ T cells), and `twist_weights_B.rda` (B cells). These models were trained on OneK1K data using the original TWiST method.

Each `.rda` file includes three objects:

* `wgtlist`: Information on genes for which the model has been trained. A data frame of five columns:
    * `ID`: Gene ID
    * `CHR`: Chromosome
    * `P0`: Gene start
    * `P1`: Gene end
    * `tss`: Transcription start site
* `weights_pred`: Pre-trained prediction models. Each entry is a gene, in the same order as `wgtlist$ID`. Each model contains:
    * `Wmat`: Coefficients of the prediction model. A matrix of (number of SNPs) x (number of B-spline bases).
    * `knots`: Internal knots for the B-spline functions used to model SNP effects on gene expression. Boundary knots 0 and 1 are not included.
    * `degree`: Degree of B-spline basis functions.
    * `n`: Number of cells in the eQTL data used for model training.
* `bim_train`: Information on model SNPs in the format of a PLINK bim file:
    * `CHR`: Chromosome
    * `SNP`: SNP ID
    * `cM`: SNP position in centimorgan
    * `BP`: SNP position in base pairs
    * `A1`: Effect allele. Coefficients in `weights_pred` are with respect to A1.
    * `A2`: Other allele

## 3. Example: Association analysis

To run this example, download the repository and run the code from its main directory. Additional datasets are provided in [example_data](example_data):

* `RA_sumstats_chr6.txt`: GWAS summary statistics for rheumatoid arthritis, chromosome 6 (Ishigaki et al., Nature Genetics 2022).
* `1000G.EUR.6.{bed,bim,fam}`: 1000 Genomes European genotype data, chromosome 6. Download all chromosomes from [here](https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2).

First, load the required packages:

```r
library(dplyr)
library(plink2R)
library(TWiST)
```

Read GWAS summary statistics into R. Compute the effective sample size defined as `ncases*ncontrols/(ncases+ncontrols)`:

```r
sumstats <- data.table::fread("example_data/RA_sumstats_chr6.txt")
ngwas <- 22350*74823/(22350+74823)
```

Load pre-trained models and subset to chromosome 6, using CD8+ T cells as an example:

```r
ctype <- "T_CD8"
load(paste0("pretrained_models/twist_weights_",ctype,".rda"))
wgtlist.chr <- wgtlist %>% filter(CHR==6)
weights_pred.chr <- weights_pred[wgtlist.chr$ID]
```

Load reference genotype data from the 1000 Genomes European sample:

```r
genos.chr <- read_plink("example_data/1000G.EUR.6")
```

Run TWiST association analysis:

```r
res <- twist_association(
    sumstat=sumstats, wgtlist=wgtlist.chr, weights_pred=weights_pred.chr,
    bim_train=bim_train, genos=genos.chr, ngwas=ngwas)
```

View results. Type `?twist_association` for definitions of the outputs:

```r
names(res)
str(res$out.tbl)
```

Create QQ plots for the global, dynamic, and nonlinear tests:

```r
library(qqman)
par(mfrow=c(1,3))
qq(res$out.tbl$p.global, main="Global test", ylim=c(0,220))
qq(res$out.tbl$p.dynamic, main="Dynamic test", ylim=c(0,220))
qq(res$out.tbl$p.nonlinear, main="Nonlinear test", ylim=c(0,220))
```

<img src="example_data/QQ_T_CD8_chr6.png" alt="QQ" width="800"/>

## 4. Example: Training prediction models

If you have your own single-cell eQTL data and would like to train a prediction model, below is an example using simulated data.

Load the example dataset:

```r
library(TWiST)
load("example_data/example_data_stage1training.rda")
```

This dataset includes real genotype data and genotype principal components (PCs) from the European subset of 1000 Genomes (chromosome 6), with simulated expression counts, pseudotime, library sizes, and covariates:

* `gene_exp_i`: Expression counts for one gene.
* `geno.cell`: Genotypes of cis-SNPs, repeated at the cell level. Cells from the same individual have the same genotypes, and row names identify the individual.
* `pt`: Pseudotime of the cells.
* `libsize`: Library size.
* `covar`: Age, sex, 10 expression PCs (`PC_1` to `PC_10`), and 10 genotype PCs (`genPC_1` to `genPC_10`).

Train a prediction model for one gene. The default method is fast-TWiST. The `donor` vector identifies the individual corresponding to each cell:

```r
donor <- rownames(geno.cell)
set.seed(1)
model <- twist_train_model(y=gene_exp_i, geno_cell=geno.cell, donor=donor,
    pt=pt, knots=c(0.25,0.5,0.75), degree=3, nlambda=10,
    libsize=libsize, covar=covar)
```

Use `method="twist"` for the original TWiST implementation or `method="fpca"` for FPCA:

```r
set.seed(1)
model.twist <- twist_train_model(y=gene_exp_i, geno_cell=geno.cell,
    pt=pt, nlambda=10, libsize=libsize, covar=covar, method="twist")

set.seed(1)
model.fpca <- twist_train_model(y=gene_exp_i, geno_cell=geno.cell, donor=donor,
    pt=pt, libsize=libsize, covar=covar, method="fpca")
```

FPCA uses 50 pseudotime bins and retains components explaining 99% of trajectory variation by default. These settings can be changed through `nbins` and `pve`. See `?twist_train_model` for parameters.

For each method, collect the complete fitted models across genes in `weights_pred`, named by gene ID, and prepare `wgtlist` and `bim_train` as described in Section 2. Pass these objects to `twist_association()` as in Section 3. For FPCA, `Wmat` contains SNP-by-component coefficients; keep the eigenfunctions and other model entries with the weights.

Global and dynamic P values from fast-TWiST and FPCA can be combined separately using `cauchy_combine()`. See `?cauchy_combine` for its arguments.

Code for simulating this dataset is provided [here](example_data/simulate_example_training.R).

## 5. Reference

Qi G, Lila E, Ji Z, Shojaie A, Battle A, Sun W. Transcriptome-wide association studies at cell state level using single-cell eQTL data. *Cell Genomics* (2026). https://www.cell.com/cell-genomics/fulltext/S2666-979X(25)00316-7.

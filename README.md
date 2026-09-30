# TWiST

TWiST (**TW**AS **i**n p**S**eudo**T**ime) is an R package for single-cell TWAS analysis of heterogeneous cell types, where gene expression and eQTL effects can vary along a continuous cell state within the cell type. Cell state is defined by pseudotime. This package implements two main analyses:

**Stage 1:** Train models to predict genetically regulated gene expression from cis-SNPs using `twist_train_model()`. The default method is fast-TWiST; TWiST and FPCA are also available.

**Stage 2:** Test the association between genetically regulated gene expression and a trait using GWAS summary statistics.

<img src="example_data/overview.png" alt="TWiST overview" width="800"/>

## 1. Installation

Install `TWiST` with fast-TWiST and FPCA from GitHub:

```r
devtools::install_github("yu000434/TWiST", ref="fast-twist-fpca")
```

In addition, install `plink2R` to read genotype data in PLINK format:

```r
devtools::install_github("gabraham/plink2R/plink2R")
```

## 2. Example data

The example below trains all three methods for one gene, HLA-A, runs association analysis, and combines the fast-TWiST and FPCA results.

Download this repository and run the example from its main directory. The required files are provided in [example_data](example_data):

* `example_data_stage1training.rda`: Real chromosome 6 genotypes from the European subset of 1000 Genomes, with simulated expression counts, pseudotime, library sizes, and covariates.
* `RA_sumstats_chr6.txt`: GWAS summary statistics for rheumatoid arthritis, chromosome 6 (Ishigaki et al., Nature Genetics 2022).
* `1000G.EUR.6.{bed,bim,fam}`: 1000 Genomes European genotype data used as the LD reference.

Load the required packages and training data:

```r
library(TWiST)
library(plink2R)
load("example_data/example_data_stage1training.rda")
```

The training dataset contains:

* `gene_exp_i`: Expression counts for one gene in 5,000 cells.
* `geno.cell`: Cell-level genotypes, with cells in rows and cis-SNPs in columns. Row names identify the individual for each cell.
* `pt`: Pseudotime of the cells, between 0 and 1.
* `libsize`: Library size of each cell.
* `covar`: Age, sex, 10 expression PCs, and 10 genotype PCs.

Prepare a genotype matrix with one row per individual and a vector giving the genotype row for each cell. These inputs can be used with all three methods:

```r
ids <- unique(rownames(geno.cell))
donor <- match(rownames(geno.cell), ids)
geno <- geno.cell[match(ids, rownames(geno.cell)), , drop=FALSE]
```

Read the GWAS summary statistics and LD reference. For this case-control study, the effective sample size is `ncases*ncontrols/(ncases+ncontrols)`:

```r
sumstats <- data.table::fread("example_data/RA_sumstats_chr6.txt")
ngwas <- 22350*74823/(22350+74823)
genos.chr <- read_plink("example_data/1000G.EUR.6")
```

## 3. Training prediction models

Train a prediction model using `twist_train_model()`. By default, fast-TWiST fits the spline-based Poisson model using individual-level summaries:

```r
set.seed(1)
model.fast <- twist_train_model(y=gene_exp_i, geno=geno, donor=donor,
    pt=pt, libsize=libsize, covar=covar, nlambda=10)
```

Use `method="twist"` for the original TWiST implementation or `method="fpca"` for FPCA:

```r
set.seed(1)
model.twist <- twist_train_model(y=gene_exp_i, geno=geno, donor=donor,
    pt=pt, libsize=libsize, covar=covar, nlambda=10, method="twist")

set.seed(1)
model.fpca <- twist_train_model(y=gene_exp_i, geno=geno, donor=donor,
    pt=pt, libsize=libsize, covar=covar, method="fpca")
```

FPCA estimates individual trajectories using Poisson mixed models in 50 pseudotime bins and retains components explaining 99% of trajectory variation by default. Component scores are predicted from cis-SNP genotypes. The number of bins and proportion of variation can be set through `nbins` and `pve`. See `?twist_train_model` for parameters.

Each model contains `Wmat`: SNP-by-spline coefficients for TWiST and fast-TWiST, or SNP-by-component coefficients for FPCA. Keep the complete model objects, including the basis information needed for association analysis.

## 4. Association analysis

Prepare the gene annotation and collect each method's fitted model in a list named by gene ID:

```r
wgtlist <- data.frame(ID="ENSG00000206503", CHR=6, P0=29909037,
                      P1=29913661, tss=29909037)
weights_pred.twist <- setNames(list(model.twist), wgtlist$ID)
weights_pred.fast <- setNames(list(model.fast), wgtlist$ID)
weights_pred.fpca <- setNames(list(model.fpca), wgtlist$ID)
```

Read the SNP information from the training genotype data. In this example, the training genotypes and LD reference come from the same 1000 Genomes files:

```r
bim_train <- read.table("example_data/1000G.EUR.6.bim",
    col.names=c("CHR", "SNP", "cM", "BP", "A1", "A2"))
```

Run association analysis separately for each method, using the GWAS data and LD reference loaded in Section 2:

```r
res.twist <- twist_association(
    sumstat=sumstats, wgtlist=wgtlist, weights_pred=weights_pred.twist,
    bim_train=bim_train, genos=genos.chr, ngwas=ngwas)

res.fast <- twist_association(
    sumstat=sumstats, wgtlist=wgtlist, weights_pred=weights_pred.fast,
    bim_train=bim_train, genos=genos.chr, ngwas=ngwas)

res.fpca <- twist_association(
    sumstat=sumstats, wgtlist=wgtlist, weights_pred=weights_pred.fpca,
    bim_train=bim_train, genos=genos.chr, ngwas=ngwas)
```

Each result contains an `out.tbl` table with gene information and P values for the global, dynamic, and nonlinear tests. View the results:

```r
results <- rbind(TWiST=res.twist$out.tbl,
                 fast.TWiST=res.fast$out.tbl, FPCA=res.fpca$out.tbl)
results[, c("ID", "p.global", "p.dynamic", "p.nonlinear")]
#                         ID     p.global    p.dynamic  p.nonlinear
# TWiST      ENSG00000206503 5.654862e-29 4.829743e-19 7.369468e-08
# fast.TWiST ENSG00000206503 4.343681e-31 9.052337e-21 5.296791e-15
# FPCA       ENSG00000206503 8.944971e-13 3.537372e-01 5.000000e-02
```

See `?twist_association` for the full output definition. For multiple genes, collect the fitted models in the same order as `wgtlist$ID` and provide their SNP information in `bim_train`.

## 5. Combining association results

The `res.fast` and `res.fpca` objects above contain the two methods' association results. Match genes by ID and combine global and dynamic P values separately using the Cauchy combination test:

```r
paired <- merge(res.fast$out.tbl, res.fpca$out.tbl,
                by="ID", suffixes=c(".fast", ".fpca"))

combined <- data.frame(
    ID=paired$ID,
    p.global=cauchy_combine(paired$p.global.fast, paired$p.global.fpca),
    p.dynamic=cauchy_combine(paired$p.dynamic.fast, paired$p.dynamic.fpca))
combined
#                ID     p.global    p.dynamic
# 1 ENSG00000206503 8.687362e-31 1.810467e-20
```

This step is optional. Nonlinear P values are not combined.

## Pre-trained OneK1K models

The [pretrained_models](pretrained_models) folder contains models trained on OneK1K data with the original TWiST method for CD4+ T cells (`twist_weights_T_CD4.rda`), CD8+ T cells (`twist_weights_T_CD8.rda`), and B cells (`twist_weights_B.rda`). These are separate from the simulated training example above.

Each file contains `wgtlist` (gene annotation), `weights_pred` (prediction models named by gene ID), and `bim_train` (training SNP information). Users of these models can skip training. The original association example, including QQ plots, is provided in [example.R](example_data/example.R).

## 6. Reference

Qi G, Lila E, Ji Z, Shojaie A, Battle A, Sun W. Transcriptome-wide association studies at cell state level using single-cell eQTL data. *Cell Genomics* (2026). https://www.cell.com/cell-genomics/fulltext/S2666-979X(25)00316-7.

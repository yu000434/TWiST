#' Train prediction model for genetically regulated gene expression
#' @description This function trains a model using single-cell expression eQTL data to predict genetically regulated gene expression (GReX) from cis-SNPs. Training uses fast-TWiST by default, with TWiST and FPCA as alternatives. This is the Stage-1 model of TWiST.
#' @param y Expression count data for one gene. Integer vector of length ncells (number of cells).
#' @param geno_cell Genotype of cis-SNPs of the gene. Matrix of dimensions ncells x nsnps (number of SNPs). Each row is a cell and each column is a SNP. Genotypes are repeated across cells from the same individual. Supply either \code{geno_cell} or individual-level \code{geno}. For fast-TWiST and FPCA, also supply \code{donor}.
#' @param pt Pseudotime of the cells. Pseudotime values learned from packages such as Slingshot or TSCAN should be scaled to rank/(number of cells) such that values are approximately uniformly distributed between 0 and 1.
#' @param knots Internal knots for B-spline basis functions in TWiST and fast-TWiST. Not including 0 and 1. Default value is \code{c(0.25,0.5,0.75)}.
#' @param degree Degree of B-spline basis functions in TWiST and fast-TWiST. Default value is 3 (cubic B-spline).
#' @param lambda Penalty parameter for TWiST and fast-TWiST. If NULL, selected by cross-validation.
#' @param nlambda Number of \code{lambda} values in cross-validation for TWiST and fast-TWiST. Default is 20.
#' @param libsize Library size of the cells. Vector of length ncells.
#' @param covar Covariates to be adjusted. Numerical matrix of ncells x (number of covariates). Not penalized.
#' @param method Training method: \code{"fast"} (default), \code{"twist"}, or \code{"fpca"}.
#' @param geno Individual-level cis-SNP genotypes, with individuals in rows and SNPs in columns. Supply \code{donor} to match cells to genotype rows.
#' @param donor Individual corresponding to each cell. With \code{geno}, these are row numbers in \code{geno}. With \code{geno_cell}, these are individual IDs; genotype rows are extracted in order of first appearance.
#' @param ... Additional arguments to \code{fast_twist_train_model} or \code{fpca_train_model}, such as \code{ncores} for fast-TWiST or \code{nbins} and \code{pve} for FPCA.
#'
#' @return A fitted prediction model containing:
#' \item{Wmat}{SNP-by-spline coefficients for TWiST and fast-TWiST, or SNP-by-component coefficients for FPCA.}
#' \item{knots, degree}{B-spline basis settings for TWiST and fast-TWiST.}
#' \item{efunctions, argvals}{Eigenfunctions and their pseudotime coordinates for FPCA.}
#' \item{n}{Number of cells used for training.}
#' With \code{method="twist"}, the model retains its \code{grpreg} class and entries, including \code{lambda.seq} and \code{cvm} when cross-validation is used. For other method-specific entries, see \code{fast_twist_train_model} and \code{fpca_train_model}. The complete model can be included in \code{weights_pred} for \code{twist_association}.
#'
#' @details
#' TWiST and fast-TWiST use a spline-based Poisson model. The unique molecular identifier (UMI) count that represents gene expression (for cell j in individual) is modeled using Poisson distribution with mean \eqn{\mu_{ij}}. The mean is further modeled using a log-linear model of genetically regulated expression (GRex) and covariates:
#' \deqn{\log(\mu_{ij})=c_0 + \log(\alpha_{ij}) + v_i(t_{ij}) +\gamma^T\boldsymbol{z}_{ij},}
#' where \eqn{c_0} is the intercept \eqn{\alpha_{ij}} is the library size, \eqn{v_i(t_{ij})} is the GReX, and \eqn{\boldsymbol{z}_{ij}} is a vector of the covariates. The GReX is modeled using a functional linear model of cis-SNPs:
#' \eqn{v_i(t)=\sum\limits_{k} w_k(t)g_{ik}}, where \eqn{g_{ik}} is SNP k of individual i, and \eqn{w_k(t)} is a B-spline function over pseudotime t.
#' A group elastic-net penalty is applied to select SNP groups and shrink their coefficients.
#'
#' The \code{"twist"} method uses the original TWiST implementation, fitting log library size and covariates jointly with the genetic effects. The \code{"fast"} and \code{"fpca"} methods first estimate covariate effects and hold their contribution fixed as an offset, including log library size.
#' Fast-TWiST uses individual-level summaries during fitting. FPCA estimates individual expression trajectories using binned Poisson mixed models, then predicts their component scores from cis-SNP genotypes. FPCA uses 50 bins and retains components explaining 99\% of trajectory variation by default.
#'
#' @import grpreg
#' @import splines
#' @import dplyr
#' @export
twist_train_model <- function(y, geno_cell=NULL, pt, knots=c(0.25,0.5,0.75), degree=3,
                              lambda=NULL, nlambda=20, libsize, covar,
                              method=c("fast", "twist", "fpca"), geno=NULL, donor, ...){
    method <- match.arg(method)
    if (method != "twist") {
        if (is.null(geno)) {
            ids <- unique(donor)
            geno <- geno_cell[match(ids, donor), , drop=FALSE]
            donor <- match(donor, ids)
        }
        if (method == "fast")
            return(fast_twist_train_model(y=y, geno=geno, donor=donor, pt=pt,
                libsize=libsize, covar=covar, knots=knots, degree=degree,
                lambda=lambda, nlambda=nlambda, ...))
        return(fpca_train_model(y=y, geno=geno, donor=donor, pt=pt,
            libsize=libsize, covar=covar, ...))
    }
    if (is.null(geno_cell)) geno_cell <- geno[donor, , drop=FALSE]

    # intercept needs to be set to TRUE to get the complete set of bases (the bases sum to 1 at each pseudotime point)
    # If intercept=FALSE, bs() will remove one basis
    ptbs <- bs(pt, knots=knots, degree=degree, intercept=TRUE, Boundary.knots=c(0,1))
    desmat <- generate_desmat(geno.temp=geno_cell, ptbs=ptbs)
    grp <- rep(1:ncol(geno_cell),each=ncol(ptbs))

    # Cross validation to select optimal lambda parameter
    if (is.null(lambda)){
        gr_cv <- cv.grpreg(X=cbind(log(libsize),covar,desmat), y=y,
                           group=c(rep(0,1+ncol(covar)), grp), penalty="grLasso",
                           family="poisson", nlambda=nlambda, alpha=0.5, nfolds=5)
        lambda.min <- gr_cv$lambda.min
    } else{
        lambda.min <- lambda
    }

    res.opt <- grpreg(X=cbind(log(libsize),covar,desmat), y=y,
                      group=c(rep(0,1+ncol(covar)), grp), penalty="grLasso",
                      family="poisson", lambda=lambda.min, alpha=0.5)
    # Different from glmnet, res.opt$beta from grpnet includes intercept
    res.opt$Wmat <- matrix(res.opt$beta[-(1:(ncol(covar)+2))], nrow=ncol(geno_cell), byrow=TRUE,
                              dimnames=list(colnames(geno_cell),paste0("bs",1:ncol(ptbs))))
    res.opt$knots <- knots
    res.opt$degree <- degree

    if (is.null(lambda)){
        res.opt$lambda.seq <- gr_cv$lambda
        res.opt$cvm <- gr_cv$cvm
    }

    return(res.opt)
}

#' Generate design matrix that combines genotypes and B-spline basis functions.
#' @param geno.temp Genotype of cis-SNPs. Matrix of ncells x (number of SNPs).
#' @param ptbs B-spline basis functions. Matrix of ncells x (number of B-spline basis functions)
generate_desmat <- function(geno.temp, ptbs){
    if (nrow(geno.temp)!=nrow(ptbs)) stop("geno.temp and ptbs should have the same number of rows.")

    desmat <- matrix(NA, nrow=nrow(geno.temp), ncol=ncol(ptbs)*ncol(geno.temp))
    for (i in 1:ncol(geno.temp)){
        desmat[,(ncol(ptbs)*(i-1)+1):(ncol(ptbs)*i)] <- ptbs*geno.temp[,i]
    }

    return(desmat)
}

# run a single regression for the gene data in gdat
# pdat: data frame with predictors
# calibrate: use Fredrickson et al.'s predictors
#' Fit a linear model to expression data for a single gene
#'
#' @param gdat Data frame containing expression data for the gene
#' @param pdat Data frame contining predictors
#' @param calibrate Use Fredrickson et al.'s predictors (\code{TRUE}) or Brown et al.'s (\code{FALSE})
#' @param joinby What field to join genetic data with psychometric data (leave as is)
#' @param resids Whether to use residualized genetic data or raw data
#' @return A data frame with field \code{Pred} naming the predictor, and \code{x} containing its value.
#' @export 
runLM <- function(gdat, pdat, calibrate=TRUE, joinby="SID", resids=FALSE) {
    cols.keep <- c("Hed", "Eud")
    if (calibrate) {
        cols.keep <- c("ZScore_Hed", "ZScore_Eud")
    } else {}
    dat <- inner_join(gdat, pdat, joinby)
    if (!resids) {
        mod.func <- as.formula(
            paste("Value ~ Male + Age + White + BMI + Alcohol + Smoke + ",
                  "Illness + CD3D + CD3E + CD4 + CD8A + CD19 + FCGR3A + ",
                  "NCAM1 + CD14 + ",
                  paste(cols.keep, collapse=" + "), sep=""))
    } else {
        mod.func <- as.formula(paste0("Value ~ ", paste(cols.keep, collapse=" + ")))
    }
    reg <- lm(mod.func, dat)
    data.frame(Pred=cols.keep, x=coef(reg)[cols.keep])
}

#' Run 53 repeated regressions as in Fredrickson et al.
#'
#' @param gedata Genetic data
#' @param preds Psychometric data (predictors)
#' @param geneSigns Dataframe with contrast codes for genes
#' @param calibrate Use eudaimonic/hedonic z-scores from Fredrickson et al. (default \code{TRUE}). Set to FALSE to use Brown et al.'s predictors.
#' @param joinby On which field to join the Psychometric / Genetic data (this should be left as is)
#' @param resids Whether gedata is the residualized gene data (\code{TRUE}) or on the raw data (\code{FALSE}=default).
#' @return A dataframe containing the estimated beta weights for the eudaimonic factor.
#' @export
rr53 <- function(gedata, preds, geneSigns, calibrate=TRUE, joinby="SID",
                 resids=FALSE) {
    lmResults <- gedata %>% group_by(Gene) %>%
        do(runLM(., preds, calibrate, joinby, resids)) %>%
        ungroup
    if (!resids) { # now modulate by the sign of the gene
        # unless we are dealing with the residuals
        # in which case contrast codes have already been applied,
        # at the very beginning.
        res <- inner_join(geneSigns, lmResults, by="Gene") %>%
            mutate(Y=x*Sign) %>% select(Gene, Pred, Y)
    } else {
        res <- lmResults %>% select(Gene, Pred, Y=x) %>% as.data.frame
    }
    return(res)
}

#' Perform a one sample t-test on the beta values
#'
#' @param x data frame, with beta values for one of the two psychometric factors in variable named \code{Y}
#' @param nbs number of bootstrap samples
#' @return data frame with results of one-sample test; \code{t} and \code{p} are the results from a parametric analysis, while \code{t.b} and \code{p.b} are the results using bootstrapped standard errors.
#' @export
oneSampTest <- function(x, nbs=10000) {
    calct <- function(xx) {
        diff.se <- sd(xx)/sqrt(length(xx))
        mean(xx)/diff.se
    }
    t.val <- calct(x$Y)
    degfreedom <- length(x$Y)-1
    p.val <- 2*pt(abs(t.val), degfreedom, lower.tail=FALSE)
    bootdist <- replicate(nbs, 
                          mean(sample(x$Y, length(x$Y), replace=TRUE)))
    se.boot <- sd(bootdist)
    t.boot <- mean(x$Y)/se.boot
    p.boot <- 2*pt(abs(t.boot), degfreedom, lower.tail=FALSE)
    data.frame(x=mean(x$Y), t=t.val, p=p.val, t.b=t.boot, p.b=p.boot)
}

#######################
# for the test of gene expression independence

#' Compute the mean absolute correlation for a correlation matrix.
#'
#' This function computes the correlation matrix using \code{\link{cor}},
#' and then calculates the mean absolute pairwise correlation in the
#' upper triangle of the matrix.
#' @param x Matrix of data
#' @return Mean absolute correlation for the matrix.
#' @export
meanAbsCor <- function(x) {
    cmx <- cor(x)    
    mean(abs(cmx[upper.tri(cmx)]))
}

#' Compute MAC after randomly shuffling columns
#'
#' Randomly shuffle the columns of the matrix and then recompute the
#' mean absolute correlation using \code{\link{meanAbsCor}}.
#' @param x Matrix of gene expression data (rows=subjects, columns=genes)
#' @return Mean absolute correlation for the shuffled matrix
#' @export
shuffleAndRecomputeMAC <- function(x) {
    smx <- apply(x, 1, sample) # shuffle genes for each subject
    meanAbsCor(smx)
}

#######################################
# for reshuffling residuals

#' Compute the residual gene expression values
#'
#' Calculates residual gene expression values for a single gene after
#' controlling for the 15 control factors in Fredrickson et al.
#'
#' @param x Data frame containing gene values for a single gene along with predictor values
#' @return Data frame with new residualized gene values.
#' @export
getGeneResids <- function(x) {
    
    # first we need to apply the contrast codes to the gene values
    # beforehand so that the direction of the prediction is not lost.
    x.lm <- lm(Value ~ Male + Age + White + BMI + Alcohol + Smoke +
                   Illness + CD3D + CD3E + CD4 + CD8A + CD19 + FCGR3A +
                       NCAM1 + CD14, x)
    res <- x$Value
    cc <- x %>% select(Value, Male, Age, White, BMI, Alcohol, Smoke,
                       Illness, CD3D, CD3E, CD4, CD8A, CD19, FCGR3A,
                       NCAM1, CD14) %>% complete.cases
    res[cc] <- residuals(x.lm)
    res[!cc] <- NA
    data.frame(SID=x$SID, Gene=x$Gene, Value=res)
}

#' Run the RR53 analysis on the gene expression residuals
#'
#' @param gdat Data frame with residualized gene expression data
#' @param pdat Data frame with predictors
#' @param gSigns Data frame with gene contrast codes
#' @return t-value for the one-sample test (bootstrapped SEs)
#' @export
rr53Resids <- function(gdat, pdat, gSigns) {
    res <- rr53(gdat, pdat, gSigns, resids=TRUE)
    tsum <- res %>% group_by(Pred) %>% do(oneSampTest(.))
    tsum$t.b[1]
}

#' Shuffle subject identifiers
#'
#' @param x A data frame with column \code{SID} containing subject identifiers.
#' @return Same data frame with shuffled IDs.
#' @export
shuffleSubjects <- function(x) {
    x$SID <- sample(x$SID)
    x
}

#' Generate variance-covariance matrix for gene expression data
#'
#' @param dat matrix of raw gene expression data
#' @return Estimated variance-covariance matrix
#' @export
generateVCMatrix <- function(dat) {
    # derive the variance-covariance matrix of the
    # data from the raw data (dat)
    gmx.var <- apply(dat, 2, var) # variance for each unit
    cmx <- cor(dat) # correlations between units
    # now create the full variance covariance matrix
    # using some matrix multiplication tricks
    omx <- matrix(1, nrow=dim(cmx)[1], ncol=dim(cmx)[2]) # matrix of ones
    cov1.sd <- t(sapply(sqrt(gmx.var), rep, times=ncol(cmx)))
    cov1.sd.t <- t(cov1.sd)
    cov1.mx <- omx
    cov1.mx[upper.tri(cov1.mx)] <- cov1.sd[upper.tri(cov1.sd)]
    cov1.mx[lower.tri(cov1.mx)] <- cov1.sd.t[lower.tri(cov1.sd.t)]
    cov2.sd <- sapply(sqrt(gmx.var), rep, times=nrow(cmx))
    cov2.sd.t <- t(cov2.sd)
    cov2.mx <- omx
    cov2.mx[upper.tri(cov2.mx)] <- cov2.sd[upper.tri(cov2.sd)]
    cov2.mx[lower.tri(cov2.mx)] <- cov2.sd.t[lower.tri(cov2.sd.t)]
    var.mx <- omx
    diag(var.mx) <- gmx.var
    # this is the full variance/covariance matrix
    cmx * var.mx * cov1.mx * cov2.mx
}

#' Randomly generate psychometric/demographic data
#'
#' Basically wraps a call to \code{\link{mvrnorm}}, with
#' some additional processing for categorical variables
#' @param n number of subjects (rows)
#' @param pmeans named vector of predictor means
#' @param pvcov variance-covariance matrix for predictors
#' @return a dataframe of new predictors
#' @importFrom MASS mvrnorm
#' @export
simulatePsychometricData <- function(n, pmeans, pvcov) {
    zeroOrOne <- function(x) {
        1*(x>=.5)
    }
    mx <- MASS::mvrnorm(n, pmeans, pvcov)
    colnames(mx) <- names(pmeans)
    dmx <- as.data.frame(mx)
    dmx$Male <- zeroOrOne(dmx$Male)
    dmx$Age <- round(dmx$Age)
    dmx$White <- zeroOrOne(dmx$White)
    dmx$Alcohol <- zeroOrOne(dmx$Alcohol)
    dmx$Smoke <- zeroOrOne(dmx$Smoke)
    dmx$Illness <- ifelse(dmx$Illness<0, 0, round(dmx$Illness))
    dmx$SID <- 1:nrow(dmx)
    return(dmx[,c("SID", setdiff(colnames(dmx),"SID"))])
}

# not exported
#' @importFrom MASS mvrnorm
simulateGeneMatrix <- function(n, gmeans, gvcov, indep=FALSE) {
    if (indep) {
        gvcov[upper.tri(gvcov)] <- 0
        gvcov[lower.tri(gvcov)] <- 0
    } else {}
    mvrnorm(n, gmeans, gvcov)
}

#' Generate simulated gene expresion data
#'
#' Simply wraps a call to \code{\link{mvrnorm}} in package \code{MASS}.  If \code{indep} is \code{TRUE}, zeroes out all of the covariances in the matrix.
#' @param n Number of rows (subjects) of data to generate
#' @param gmeans Mean expression for each gene
#' @param gvcov
#' @param indep whether to generate dependent (default=\code{FALSE}) or independent (\code{TRUE}) data.
#' @return A dataframe with SID (subject ID), Gene and Value as columns.
#' @seealso \code{\link{generateVCMatrix}}
#' @export
simulateGeneData <- function(n, gmeans, gvcov, indep=FALSE) {
   gmxt <- simulateGeneMatrix(n, gmeans, gvcov, indep)
   colnames(gmxt) <- names(gmeans)
   gdf <- gmxt %>% as.data.frame %>% mutate(SID=row_number()) %>%
       gather("Gene", "Value", -SID) %>% arrange(SID)
   gdf
}

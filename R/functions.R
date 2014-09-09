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
#' @param resids Whether to run the RR53 on the residualized gene data (\code{TRUE}) or on the raw data (\code{FALSE}=default).
#' @return A dataframe containing the estimated beta weights for the eudaimonic factor.
#' @export
rr53 <- function(gedata, preds, geneSigns, calibrate=TRUE, joinby="SID",
                 resids=FALSE) {
    gedata %>% group_by(Gene) %>%
        do(runLM(., preds, calibrate, joinby, resids)) %>%
        ungroup %>%
        inner_join(geneSigns, lmResults, by="Gene") %>%
        mutate(Y=x*Sign) %>% select(Gene, Pred, Y)
}

#' Perform a one sample t-test on the beta values
#'
#' @param x data frame, with beta values for one of the two psychometric factors in variable named \code{Y}
#' @return data frame with results of one-sample test.
#' @export
oneSampTest <- function(x) {
    diff.se <- sd(x$Y)/sqrt(length(x$Y))
    t.val <- mean(x$Y)/diff.se
    degfreedom <- length(x$Y)-1
    p.val <- 2*pt(abs(t.val), degfreedom, lower.tail=FALSE)
    data.frame(x=mean(x$Y), t=t.val, p=p.val)
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
#' @return t-value for the one-sample test
#' @export
rr53Resids <- function(gdat, pdat, gSigns) {
    res <- rr53(gdat, pdat, gSigns, resids=TRUE)
    tsum <- res %>% group_by(Pred) %>% do(oneSampTest(.))
    tsum$t[1]
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

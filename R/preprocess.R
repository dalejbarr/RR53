#' Read in the data from Fredrickson et al.
#'
#' This function reads in the data from the PNAS website
#' or from a downloaded CSV file and tags each row with
#' a subject ID number (SID).
#'
#' @param file Name of the local CSV file.  If left unspecified, defaults to download the data from PNAS.
#' @return A dataframe, with fields as described by Brown et al.
#' @export
readData <- function(file="http://www.pnas.org/content/suppl/2014/08/21/1407057111.DCSupplemental/pnas.1407057111.sd02.csv") {
    read.csv(file, header=TRUE) %>% mutate(., SID=row_number())
}

#' Download the data from PNAS and save it locally as a file
#'
#' @param outfile name of the file to be stored locally
#' @return No return value.
#' @export
downloadRR53Data <- function(outfile) {
    dat <- readData()
    write.csv(dat, file=outfile, row.names=FALSE)
    return()
}

#' Read in and tidy the data from Fredrickson et al.
#'
#' This downloads the data from the PNAS website
#' (or from a local CSV file) and performs some basic
#' preprocessing to create 'tidy' data.frames for
#' easy analysis.
#' @param file Name of the local CSV file.  If left unspecified, defaults to download the data from PNAS.
#' @return a list with four elements: Predictors, Genetics, Psychometrics, FactorAssignment, GeneSigns.
#' @export
readTidyData <- function(file="http://www.pnas.org/content/suppl/2014/08/21/1407057111.DCSupplemental/pnas.1407057111.sd02.csv") {
    # calculate zscores for the psychometric data
    zscorePsychometrics <- function(dat, fassign) {
        inner_join(pdata, pfactor, by="SFItem") %>% group_by(SID, Hedonic) %>%
            summarize(xbar=mean(Score, na.rm=TRUE)) %>% ungroup %>%
                group_by(Hedonic) %>% mutate(zz=scale(xbar)) %>% ungroup %>%
                    select(-xbar) %>% spread(Hedonic, zz) %>%
                        select(SID, Eud=`0`, Hed=`1`)
    }
    rr53data <- readData(file)
    # psychometric data
    pdata <- rr53data %>% select(SID, SF1:SF14) %>%
        gather(SFItem, Score, SF1:SF14, -SID) %>%
            arrange(SID, SFItem)
    # this is the assignment of scores to factors
    # by Fredrickson et al., change variable
    # "hedonic" to try other assignments
    pfactor <- pdata %>% select(SFItem) %>% unique %>%
        mutate(Hedonic=rep(c(1,0),c(3,11)))
    
    # z-score the psychometric data
    pscaled <- zscorePsychometrics(pdata, pfactor)

    # pull together all predictors including z-scored preds
    preds <- rr53data %>% select(SID, Male:ZScore_Eud) %>%
        inner_join(pscaled, by="SID")

    # pull out gene expression data
    gedata <- rr53data %>% select(SID, IFIT1L:TNF) %>%
        gather(Gene, Value, IFIT1L:TNF, -SID) %>%
            arrange(SID, Gene)

    geneSigns <- gedata %>% select(Gene) %>% unique %>%
        mutate(Sign=rep(c(-1,1), c(34,19)))

    return(list(Predictors=preds,
                Genetics=gedata,
                Psychometrics=pdata,        
                FactorAssignment=pfactor,
                GeneSigns=geneSigns))
}

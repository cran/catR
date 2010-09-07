startItems<-function (itemBank, fixItems = NULL, seed = NULL, nrItems = 1, 
    theta = 0, halfRange = 2, startSelect = "bOpt") 
{
    it <- itemBank$itemPar
    if (is.null(fixItems) == FALSE) {
        items <- fixItems
        par <- it[fixItems, ]
        thStart <- NA
    }
    else {
        if (is.numeric(seed) == TRUE) {
            set.seed(seed)
            repeat {
                items <- round(runif(nrItems, 0, 1) * nrow(it)) + 
                  1
                if (length(unique(items)) == nrItems) 
                  break
            }
            par <- it[items, ]
            thStart <- NA
        }
        else {
            if (nrItems == 1) 
                thStart <- theta
            else {
                thStart <- seq(from = theta - halfRange, to = theta + 
                  halfRange, length = nrItems)
                ra <- rnorm(nrItems)
                thStart <- thStart[rank(ra)]
            }
            if (startSelect != "bOpt" & startSelect != "MFI") 
                stop("'startSelect' must be either 'bOpt' or 'MFI'", 
                  call. = FALSE)
            if (startSelect == "bOpt") {
                items = NULL
                ind <- rep(1, nrow(it))
                for (i in 1:length(thStart)) {
                  ind[items] <- 0
                  ra<-rank(abs(thStart[i]-it[, 2]))
                  minRA<-min(ra[ind>0])
                  items <- c(items, min((1:length(ind))[ra==minRA]))
                }
            }
            if (startSelect == "MFI") {
                items = NULL
                ind <- rep(1, nrow(it))
                for (i in 1:length(thStart)) {
                  ind[items] <- 0
                  seqZE <- (1:length(itemBank$theta))
                  zeTheta <- min(seqZE[abs(itemBank$theta - thStart[i]) == 
                    min(abs(itemBank$theta - thStart[i]))])
ra<-rank(itemBank$infoTab[zeTheta, ])
                  maxRA<-max(ra[ind>0])
                  items <- c(items, min((1:length(ind))[ra==maxRA]))
                }
            }
            par = it[items, ]
        }
    }
    res <- list(items = items, par = par, thStart = thStart, 
        startSelect = startSelect)
    return(res)
}
library(sfsmisc)

MWI<-function (itemBank, item, x, it, lower = -4, upper = 4, nqp = 33, 
    type = "MLWI", priorDist="norm",priorPar = c(0, 1)) 
{
    if (type != "MLWI" & type != "MPWI") 
        stop("'type' must be either 'MLWI' or 'MPWI'", call. = FALSE)
    L <- function(th, x, par) prod(Pi(th, par)$Pi^x * (1 - Pi(th, 
        par)$Pi)^(1 - x))
    X <- seq(from = lower, to = upper, length = nqp)
    ITEMS <- NULL
    for (i in 1:length(X)) ITEMS[i] <- min((1:length(itemBank$theta))[(itemBank$theta - 
        X[i]) == min(itemBank$theta - X[i])])
    if (type == "MPWI"){
       Z<-switch(priorDist,norm=dnorm(X, priorPar[1], priorPar[2]),unif=dunif(X, priorPar[1], priorPar[2]))
        Y <- itemBank$infoTab[ITEMS, item] * sapply(X, L, x, 
            it) * Z
}
    else Y <- itemBank$infoTab[ITEMS, item] * sapply(X, L, x, 
        it)
    res <- integrate.xy(X, Y)
    return(res)
}

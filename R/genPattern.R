genPattern<-function (th, it, model = NULL, D = 1, seed = NULL) 
{
    it <- rbind(it)
if (!is.null(seed)) set.seed(seed)
    if (is.null(model)) 
        res <- rbinom(nrow(it), 1, Pi(th, it, model = model, 
            D = D)$Pi)
    else {
        pr <- Pi(th, it, model = model, D = D)$Pi
        res <- NULL
        for (i in 1:nrow(pr)) {
            pp <- pr[i, ][!is.na(pr[i, ])]
            vec <- rmultinom(n = 1, size = 1, prob = pp)
            res[i] <- (1:nrow(vec))[vec[, 1] == 1] - 1
        }
    }
    return(res)
}

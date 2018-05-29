semTheta<-function (thEst, it, x = NULL, model = NULL, D = 1, method = "BM", 
    priorDist = "norm", priorPar = c(0, 1), weight = "Huber", 
    tuCo = 1, sem.type = "classic", parInt = c(-4, 4, 33), constantPatt = NULL, 
    sem.exact = FALSE,trueTh = NULL, range = c(-4,4)) 
{
    constantPattern <- function(t) ifelse(sum(t) == 0 | sum(t) == 
        length(t), TRUE, FALSE)
    METHOD <- NULL
    if (!constantPattern(x) | !is.null(model) | is.null(constantPatt)) 
        METHOD <- method
    else {
        if (sum(constantPatt == c("BM", "EAP", "WL")) == 1) 
            METHOD <- constantPatt
        else RES <- Inf
    }
    if (!is.null(x)) {
        ind <- which(!is.na(x))
        it <- it[ind, ]
        x <- x[ind]
    }
    if (!is.null(METHOD)) {
if (sem.exact & is.null(model)){
if (!is.null(trueTh)) TH<-c(thEst, trueTh)
else TH<-thEst
 dist <- fullDist(TH, rbind(it), method = method,
                   priorDist = priorDist, priorPar = priorPar,
                   weight = weight, tuCo = tuCo, range = range, parInt=parInt)
  formula <- function(t, u) sqrt(sum(t*(u-sum(t*u))^2))
  exactSE <- formula(dist[, 2], dist[, 1])
 if (length(TH)>1) trueSE <- formula(dist[, 3], dist[, 1]) 
else trueSE<-NULL
RES<-c(exactSE,trueSE)
}
else{
        if (method == "EAP") {
            RES <- eapSem(thEst, it, x = x, model = model, D = D, 
                priorDist = priorDist, priorPar = priorPar, lower = parInt[1], 
                upper = parInt[2], nqp = parInt[3])
        }
        else {
            if (is.null(model)) {
                info <- sum(Ii(thEst, it, D = D)$Ii)
                dr0 <- function(th, it, D = 1, method = "BM", 
                  priorDist = "norm", priorPar = c(0, 1)) {
                  if (method == "BM") 
                    res <- switch(priorDist, norm = -1/priorPar[2]^2, 
                      unif = 0, Jeffreys = (sum(Ii(th, it, D = D)$d2Ii) * 
                        sum(Ii(th, it, D = D)$Ii) - sum(Ii(th, 
                        it, D = D)$dIi)^2)/(2 * sum(Ii(th, it, 
                        D = D)$Ii)^2))
                  else {
                    res <- (sum(Ji(th, it, D = D)$dJi) * sum(Ii(th, 
                      it, D = D)$Ii) - sum(Ji(th, it, D = D)$Ji) * 
                      sum(Ii(th, it, D = D)$dIi))/(2 * sum(Ii(th, 
                      it, D = D)$Ii)^2)
                  }
                  return(res)
                }
                if (METHOD == "ROB") {
                  if (sum(weight == c("Huber", "Tukey")) == 0) {
                    stop("'weight' must be either 'Huber' or 'Tukey'", 
                      call. = FALSE)
                  }
it<-rbind(it)
                  ri <- it[, 1] * (thEst - it[, 2])
                  if (weight == "Huber") {
                    wi <- tuCo/abs(ri)
                    wi[abs(ri) <= tuCo] <- 1
                  }
                  if (weight == "Tukey") {
                    wi <- (1 - (ri/tuCo)^2)^2
                    wi[abs(ri) > tuCo] <- 0
                  }
                  info <- Ii(thEst, it, D = D)$Ii
                  NUM <- sqrt(sum(wi^2 * info))
                  DEN <- sum(wi * info)
                }
                if (METHOD == "ML") {
                  NUM <- 1
                  DEN <- sqrt(info)
                }
                if (METHOD == "BM") {
                  NUM <- switch(sem.type, classic = 1, new = sqrt(info))
                  DEN <- switch(sem.type, classic = sqrt(info - 
                    dr0(thEst, it, method = "BM", priorDist = priorDist, 
                      priorPar = priorPar)), new = abs(info - 
                    dr0(thEst, it, method = "BM", priorDist = priorDist, 
                      priorPar = priorPar)))
                }
                if (METHOD == "WL") {
                  NUM <- switch(sem.type, classic = 1, new = sqrt(info))
                  DEN <- switch(sem.type, classic = sqrt(info), 
                    new = abs(info - dr0(thEst, it, method = "WL", 
                      priorDist = priorDist, priorPar = priorPar)))
                }
                RES <- NUM/DEN
            }
            else {
                met <- switch(method, ML = 1, BM = 2, WL = 3, 
                  EAP = 4)
                pd <- switch(priorDist, norm = 1, unif = 2, Jeffreys = 3)
                if (met == 1 | (met == 2 & pd == 2)) 
                  optI <- sum(Ii(thEst, it, model = model, D = D)$Ii)
                if (met == 2 & pd == 1) 
                  optI <- sum(Ii(thEst, it, model = model, D = D)$Ii) + 
                    1/priorPar[2]^2
                if ((met == 2 & pd == 3) | met == 3) {
                  prI <- Ii(thEst, it, model = model, D = D)
                  prJ <- Ji(thEst, it, model = model, D = D)
                  if (met == 2) 
                    optI <- sum(prI$Ii) + (sum(prI$dIi)^2 - sum(prI$d2Ii) * 
                      sum(prI$Ii))/(2 * sum(prI$Ii)^2)
                  else optI <- sum(prI$Ii)
                }
                RES <- 1/sqrt(optI)
            }
        }
    }

}
    return(RES)
}

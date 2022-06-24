\name{semTheta}
\alias{semTheta}

\title{Standard error of ability estimation (dichotomous and polytomous models)}

\description{
 This command returns the estimated standard error of the ability estimate, for a given response pattern and a given matrix of item parameters, either under the 4PL model or any suitable polytomous IRT model. Exact standard errors are available for dichotomous models.
 }

\usage{
semTheta(thEst, it, x = NULL, model = NULL, D = 1, method = "BM", 
  	priorDist = "norm", priorPar = c(0, 1), weight = "Huber", tuCo = 1, 
  	sem.type = "classic", parInt = c(-4, 4, 33), constantPatt = NULL,
  	sem.exact = FALSE, trueTh = NULL, range = c(-4, 4))
 }

\arguments{
 \item{thEst}{numeric: the ability estimate.}
 \item{it}{numeric: a suitable matrix of item parameters. See \bold{Details}.}
 \item{x}{numeric: a vector of item responses (default is \code{NULL}). Can also hold missing data (coded as \code{NA}s). Ignored if \code{method} is not \code{"EAP"} but must be provided in presence of missing responses. See \bold{Details}.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952).} 
 \item{method}{character: the ability estimator. Possible values are \code{"BM"} (default), \code{"ML"}, \code{"WL"}, \code{"EAP"} and \code{"ROB"}. See \bold{Details}.}
 \item{priorDist}{character: specifies the prior distribution. Possible values are \code{"norm"} (default), \code{"unif"} and \code{"Jeffreys"}.
 Ignored if \code{method} is neither \code{"BM"} nor \code{"EAP"}. See \bold{Details}.}
 \item{priorPar}{numeric: vector of two components specifying the prior parameters (default is \code{c(0,1)}) of the prior ability distribution. Ignored if \code{method} is neither \code{"BM"}
 nor \code{"EAP"}, or if \code{priorDist="Jeffreys"}. See \bold{Details}.}
\item{weight}{character: the type of weight function for the robust estimator. Possible values are \code{"Huber"} (default) and \code{"Tukey"}. Ignored if \code{method} is not \code{"ROB"} or if \code{model} is not \code{NULL}. See \bold{Details}.}
\item{tuCo}{numeric: the value of the tuning constant for the weight function (default is 1, suitable with \code{"Huber"} weight). Ignored if \code{method} is not \code{"ROB"} or if \code{model} is not \code{NULL}. See \bold{Details}.}
\item{sem.type}{character: the type of ASE formula to be used, either \code{"classic"} (default) or \code{"new"}. Ignored if \code{method} is neither \code{"BM"} nor \code{"WL"}, or if \code{model} is not \code{NULL}. See \bold{Details}.}
 \item{parInt}{numeric: vector of three components, holding respectively the values of the arguments \code{lower}, \code{upper} and \code{nqp} of the
 \code{\link{eapEst}} command. Default vector is (-4, 4, 33). Ignored if \code{method} is not \code{"EAP"}.}
\item{constantPatt}{character: the method to estimate ability in case of constant pattern (i.e. only correct or only incorrect responses). Can be eitehr \code{NULL} (default), \code{"BM"}, \code{"EAP"}, \code{"WL"}, \code{"fixed4"}, \code{"fixed7"} or \code{"var"}. \emph{Currently only implemented for dichotomous IRT models}. See \bold{Details}.}
\item{sem.exact}{logical: should \emph{exact} standard error be computed instead of asymptotic standard error? (default is \code{FALSE}). Ignored if \code{model} is not \code{NULL}. See \bold{Details}.}
\item{trueTh}{(For simulation study purposes only) either \code{NULL} (default) or the true ability level of interest. Ignored if \code{sem.exact} is \code{FALSE} or if \code{model} is not \code{NULL}. See \bold{Details}.}
 \item{range}{numeric: vector of two components specifying the range wherein the ability estimate must be looked for (default is \code{c(-4,4)}). Used only to compute exact standard errors. Ignored if \code{sem.exact} is \code{FALSE} or if \code{model} is not \code{NULL}.}
}

\value{
 The estimated standard error of the ability level (or \code{Inf} if the response pattern is constant and \code{constantPatt} is not \code{NULL}). 

If exact standard error is computed (\code{sem.exact = TRUE}) and \code{trueTh} argument takes some numeric value, then output is a vector of two components with the two exact standard errors, the first SE corresponding to the \code{thEst} value and the second SE to the \code{trueTh} value.
 }
 
\details{
Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{it} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model
 (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model, \code{"PCM"} for Partical Credit Model, \code{"GPCM"} for Generalized Partial Credit Model, \code{"RSM"} for Rating Scale Model and \code{"NRM"} for Nominal Response Model. The \code{it} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.

The vector of response patterns \code{x} can also hold missing responses (more useful in linear testing, not in CAT). In this case the missing responses must be coded as \code{NA} values. They are discarded from the ability estimation process. Note that in presence of missing responses, the pattern \code{x} must be provided indepedently of the estimation method (to discard items with missing responses from the computation). 

 Five ability estimators are available: the maximum likelihood (ML) estimator (Lord, 1980), the Bayes modal (BM) estimator (Birnbaum, 1969), the expected a posteriori (EAP) estimator (Bock and Mislevy, 1982), the weighted likelihood (WL) estimator (Warm, 1989) and the robust estimator (Schuster & Yuan, 2011). The selected estimator is specified by the \code{method} argument, with values \code{"ML"}, \code{"BM"}, \code{"EAP"}, \code{"WL"} and \code{"ROB"} respectively.

 For the BM and EAP estimators, three prior distributions are available: the normal distribution, the uniform distribution and the Jeffreys' prior distribution (Jeffreys, 1939, 1946). The prior distribution is specified by the argument \code{priorPar}, with values \code{"norm"}, \code{"unif"} and \code{"Jeffreys"}, respectively. The \code{priorPar} argument is ignored if \code{method="ML"} or \code{method="WL"}.
 
 The argument \code{priorPar} determines either: the prior mean and standard deviation of the normal prior distribution (if \code{priorDist="norm"}), or the range for defining the prior uniform distribution (if \code{priorDist="unif"}). This argument is ignored if \code{priorDist="Jeffreys"}.
 
 The \code{eapPar} argument sets the range and the number of quadrature points for numerical integration in the EAP process. By default, it takes the vector value (-4, 4, 33), that is, 33 quadrature points on the range [-4; 4] (or, by steps of 0.25). See \code{\link{eapEst}} for further details. 

Robust estimation requires an appropriate weight function that depends on an accurate tuning constant. Suggested functions are the Huber weight (Schuester and Yuan, 2011) and the Tukey weight (Mosteller and Tukey, 1977). Both can be set by the \code{weight} argument, with respective values \code{"Huber"} and \code{"Tukey"}. Default function is Huber. Moreover, the \code{tuCo} argument specifies the tuning constant for the weight function. Default value is 1 and suggested for Huber weight (also by default), and value 4 is suggested for Tukey weight (Schuester and Yuan, 2011).

New ASE formulas, proposed by Magis (2016), are now available. They can be supplied by the \code{sem.type} argument, which takes the default value \code{"classic"} (so usual ASEs are computed) or \code{"new"} (for the newly suggested formulas). Note that new ASEs are available only for BM and WL estimators, as for other estimators the classic and new versions are identical. Note also that these new ASE formulas are available for dichotomous IRT models only.

 Note that in the current version, the ability estimate must be specified through the \code{thEst} argument. Moreover, the response pattern must be specified through the \code{x} argument to compute the standard error of the EAP estimate. For the other estimation methods, this is not  necessary, and \code{x} is set to \code{NULL} by default for this purpose.

Note also that if specific stepsize adjustment was required for constant patterns with the \code{constantPatt} argument (that is, if it takes value \code{"fixed4"}, \code{"fixed7"} or \code{"var"}) then an infinite value \code{Inf} is being returned. 

Finally, exact standard errors can be computed with dichotomous IRT models only (Magis, 2014). This is requested by setting argument \code{sem.exact} to \code{TRUE} (default is \code{FALSE} so regular, asymptotic formula is considered). The exact standard error is computed using the full ability distribution provided by the \code{link{fullDist}} function. Two additional arguments can be set up: \code{range} to define the range of estimated ability levels during the derivation of the full distribution, and \code{trueTh} that can be used to specify the true ability level in addition to the estimated one (through \code{thEst} argument). The latter argument is for simulation study purposes, as the set of all observable ability estimates in \code{\link{fullDist}} must then be computed only once. If \code{trueTh} is provided some real value, then two SEs are returned: the \emph{exact} SE obtained with the estimated \code{thEst} level, and the \emph{true} SE derived as the exact SE with true \code{trueTh} level.

\bold{Important}: The computation of the exact standard error is very efficient under the Rasch (1PL) model due to the ue of the Lord-Wingersky algorithm (Lord and Wingersky, 1984).  Therefore, any test length can be considered for computing the exact SE under this model. But for other dichotomous IRT models, the computatinal effort becomes very demanding even with a limited set of \eqn{n} items, since \eqn{2^n} patterns must be generated and ability estimated with each such pattern. For this reason, it is recommended not to consider exact SE with models other than the Rasch (1PL) with more than 10 items.
 }


\note{
The classic asymptotic standard error of the WL estimator is computed with the same formula as that of the ML estimator (up to the plug-in of the WL estimate instead of the ML estimate). Note however that versions of \code{catR} prior to 3.0 hold a different formula mentioned in Magis and raiche (2012), but it appeared that this formula can lead to negative values of the square of the standard error. So the usual suggestion by Warm (1989) of using the same asymptotic formulas for ML and WL is the currently used formula for classic asymptotic computation of the standard error.
}

\references{
 Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}.
 Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Birnbaum, A. (1969). Statistical theory for logistic mental test models with a prior distribution of ability. \emph{Journal of Mathematical Psychology, 6}, 258-276. \doi{10.1016/0022-2496(69)90005-4}

 Bock, R. D., and Mislevy, R. J. (1982). Adaptive EAP estimation of ability in a microcomputer environment. \emph{Applied Psychological Measurement, 6}, 431-444. \doi{10.1177/014662168200600405}

Dodd, B. G., De Ayala, R. J., and Koch, W. R. (1995). Computerized adaptive testing with polytomous items. \emph{Applied Psychological Measurement, 19}, 5-22. \doi{10.1177/014662169501900103}

 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error}. 
 Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 
 
 Jeffreys, H. (1939). \emph{Theory of probability}. Oxford, UK: Oxford University Press.

 Jeffreys, H. (1946). An invariant form for the prior probability in estimation problems. \emph{Proceedings of the Royal Society of London. Series A, Mathematical and Physical Sciences, 186}, 453-461.

 Lord, F.M. (1980). \emph{Applications of item response theory to practical testing problems}. Hillsdale, NJ: Lawrence Erlbaum.

Lord, F. M., and Wingersky, M. S. (1984). Comparison of IRT true-score and equipercentile observed-score equatings. \emph{Applied Psychological Measurement, 8}, 453-461.  \doi{10.1177/014662168400800409}

Magis, D. (2014). Accuracy of asymptotic standard errors of the maximum and weighted likelihood estimators of proficiency levels with short tests. \emph{Applied Psychological Measurement, 38}, 105-121. \doi{10.1177/0146621613496890}

Magis, D. (2016). Efficient standard errors formulas of ability estimators with dichotomous item response models. \emph{Psychometrika, 81}, 184-200. \doi{10.1007/s11336-015-9443-3}

Magis, D. and Barrada, J. R. (2017). Computerized Adaptive Testing with R: Recent Updates of the Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{Code Snippets}, \emph{76(1)}, 1-18. \doi{10.18637/jss.v076.c01}

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. \doi{10.18637/jss.v048.i08}

Mosteller, F., and Tukey, J. (1977). \emph{Exploratory data analysis and regression}. Reading, MA: Addison-Wesley.

Schuester, C., and Yuan, K.-H. (2011). Robust estimation of latent ability in item response models. \emph{Journal of Educational and Behavioral Statistics, 36}, 720)735.\doi{10.3102/1076998610396890}

 Warm, T.A. (1989). Weighted likelihood estimation of ability in item response models. \emph{Psychometrika, 54}, 427-450. \doi{10.1007/BF02294627}
}

\author{
    David Magis \cr
    Department of Psychology, University of Liege, Belgium \cr
    \email{david.magis@uliege.be}

    Lianne Ippel \cr
    Department of Psychology, University of Liege, Belgium \cr
    \email{g.j.e.ippel@gmail.com}
 }

\seealso{
\code{\link{eapSem}}, \code{\link{thetaEst}}, \code{\link{genPolyMatrix}}, \code{\link{fullDist}}

} 

\examples{

## Dichotomous models ##

 # Loading the 'tcals' parameters 
 data(tcals)
 
 # Selecting item parameters only
 tcals <- as.matrix(tcals[,1:4])

 # Creation of a response pattern (tcals item parameters, true ability level 0)
 x <- genPattern(0, tcals, seed = 1)

 # ML estimation
 th <- thetaEst(tcals, x, method = "ML")
 c(th, semTheta(th, tcals, method = "ML"))

 # ML estimation, new ASE formula (=> yields the same result)
 c(th, semTheta(th, tcals, method = "ML", sem.type = "new"))

 # With first two responses missing
 x.mis <- x
 x.mis[1:2] <- NA
 th <- thetaEst(tcals, x.mis, method = "ML")
 c(th, semTheta(th, tcals, x.mis, method = "ML"))

 # BM estimation, standard normal prior distribution
 th <- thetaEst(tcals, x)
 c(th, semTheta(th, tcals))

 # BM estimation and new ASE formula
 c(th, semTheta(th, tcals, sem.type = "new"))

 # BM estimation, uniform prior distribution upon range [-2,2]
 th <- thetaEst(tcals, x, method = "BM", priorDist = "unif",
                priorPar = c(-2, 2))
 c(th, semTheta(th, tcals, method = "BM", priorDist = "unif",
		    priorPar = c(-2, 2)))

 # BM estimation, Jeffreys' prior distribution  
 th <- thetaEst(tcals, x, method = "BM", priorDist = "Jeffreys")
 c(th, semTheta(th, tcals, method = "BM", priorDist = "Jeffreys"))

 # EAP estimation, standard normal prior distribution
 th <- thetaEst(tcals, x, method = "EAP")
 c(th, semTheta(th, tcals, x, method = "EAP"))

\dontrun{

 # EAP estimation, uniform prior distribution upon range [-2,2]
 th <- thetaEst(tcals, x, method = "EAP", priorDist = "unif",
                priorPar = c(-2, 2))
 c(th, semTheta(th, tcals, x, method = "EAP", priorDist = "unif",
		    priorPar = c(-2, 2)))

 # EAP estimation, Jeffreys' prior distribution  
 th <- thetaEst(tcals, x, method = "EAP", priorDist = "Jeffreys")
 c(th, semTheta(th, tcals, x, method = "EAP", priorDist = "Jeffreys"))

 # WL estimation
 th <- thetaEst(tcals, x, method = "WL")
 c(th, semTheta(th, tcals, method = "WL"))

 # WL estimation, new ASE formula
 c(th, semTheta(th, tcals, method = "WL", sem.type = "new"))

 # 'fixed4' adjustment for constant pattern
 th <- thetaEst(tcals, rep(0, nrow(tcals)), constantPatt = "fixed4")
 c(th, semTheta(th, tcals, constantPatt = "fixed4"))

 # Robust estimation
 th <- thetaEst(tcals, x, method = "ROB")
 c(th, semTheta(th, tcals, method = "ROB"))

 # Robust estimation, Huber weight and tuning constant 2
 th <- thetaEst(tcals, x, method = "ROB", tuCo = 2)
 c(th, semTheta(th, tcals, method = "ROB", tuCo = 2))

 # Robust estimation, Tukey weight and tuning constant 4
 th <- thetaEst(tcals, x, method = "ROB", weight = "Tukey", tuCo = 4)
 c(th, semTheta(th, tcals, method = "ROB", weight = "Tukey", tuCo = 4))

 ## Exact SE computation under 1PL model: 
 # Creation of a 1PL item bank with difficulties from 'tcals' (85 items)
 tcals2 <- cbind(1, tcals[, 2],  0, 1)

 # Pattern generation for true ability level 1
 x2 <- genPattern(1, tcals2, seed = 1)

 # ML estimation
 th2 <- thetaEst(tcals2, x2, method = "ML")
 c(th2, semTheta(th2, tcals2, x2, method = "ML", sem.exact = TRUE))

 # ML estimation, true SE in addition
 c(th2, semTheta(th2, tcals2, x2, method = "ML", sem.exact = TRUE, trueTh = 1))

 ## Exact SE computation under 2PL model: 
 # Creation of a 2PL item bank with ten items
 it <- genDichoMatrix(10, model = "2PL", seed = 1)

 # Pattern generation for true ability level 1
 x3 <- genPattern(1, it, seed = 1)

 # ML estimation
 th3 <- thetaEst(it, x3, method = "ML")
 c(th3, semTheta(th3, it, x3, method = "ML", sem.exact = TRUE))

 # ML estimation, true SE in addition
 c(th3, semTheta(th3, it, x3, method = "ML", sem.exact = TRUE, trueTh = 1))
}

## Polytomous models ##

 # Generation of an item bank under GRM with 100 items and at most 4 categories
 m.GRM <- genPolyMatrix(100, 4, "GRM")
 m.GRM <- as.matrix(m.GRM)

 # Creation of a response pattern (true ability level 0)
 set.seed(1)
 x <- genPattern(0, m.GRM, model = "GRM")

# ML estimation
 th <- thetaEst(m.GRM, x, model = "GRM", method = "ML")
 c(th, semTheta(th, m.GRM, model = "GRM", method = "ML"))

 # BM estimation, standard normal prior distribution
 th <- thetaEst(m.GRM, x, model = "GRM")
 c(th, semTheta(th, m.GRM, model = "GRM"))

 # BM estimation, uniform prior distribution upon range [-2,2]
 th <- thetaEst(m.GRM, x, model = "GRM", method = "BM", priorDist = "unif", 
    priorPar = c(-2, 2))
 c(th, semTheta(th, m.GRM, model = "GRM", method = "BM", priorDist = "unif", 
  priorPar = c(-2, 2)))

 # BM estimation, Jeffreys' prior distribution  
 th <- thetaEst(m.GRM, x, model = "GRM", method = "BM", priorDist = "Jeffreys")
 c(th, semTheta(th, m.GRM, model = "GRM", method = "BM", priorDist = "Jeffreys"))

 # EAP estimation, standard normal prior distribution
 th <- thetaEst(m.GRM, x, model = "GRM", method = "EAP")
 c(th, semTheta(th, m.GRM, x, model = "GRM", method = "EAP") )

\dontrun{

 # EAP estimation, uniform prior distribution upon range [-2,2]
 th <- thetaEst(m.GRM, x, model = "GRM", method = "EAP", priorDist = "unif", 
    priorPar = c(-2, 2))
 c(th, semTheta(th, m.GRM, x, model = "GRM", method = "EAP", priorDist = "unif", 
  priorPar = c(-2, 2)))

 # EAP estimation, Jeffreys' prior distribution  
 th <- thetaEst(m.GRM, x, model = "GRM", method = "EAP", priorDist = "Jeffreys")
 c(th, semTheta(th, m.GRM, x, model = "GRM", method = "EAP", priorDist = "Jeffreys"))

 # WL estimation
 th <- thetaEst(m.GRM, x, model = "GRM", method = "WL")
 c(th, semTheta(th, m.GRM, model = "GRM", method = "WL"))


 # Loading the cat_pav data
 data(cat_pav)
 cat_pav <- as.matrix(cat_pav)

 # Creation of a response pattern (true ability level 0)
 set.seed(1)
 x <- genPattern(0, cat_pav, model = "GPCM")

# ML estimation
 th <- thetaEst(cat_pav, x, model = "GPCM", method = "ML")
 c(th, semTheta(th, cat_pav, model = "GPCM", method = "ML"))

 # BM estimation, standard normal prior distribution
 th <- thetaEst(cat_pav, x, model = "GPCM")
 c(th, semTheta(th, cat_pav, model = "GPCM"))

 # BM estimation, uniform prior distribution upon range [-2,2]
 th <- thetaEst(cat_pav, x, model = "GPCM", method = "BM", priorDist = "unif", 
    priorPar = c(-2, 2))
 c(th, semTheta(th, cat_pav, model = "GPCM", method = "BM", priorDist = "unif", 
  priorPar = c(-2, 2)))

 # BM estimation, Jeffreys' prior distribution  
 th <- thetaEst(cat_pav, x, model = "GPCM", method = "BM", priorDist = "Jeffreys")
 c(th, semTheta(th, cat_pav, model = "GPCM", method = "BM", priorDist = "Jeffreys"))

 # EAP estimation, standard normal prior distribution
 th <- thetaEst(cat_pav, x, model = "GPCM", method = "EAP")
 c(th, semTheta(th, cat_pav, x, model = "GPCM", method = "EAP"))

 # EAP estimation, uniform prior distribution upon range [-2,2]
 th <- thetaEst(cat_pav, x, model = "GPCM", method = "EAP", priorDist = "unif", 
    priorPar = c(-2, 2))
 c(th, semTheta(th, cat_pav, x, model = "GPCM", method = "EAP", priorDist = "unif", 
  priorPar = c(-2, 2)))

 # EAP estimation, Jeffreys' prior distribution  
 th <- thetaEst(cat_pav, x, model = "GPCM", method = "EAP", priorDist = "Jeffreys")
 c(th, semTheta(th, cat_pav, x, model = "GPCM", method = "EAP", priorDist = "Jeffreys"))

 # WL estimation
 th <- thetaEst(cat_pav, x, model = "GPCM", method = "WL")
 c(th, semTheta(th, cat_pav, model = "GPCM", method = "WL"))
}
 }

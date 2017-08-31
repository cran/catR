\name{semTheta}
\alias{semTheta}

\title{Standard error of ability estimation (dichotomous and polytomous models)}

\description{
 This command returns the estimated standard error of the ability estimate, for a given response pattern and a given matrix of item parameters, either under the 4PL model or any suitable polytomous IRT model.
 }

\usage{
semTheta(thEst, it, x = NULL, model = NULL, D = 1, method = "BM", 
  	priorDist = "norm", priorPar = c(0, 1), parInt = c(-4, 4, 33),
 	 constantPatt = NULL)
 }

\arguments{
 \item{thEst}{numeric: the ability estimate.}
 \item{it}{numeric: a suitable matrix of item parameters. See \bold{Details}.}
 \item{x}{numeric: a vector of item responses (default is \code{NULL}). Can also hold missing data (coded as \code{NA}s). Ignored if \code{method} is not \code{"EAP"} but must be provided in presence of missing responses. See \bold{Details}.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952).} 
 \item{method}{character: the ability estimator. Possible values are \code{"BM"} (default), \code{"ML"}, \code{"WL"} and \code{"EAP"}. See \bold{Details}.}
 \item{priorDist}{character: specifies the prior distribution. Possible values are \code{"norm"} (default), \code{"unif"} and \code{"Jeffreys"}.
 Ignored if \code{method} is neither \code{"BM"} nor \code{"EAP"}. See \bold{Details}.}
 \item{priorPar}{numeric: vector of two components specifying the prior parameters (default is \code{c(0,1)}) of the prior ability distribution. Ignored if \code{method} is neither \code{"BM"}
 nor \code{"EAP"}, or if \code{priorDist="Jeffreys"}. See \bold{Details}.}
 \item{parInt}{numeric: vector of three components, holding respectively the values of the arguments \code{lower}, \code{upper} and \code{nqp} of the
 \code{\link{eapEst}} command. Default vector is (-4, 4, 33). Ignored if \code{method} is not \code{"EAP"}.}
\item{constantPatt}{character: the method to estimate ability in case of constant pattern (i.e. only correct or only incorrect responses). Can be eitehr \code{NULL} (default), \code{"BM"}, \code{"EAP"}, \code{"WL"}, \code{"fixed4"}, \code{"fixed7"} or \code{"var"}. \emph{Currently only implemented for dichotomous IRT models}. See \bold{Details}.}
}

\value{
 The estimated standard error of the ability level (or \code{Inf} if the response pattern is constant and \code{constantPatt} is not \code{NULL}).
 }
 
\details{
Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{it} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model
 (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model, \code{"PCM"} for Partical Credit Model, \code{"GPCM"} for Generalized Partial Credit Model, \code{"RSM"} for Rating Scale Model and \code{"NRM"} for Nominal Response Model. The \code{it} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.

The vector of response patterns \code{x} can also hold missing responses (more useful in linear testing, not in CAT). In this case the missing responses must be coded as \code{NA} values. They are discarded from the ability estimation process. Note that in presence of missing responses, the pattern \code{x} must be provided indepedently of the estimation method (to discard items with missing responses from the computation). 

Four ability estimators are available: the maximum likelihood (ML) estimator (Lord, 1980), the Bayes modal (BM) estimator (Birnbaum, 1969), the expected a posteriori (EAP) estimator (Bock and Mislevy, 1982)and the weighted likelihood (WL) estimator (Warm, 1989). The selected estimator is specified by the \code{method} argument, with values \code{"ML"}, \code{"BM"}, \code{"EAP"} and \code{"WL"} respectively.

 For the BM and EAP estimators, three prior distributions are available: the normal distribution, the uniform distribution and the Jeffreys' prior distribution (Jeffreys, 1939, 1946). The prior distribution is specified by the argument \code{priorPar}, with values \code{"norm"}, \code{"unif"} and \code{"Jeffreys"}, respectively. The \code{priorPar} argument is ignored if \code{method="ML"} or \code{method="WL"}.
 
 The argument \code{priorPar} determines either: the prior mean and standard deviation of the normal prior distribution (if \code{priorDist="norm"}), or the range for defining the prior uniform distribution (if \code{priorDist="unif"}). This argument is ignored if \code{priorDist="Jeffreys"}.
 
 The \code{eapPar} argument sets the range and the number of quadrature points for numerical integration in the EAP process. By default, it takes the vector value (-4, 4, 33), that is, 33 quadrature points on the range [-4; 4] (or, by steps of 0.25). See \code{\link{eapEst}} for further details. 

 Note that in the current version, the ability estimate must be specified through the \code{thEst} argument. Moreover, the response pattern must be specified through the \code{x} argument to compute the standard error of the EAP estimate. For the other estimation methods, this is not  necessary, and \code{x} is set to \code{NULL} by default for this purpose.

Note also that if specific stepsize adjustment was required for constant patterns with the \code{constantPatt} argument (that is, if it takes value \code{"fixed4"}, \code{"fixed7"} or \code{"var"}) then an infinite value \code{Inf} is being returned. 
 }


\note{
Currently the standard error of the WL estimator is computed with the same formula as that of the ML estimator (up to the plug-in of the WL estimate instead of the ML estimate). Version of \code{catR} prior to 3.0 holds a different formula mentioned in Magis and raiche (2012), but it appeared that this formula can lead to negative values of the square of the standard error. So the usual suggestion by Warm (1989) of using the same asymptotic formulas for ML and WL is currently in application until a corrected formula can be provided.
}

\references{
 Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}.
 Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Birnbaum, A. (1969). Statistical theory for logistic mental test models with a prior distribution of ability. \emph{Journal of Mathematical Psychology, 6}, 258-276. \doi{10.1016/0022-2496(69)90005-4}

 Bock, R. D., and Mislevy, R. J. (1982). Adaptive EAP estimation of ability in a microcomputer environment. \emph{Applied Psychological Measurement, 6}, 431-444. \doi{10.1177/014662168200600405}

Dodd, B. G., De Ayala, R. J., and Koch, W. R. (1995) . Computerized adaptive testing with polytomous items. \emph{Applied Psychological Measurement, 19}, 5-22. \doi{10.1177/014662169501900103}

 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error}. 
 Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 
 
 Jeffreys, H. (1939). \emph{Theory of probability}. Oxford, UK: Oxford University Press.

 Jeffreys, H. (1946). An invariant form for the prior probability in estimation problems. \emph{Proceedings of the Royal Society of London. Series A, Mathematical and Physical Sciences, 186}, 453-461.

 Lord, F.M. (1980). \emph{Applications of item response theory to practical testing problems}. Hillsdale, NJ: Lawrence Erlbaum.

Magis, D. and Barrada, J. R. (2017). Computerized Adaptive Testing with R: Recent Updates of the Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{Code Snippets}, \emph{76(1)}, 1-18. \doi{10.18637/jss.v076.c01}

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. \doi{10.18637/jss.v048.i08}

 Warm, T.A. (1989). Weighted likelihood estimation of ability in item response models. \emph{Psychometrika, 54}, 427-450. \doi{10.1007/BF02294627}
}

\author{
    David Magis \cr
    Department of Education, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{eapSem}}, \code{\link{thetaEst}}, \code{\link{genPolyMatrix}}
} 

\examples{

## Dichotomous models ##

 # Loading the 'tcals' parameters 
 data(tcals)
 
 # Selecting item parameters only
 tcals <- as.matrix(tcals[,1:4])

 # Creation of a response pattern (tcals item parameters, true ability level 0)
 set.seed(1)
 x <- genPattern(0, tcals)

 # ML estimation
 th <- thetaEst(tcals, x, method = "ML")
 c(th, semTheta(th, tcals, method = "ML"))

 # With first two responses missing
 x.mis <- x
 x.mis[1:2] <- NA
 th <- thetaEst(tcals, x.mis, method = "ML")
 c(th, semTheta(th, tcals, x.mis, method = "ML"))

 # BM estimation, standard normal prior distribution
 th <- thetaEst(tcals, x)
 c(th, semTheta(th, tcals))

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

 # 'fixed4' adjustment for constant pattern
 th <- thetaEst(tcals, rep(0, nrow(tcals)), constantPatt = "fixed4")
 c(th, semTheta(th, tcals, constantPatt = "fixed4"))
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

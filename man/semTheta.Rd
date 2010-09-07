\name{semTheta}
\alias{semTheta}

\title{Standard error of ability estimation under the 4PL model}

\description{
 This command returns the estimated standard error of the ability estimate, for a given matrix of
 item parameters of the 4PL model, an ability estimate and a specified estimator. 
 }

\usage{
 semTheta(thEst, it, x=NULL, D=1, method="BM", priorDist="norm",
 	priorPar=c(0,1), parInt=c(-4,4,33))
 }

\arguments{
 \item{thEst}{numeric: the ability estimate.}
 \item{it}{numeric: a matrix with one row per item and four columns, with the values of the discrimination, the difficulty,
           the pseudo-guessing and the inattention parameters (in this order).}
 \item{x}{numeric: a vector of dichotomous item responses (default is \code{NULL}). Ignored if \code{method} is not \code{"EAP"}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952).} 
 \item{method}{character: the ability estimator. Possible values are \code{"BM"} (default), \code{"ML"} and \code{"WL"}. See \bold{Details}.}
 \item{priorDist}{character: specifies the prior distribution. Possible values are \code{"norm"} (default), \code{"unif"} and \code{"Jeffreys"}.
 Ignored if \code{method} is neither \code{"BM"} nor \code{"EAP"}. See \bold{Details}.}
 \item{priorPar}{numeric: vector of two components specifying the prior parameters (default is \code{c(0,1)}) of the prior ability distribution. Ignored if \code{method} is neither \code{"BM"}
 nor \code{"EAP"}, or if \code{priorDist="Jeffreys"}. See \bold{Details}.}
 \item{parInt}{numeric: vector of three components, holding respectively the values of the arguments \code{lower}, \code{upper} and \code{nqp} of the
 \code{\link{eapEst}} command. Default vector is (-4, 4, 33). Ignored if \code{method} is not \code{"EAP"}.}
}

\value{
 The estimated standard error of the ability level.
 }
 
\details{
 Four ability estimators are available: the maximum likelihood (ML) estimator (Lord, 1980), the Bayes modal (BM) estimator (Birnbaum, 1969),
 the expected a posteriori (EAP) estimator (Bock and Mislevy, 1982)and the weighted likelihood (WL) estimator (Warm, 1989). The selected 
 estimator is specified by the \code{method} argument, with values \code{"ML"}, \code{"BM"}, \code{"EAP"} and \code{"WL"} respectively.

 For the BM and EAP estimators, three prior distributions are available: the normal distribution, the uniform distribution and the Jeffreys' prior
 distribution (Jeffreys, 1939, 1946). The prior distribution is specified by the argument \code{priorPar}, with values
 \code{"norm"}, \code{"unif"} and \code{"Jeffreys"}, respectively. The \code{priorPar} argument is ignored if \code{method="ML"}
 or \code{method="WL"}.
 
 The argument \code{priorPar} determines either: the prior mean and standard deviation of the normal prior distribution (if 
 \code{priorDist="norm"}), or the range for defining the prior uniform distribution (if \code{priorDist="unif"}). This argument
 is ignored if \code{priorDist="Jeffreys"}.
 
 The \code{eapPar} argument sets the range and the number of quadrature points for numerical integration in the EAP process. By default, it takes
 the vector value (-4, 4, 33), that is, 33 quadrature points on the range [-4; 4] (or, by steps of 0.25). See \code{\link{eapEst}} for further details. 

 Note that in the current version, the ability estimate must be specified through the \code{thEst} argument. Moreover, the response pattern must
 be specified through the \code{x} argument to compute the standard error of the EAP estimate. For the other estimation methods, thi is not 
 necessary, and \code{x} is set to \code{NULL} by default for this purpose.
 }

\references{
 Barton, M.A., and Lord, F.M. (1981). An upper asymptote for the three-parameter logistic item-response model.
 Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Birnbaum, A. (1969). Statistical theory for logistic mental test models with a prior distribution of ability.
 \emph{Journal of Mathematical Psychology, 6}, 258-276.

 Bock, R. D., and Mislevy, R. J. (1982). Adaptive EAP estimation of ability in a microcomputer environment. 
 \emph{Applied Psychological Measurement, 6}, 431-444.

 Haley, D.C. (1952). Estimation of the dosage mortality relationship when the dose is subject to error. 
 Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 
 
 Jeffreys, H. (1939). \emph{Theory of probability}. Oxford, UK: Oxford University Press.

 Jeffreys, H. (1946). An invariant form for the prior probability in estimation problems. \emph{Proceedings of the Royal
 Society of London. Series A, Mathematical and Physical Sciences, 186}, 453-461.

 Lord, F.M. (1980). \emph{Applications of item response theory to practical testing problems}.
 Hillsdale, NJ: Lawrence Erlbaum.

 Warm, T.A. (1989). Weighted likelihood estimation of ability in item response models. \emph{Psychometrika, 54}, 427-450.
}

\author{
    David Magis \cr
    Post-doctoral researcher, FNRS (Fonds National de la Recherche Scientifique) \cr
    Department of Mathematics, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{eapSem}}, \code{\link{thetaEst}}
} 

\examples{
\dontrun{

 # Loading the 'tcals' parameters 
 data(tcals)
 tcals <- as.matrix(tcals)

 # Creation of a response pattern (tcals item parameters,
 # true ability level 0)
 set.seed(1)
 x <- rbinom(85, 1, Pi(0, tcals)$Pi)

 # ML estimation
 th <- thetaEst(tcals, x, method="ML")
 c(th, semTheta(th, tcals, method="ML"))

 # BM estimation, standard normal prior distribution
 th <- thetaEst(tcals, x)
 c(th, semTheta(th, tcals))

 # BM estimation, uniform prior distribution upon range [-2,2]
 th <- thetaEst(tcals, x, method="BM", priorDist="unif",
                priorPar=c(-2,2))
 c(th, semTheta(th, tcals, method="BM", priorDist="unif",
		    priorPar=c(-2,2)))

 # BM estimation, Jeffreys' prior distribution  
 th <- thetaEst(tcals, x, method="BM", priorDist="Jeffreys")
 c(th, semTheta(th, tcals, method="BM", priorDist="Jeffreys"))

 # EAP estimation, standard normal prior distribution
 th <- thetaEst(tcals, x, method="EAP")
 c(th, semTheta(th, tcals, x, method="EAP"))

 # EAP estimation, uniform prior distribution upon range [-2,2]
 th <- thetaEst(tcals, x, method="EAP", priorDist="unif",
                priorPar=c(-2,2))
 c(th, semTheta(th, tcals, x, method="EAP", priorDist="unif",
		    priorPar=c(-2,2)))

 # EAP estimation, Jeffreys' prior distribution  
 th <- thetaEst(tcals, x, method="EAP", priorDist="Jeffreys")
 c(th, semTheta(th, tcals, x, method="EAP", priorDist="Jeffreys"))

 # WL estimation
 th <- thetaEst(tcals, x, method="WL")
 c(th, semTheta(th, tcals, method="WL"))
 }
 }

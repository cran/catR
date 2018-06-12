\name{eapSem}
\alias{eapSem}

\title{Standard error of EAP ability estimation (dichotomous and polytomous IRT models)}

\description{
 This command returns the estimated standard error of the ability estimate, for a given response pattern and a given matrix of item parameters, either under the 4PL model or any suitable polytomous IRT model.
 }

\usage{
eapSem(thEst, it, x, model = NULL, D = 1, priorDist = "norm", 
 	priorPar = c(0, 1), lower = -4, upper = 4, nqp = 33)
 }

\arguments{
 \item{thEst}{numeric: the EAP ability estimate.}
 \item{it}{numeric: a suitable matrix of item parameters. See \bold{Details}.}
 \item{x}{numeric: a vector of item responses.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952). Ignored if \code{model} is not \code{NULL}.} 
 \item{priorDist}{character: specifies the prior distribution. Possible values are \code{"norm"} (default), \code{"unif"} and \code{"Jeffreys"}.}
 \item{priorPar}{numeric: vector of two components specifying the prior parameters (default is \code{c(0,1)}). Ignored if \code{priorDist="Jeffreys"}.
               See \bold{Details}.}
 \item{lower}{numeric: the lower bound for numercal integration (default is -4).}
 \item{upper}{numeric: the upper bound for numercal integration (default is 4).}
 \item{nqp}{numeric: the number of quadrature points (default is 33).}
 }

\value{
 The estimated standard error of the EAP ability level.
 }
 
\details{
 This command computes the standard error of the EAP (expected a posteriori) ability estimator (Bock and Mislevy, 1982).

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{it} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model
 (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model, \code{"PCM"} for Partical Credit Model, \code{"GPCM"} for Generalized Partial Credit Model, \code{"RSM"} for Rating Scale Model and \code{"NRM"} for Nominal Response Model. The \code{it} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.

 Three prior distributions are available: the normal distribution, the uniform distribution and Jeffreys' prior distribution (Jeffreys, 1939, 1946).
 The prior distribution is specified by the argument \code{priorPar}, with values \code{"norm"}, \code{"unif"} and \code{"Jeffreys"}, respectively. 
 
 The argument \code{priorPar} determines either the prior mean and standard deviation of the normal prior distribution (if 
 \code{priorDist="norm"}), or the range for defining the prior uniform distribution (if \code{priorDist="unif"}). This argument is ignored if \code{priorDist="Jeffreys"}.
 
 The required integrals are approximated by numerical adaptive quadrature. This is achieved by using the \code{\link{integrate.catR}} function. Arguments \code{lower}, \code{upper} and \code{nqp} define respectively the lower and upper bounds for numerical integration, and the number
 of quadrature points. By default, the numerical integration runs with 33 quadrature points on the range [-4; 4], that is, a sequence of values from -4 to 4 by steps of 0.25.
 
 Note that in the current version, the EAP ability estimate must be specified through the \code{thEst} argument.
 }

\references{
  Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}. Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

Bock, R. D., and Mislevy, R. J. (1982). Adaptive EAP estimation of ability in a microcomputer environment. \emph{Applied Psychological Measurement, 6}, 431-444. \doi{10.1177/014662168200600405}

 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error}. Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 
 
 Jeffreys, H. (1939). \emph{Theory of probability}. Oxford, UK: Oxford University Press.

 Jeffreys, H. (1946). An invariant form for the prior probability in estimation problems. \emph{Proceedings of the Royal Society of London. Series A, Mathematical and Physical Sciences, 186}, 453-461.

Magis, D. and Barrada, J. R. (2017). Computerized Adaptive Testing with R: Recent Updates of the Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{Code Snippets}, \emph{76(1)}, 1-18. \doi{10.18637/jss.v076.c01}

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. \doi{10.18637/jss.v048.i08}
}

\author{
    David Magis \cr
    Department of Psychology, University of Liege, Belgium \cr
    \email{david.magis@uliege.be}
 }

\seealso{
\code{\link{thetaEst}}, \code{\link{genPolyMatrix}}, \code{\link{integrate.catR}}
} 

\examples{

## Dichotomous models ##

 # Loading the 'tcals' parameters 
 data(tcals)

 # Selecting item parameters only
 tcals <- as.matrix(tcals[,1:4])

 # Creation of a response pattern (tcals item parameters,
 # true ability level 0)
 set.seed(1)
 x <- genPattern(0, tcals)

 # EAP estimation, standard normal prior distribution
 th <- eapEst(tcals, x)
 c(th, eapSem(th, tcals, x))

 # EAP estimation, uniform prior distribution upon range [-2,2]
 th <- eapEst(tcals, x, priorDist = "unif", priorPar = c(-2, 2))
 c(th, eapSem(th, tcals, x, priorDist = "unif", priorPar=c(-2, 2)))

 # EAP estimation, Jeffreys' prior distribution  
 th <- eapEst(tcals, x, priorDist = "Jeffreys")
 c(th, eapSem(th, tcals, x, priorDist = "Jeffreys"))

\dontrun{

## Polytomous models ##

 # Generation of an item bank under GRM with 100 items and at most 4 categories
 m.GRM <- genPolyMatrix(100, 4, "GRM")
 m.GRM <- as.matrix(m.GRM)

 # Creation of a response pattern (true ability level 0)
 set.seed(1)
 x <- genPattern(0, m.GRM, model = "GRM")

 # EAP estimation, standard normal prior distribution
 th <- eapEst(m.GRM, x, model = "GRM")
 c(th, eapSem(th, m.GRM, x, model = "GRM"))

 # EAP estimation, uniform prior distribution upon range [-2,2]
 th <- eapEst(m.GRM, x, model = "GRM", priorDist = "unif", priorPar = c(-2, 2))
 c(th, eapSem(th, m.GRM, x, model = "GRM", priorDist = "unif", priorPar = c(-2, 2)))

 # EAP estimation, Jeffreys' prior distribution  
 th <- eapEst(m.GRM, x, model = "GRM", priorDist = "Jeffreys")
 c(th, eapSem(th, m.GRM, x, model = "GRM", priorDist = "Jeffreys"))


# Loading the cat_pav data
 data(cat_pav)
 cat_pav <- as.matrix(cat_pav)

 # Creation of a response pattern (true ability level 0)
 set.seed(1)
 x <- genPattern(0, cat_pav, model = "GPCM")

 # EAP estimation, standard normal prior distribution
 th <- eapEst(cat_pav, x, model = "GPCM")
 c(th, eapSem(th, cat_pav, x, model = "GPCM"))

 # EAP estimation, uniform prior distribution upon range [-2,2]
 th <- eapEst(cat_pav, x, model = "GPCM", priorDist = "unif", priorPar = c(-2, 2))
 c(th, eapSem(th, cat_pav, x, model = "GPCM", priorDist = "unif", priorPar = c(-2, 2)))

 # EAP estimation, Jeffreys' prior distribution  
 th <- eapEst(cat_pav, x, model = "GPCM", priorDist = "Jeffreys")
 c(th, eapSem(th, cat_pav, x, model = "GPCM", priorDist = "Jeffreys"))
 }
}

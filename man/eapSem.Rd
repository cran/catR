\name{eapSem}
\alias{eapSem}

\title{Standard error of EAP ability estimation under the 4PL model}

\description{
 This command returns the estimated standard error of the ability estimate, for a given matrix of
 item parameters of the 4PL model, an ability estimate and a specified estimator. 
 }

\usage{
 eapSem(thEst, it, x, D=1, priorDist="norm", priorPar=c(0,1),
 	lower=-4, upper=4, nqp=33)
 }

\arguments{
 \item{thEst}{numeric: the EAP ability estimate.}
 \item{it}{numeric: a matrix with one row per item and four columns, with the values of the discrimination, the difficulty,
           the pseudo-guessing and the inattention parameters (in this order).}
 \item{x}{numeric: a vector of dichotomous item responses.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952).} 
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

 Three prior distributions are available: the normal distribution, the uniform distribution and Jeffreys' prior distribution (Jeffreys, 1939, 1946).
 The prior distribution is specified by the argument \code{priorPar}, with values \code{"norm"}, \code{"unif"} and \code{"Jeffreys"}, respectively. 
 
 The argument \code{priorPar} determines either the prior mean and standard deviation of the normal prior distribution (if 
 \code{priorDist="norm"}), or the range for defining the prior uniform distribution (if \code{priorDist="unif"}). This argument
 is ignored if \code{priorDist="Jeffreys"}.
 
 The required integrals are approximated by numerical adaptive quadrature. This is achieved by using the \code{\link{integrate.xy}} function of the package \code{sfsmisc}.
 Arguments \code{lower}, \code{upper} and \code{nqp} define respectively the lower and upper bounds for numerical integration, and the number
 of quadrature points. By default, the numerical integration runs with 33 quadrature points on the range [-4; 4], that is, a sequence of values from 
 -4 to 4 by steps of 0.25.
 
 Note that in the current version, the EAP ability estimate must be specified through the \code{thEst} argument.
 }

\references{
 Bock, R. D., and Mislevy, R. J. (1982). Adaptive EAP estimation of ability in a microcomputer environment. 
 \emph{Applied Psychological Measurement, 6}, 431-444.

 Haley, D.C. (1952). Estimation of the dosage mortality relationship when the dose is subject to error. 
 Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 
 
 Jeffreys, H. (1939). \emph{Theory of probability}. Oxford, UK: Oxford University Press.

 Jeffreys, H. (1946). An invariant form for the prior probability in estimation problems. \emph{Proceedings of the Royal
 Society of London. Series A, Mathematical and Physical Sciences, 186}, 453-461.

 Magis, D., and Raiche, G. (in press). Random generation of response patterns under computerized adaptive testing
 with the R package \emph{catR}. \emph{Journal of Statistical Software}.
}

\author{
    David Magis \cr
    Post-doctoral researcher, FNRS (Fonds National de la Recherche Scientifique) \cr
    Department of Mathematics, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{thetaEst}}, \code{\link{integrate.xy}}
} 

\examples{
\dontrun{

 # Loading the 'tcals' parameters 
 data(tcals)

 # Selecting item parameters only
 tcals <- as.matrix(tcals[,1:4])

 # Creation of a response pattern (tcals item parameters,
 # true ability level 0)
 set.seed(1)
 x <- rbinom(85, 1, Pi(0, tcals)$Pi)

 # EAP estimation, standard normal prior distribution
 th <- eapEst(tcals, x)
 c(th, eapSem(th, tcals, x))

 # EAP estimation, uniform prior distribution upon range [-2,2]
 th <- eapEst(tcals, x, priorDist="unif", priorPar=c(-2,2))
 c(th, eapSem(th, tcals, x, priorDist="unif", priorPar=c(-2,2)))

 # EAP estimation, Jeffreys' prior distribution  
 th <- eapEst(tcals, x, priorDist="Jeffreys")
 c(th, eapSem(th, tcals, x, priorDist="Jeffreys"))
 }
 }

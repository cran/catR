\name{eapEst}
\alias{eapEst}

\title{EAP ability estimation under the 4PL model}

\description{
 This command returns the EAP (expected a posteriori) ability estimate for a given matrix of item parameters of the 4PL model
 and a given response pattern.
 }

\usage{
 eapEst(it, x, D=1, priorDist="norm", priorPar=c(0,1),
 	lower=-4, upper=4, nqp=20)
 }

\arguments{
 \item{it}{numeric: a matrix with one row per item and four columns, with the values of the discrimination, the difficulty,
           the pseudo-guessing and the inattention parameters (in this order).}
 \item{x}{numeric: a vector of dichotomous item responses.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952).} 
 \item{priorDist}{character: specifies the prior distribution. Possible values are \code{"norm"} (default), \code{"unif"} and \code{"Jeffreys"}.}
 \item{priorPar}{numeric: vector of two components specifying the prior parameters (default is \code{c(0,1)}). Ignored if \code{priorDist="Jeffreys"}.
               See \bold{Details}.}
 \item{lower}{numeric: the lower bound for numercal integration (default is -4).}
 \item{upper}{numeric: the upper bound for numercal integration (default is 4).}
 \item{nqp}{numeric: the number of quadrature points (default is 20).}
}

\value{
 The estimated EAP ability level.
 }
 
\details{
 The EAP (expected a posteriori) ability estimator (Bock and Mislevy, 1982) is obtained by computing the average of the posterior distribution
 of ability, set as the prior distribution times the likelihood function. 

 Three prior distributions are available: the normal distribution, the uniform distribution and Jeffreys' prior distribution (Jeffreys, 1939, 1946).
 The prior distribution is specified by the argument \code{priorPar}, with values \code{"norm"}, \code{"unif"} and \code{"Jeffreys"}, respectively. 
 
 The argument \code{priorPar} determines either the prior mean and standard deviation of the normal prior distribution (if 
 \code{priorDist="norm"}), or the range for defining the prior uniform distribution (if \code{priorDist="unif"}). This argument
 is ignored if \code{priorDist="Jeffreys"}.
 
 The required integrals are approximated by numerical adaptive quadarture using the \code{\link{integrate}} function. Arguments \code{lower},
 \code{upper} and \code{nqp} define respectively the lower and upper bounds for numerical integration, and the number of quadrature points. By
 default, the numerical integration runs with 20 quadrature points on the range [-4; 4].
 }

\references{
 Bock, R. D., and Mislevy, R. J. (1982). Adaptive EAP estimation of ability in a microcomputer environment. 
 \emph{Applied Psychological Measurement, 6}, 431-444.

 Haley, D.C. (1952). Estimation of the dosage mortality relationship when the dose is subject to error. 
 Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 
 
 Jeffreys, H. (1939). \emph{Theory of probability}. Oxford, UK: Oxford University Press.

 Jeffreys, H. (1946). An invariant form for the prior probability in estimation problems. \emph{Proceedings of the Royal
 Society of London. Series A, Mathematical and Physical Sciences, 186}, 453-461.
}

\author{
    David Magis \cr
    Post-doctoral researcher, FNRS (Fonds National de la Recherche Scientifique) \cr
    Department of Mathematics, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{thetaEst}}
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

 # EAP estimation, standard normal prior distribution
 eapEst(tcals, x)

 # EAP estimation, uniform prior distribution upon range [-2,2]
 eapEst(tcals, x, priorDist="unif", priorPar=c(-2,2))

 # EAP estimation, Jeffreys' prior distribution  
 eapEst(tcals, x, priorDist="Jeffreys")

 # Changing the integration settings
 eapEst(tcals, x, nqp=100)
 }
 }

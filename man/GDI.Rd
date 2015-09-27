\name{GDI}
\alias{GDI}

\title{Global-discrimination index (GDI) and posterior global-discrimination index (GDIP) for item selection}

\description{
 This command returns the value of the global-discrimination index (GDI) and posterior global-discrimination index (GDIP) for a given item, an item bank and a set of previously administered items.}

\usage{
GDI(itemBank, item, x, it.given, model = NULL, lower = -4, upper = 4, nqp = 33, 
 	type = "GDI", priorDist="norm", priorPar = c(0, 1), D = 1, X = NULL, 
 	lik = NULL)}

\arguments{
 \item{itemBank}{numeric: a suitable matrix of item parameters. See \bold{Details}.}
 \item{item}{numeric: the item (referred to as its rank in the item bank) for which the GDI or GDIP must be computed.}
 \item{x}{numeric: a vector of item responses, coded as 0 or 1 only (for dichotomous items) or from 0 to the number of response categories minus one (for polytomous items).}
 \item{it.given}{numeric: a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). The number of rows of \code{it} must be equal to the length of \code{x}.}
\item{model}{now only \code{NULL} (default) for dichotomous models is allowed.}
 \item{lower}{numeric: the lower bound for numercal integration (default is -4).}
 \item{upper}{numeric: the upper bound for numercal integration (default is 4).}
 \item{nqp}{numeric: the number of quadrature points (default is 33).}
 \item{type}{character: the type of index to be computed. Possible values are \code{"GDI"} (default) and \code{"GDIP"}. See \bold{Details}.}
 \item{priorDist}{character: the prior ability distribution. Possible values are \code{"norm"} (default) for the normal distribution, and \code{"unif"} for the uniform distribution. Ignored if \code{type} is \code{"GDI"}.}
 \item{priorPar}{numeric: a vector of two components with the prior parameters. If \code{priorDist} is \code{"norm"}, then \code{priorPar} contains the mean and the standard deviation of the normal distribution. If \code{priorDist} is \code{"unif"}, then \code{priorPar} contains the bounds of the uniform distribution. The default values are 0 and 1 respectively. Ignored if \code{type} is \code{"GDI"}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952).} 
\item{X}{either a vector of numeric values or \code{NULL} (default). See \bold{Details}.}
\item{lik}{either a vector of numeric values or \code{NULL} (default). See \bold{Details}.}
}

\value{
 The required GDI or GDIP value for the selected item.
 }
 
\details{
 Global-discrimination index  can be used as a rule for selecting the next item in the CAT process (Kaplan, de la Torre, and Barrada, 2015). This command serves as a subroutine for the \code{\link{nextItem}} function.

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{itemBank} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model (Barton and Lord, 1981). 

Currently both GDI and GDIP are not implemented for polytomous IRT models.

The integrals within GDI and GDIP are approximated by the \code{\link{integrate.catR}} function. The range of integration is set up by the arguments \code{lower}, \code{upper} and \code{nqp}, giving respectively the lower bound, the upper bound and the number of quadrature points. The default range goes from -4 to 4 with length 33 (that is, by steps of 0.25).

To speed up the computation, both the range of integration of values \eqn{\theta} and the values of the likelihood function \eqn{L(\theta)} can be directly provided to the function through the arguments \code{X} and \code{lik}. If \code{X} is set to \code{NULL} (default), the sequence of ability values for integration is determined by the arguments \code{lower}, \code{upper} and \code{nqp} as explained above. If \code{lik} is \code{NULL} (default), it is also internally computed from an implementation of the likelihood function.

The provisional response pattern and the related item parameters are provided by the arguments \code{x} and \code{it.given} respectively. The target item (for which the KL information is computed) is given by its rank number in the item bank, through the \code{item} argument. 

 The argument \code{type} defines the type of KL information to be computed. The default value, \code{"GDI"}, computes the GDI indexinformation, while the posterior GDI index is obtained with \code{type="GDIP"}. For the latter, the \code{priorDist} and \code{priorPar} arguments fix the prior ability distribution. The normal distribution is set up by \code{priorDist="norm"} and then, \code{priorPar} contains the mean and the standard deviation of the normal distribution. If \code{priorDist} is \code{"unif"}, then the uniform distribution is considered, and \code{priorPar} fixes the lower and upper bounds of that uniform distribution. By default, the standard normal prior distribution is assumed. These arguments are ignored whenever \code{method} is \code{"GDI"}.
}

\references{
 Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}.
 Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error}. Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 

Kaplan, M., de la Torre, J., and Barrada, J. R. (2015). New item selection methods for cognitive diagnosis computerized adaptive testing. \emph{Applied Psychological Measurement, 39}, 167-188. doi: 10.1177/0146621614554650

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R 
 Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. URL \url{http://www.jstatsoft.org/v48/i08/}
}

\author{
    David Magis \cr
    Department of Education, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}

 Juan Ramon Barrada \cr
Department of Psychology and Sociology, Universidad Zaragoza, Spain \cr
\email{barrada@unizar.es}
 }

\seealso{
\code{\link{integrate.catR}}, \code{\link{nextItem}}, \code{\link{genPolyMatrix}} 
} 

\examples{

## Dichotomous models ##

 # Loading the 'tcals' parameters 
 data(tcals)

 # Selecting item parameters only
 bank <- as.matrix(tcals[,1:4])
 
 # Selection of two arbitrary items (15 and 20) of the
 # 'tcals' data set
 it.given <- bank[c(15, 20),]

 # Creation of a response pattern
 x <- c(0, 1)

 # GDI for item 1
 GDI(bank, 1, x, it.given)

 # GDIP for item 1
 GDI(bank, 1, x, it.given, type = "GDIP")

 # GDIP for item 1, different integration range
 GDI(bank, 1, x, it.given, type = "GDIP", lower = -2, upper = 2, nqp = 20)

 # GDIP for item 1, uniform prior distribution on the range [-2,2]
 GDI(bank, 1, x, it.given, type = "GDIP", priorDist = "unif", 
    priorPar = c(-2, 2))

 # Computation of likelihood function beforehand
 L <- function(th, r, param) 
  prod(Pi(th, param)$Pi^r * (1 - Pi(th,param)$Pi)^(1 - r))
 xx <- seq(from = -4, to = 4, length = 33)
 y <- sapply(xx, L, x, it.given) 
 GDI(bank, 1, x, it.given, X = xx, lik = y)
 }

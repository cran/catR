\name{KL}
\alias{KL}

\title{Kullback-Leibler (KL) and posterior Kullback-Leibler (KLP) values for item selection}

\description{
 This command returns the value of the Kullback-Leibler (KL) or posterior Kullback-Leibler (KLP) for a given item, an item bank and a set of previously administered items.
 }

\usage{
KL(itemBank, item, x, it.given, model = NULL, theta = NULL, lower = -4, 
  upper = 4, nqp = 33, type = "KL", priorDist = "norm", priorPar = c(0, 1), 
  D = 1, X = NULL, lik = NULL) 
 }

\arguments{
 \item{itemBank}{numeric: a suitable matrix of item parameters. See \bold{Details}.}
 \item{item}{numeric: the item (referred to as its rank in the item bank) for which the KL or KLP must be computed.}
 \item{x}{numeric: a vector of item responses, coded as 0 or 1 only (for dichotomous items) or from 0 to the number of response categories minus one (for polytomous items).}
 \item{it.given}{numeric: a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). The number of rows of \code{it} must be equal to the length of \code{x}.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
\item{theta}{either a numeric value for provisional ability estimate or \code{NULL}. See \bold{Details}.}
 \item{lower}{numeric: the lower bound for numercal integration (default is -4).}
 \item{upper}{numeric: the upper bound for numercal integration (default is 4).}
 \item{nqp}{numeric: the number of quadrature points (default is 33).}
 \item{type}{character: the type of information to be computed. Possible values are \code{"KL"} (default) and \code{"KLP"}. See \bold{Details}.}
 \item{priorDist}{character: the prior ability distribution. Possible values are \code{"norm"} (default) for the normal distribution, and \code{"unif"} for the uniform distribution. Ignored if \code{type} is \code{"KL"}.}
 \item{priorPar}{numeric: a vector of two components with the prior parameters. If \code{priorDist} is \code{"norm"}, then \code{priorPar} contains the mean and the standard deviation of the normal distribution. If \code{priorDist} is \code{"unif"}, then \code{priorPar} contains the bounds of the uniform distribution. The default values are 0 and 1 respectively. Ignored if \code{type} is \code{"KL"}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952).} 
\item{X}{either a vector of numeric values or \code{NULL} (default). See \bold{Details}.}
\item{lik}{either a vector of numeric values or \code{NULL} (default). See \bold{Details}.}
}

\value{
 The required KL or KLP value for the selected item.
 }
 
\details{
 Kullback-Leibler information can be used as a rule for selecting the next item in the CAT process (Barrada, Olea, Ponsoda and Abad, 2010; Chang and Ying, 1996), both with dichotomous and polytomous IRT models. This command serves as a subroutine for the \code{\link{nextItem}} function.

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{itemBank} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model, \code{"PCM"} for Partical Credit Model, \code{"GPCM"} for Generalized Partial Credit Model, \code{"RSM"} for Rating Scale Model and \code{"NRM"} for Nominal Response Model. The \code{itemBank} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.

Under polytomous IRT models, let \emph{k} be the number of administered items, and set \eqn{x_1, ..., x_k} as the provisional response pattern (where each response \eqn{x_l} takes values in \eqn{\{0, 1, ..., g_l\}}). Set \eqn{\hat{\theta}_k} as the provisional ability estimate (with the first \emph{k} responses) and let \emph{j} be the item of interest (not previously administered). Set also \eqn{L(\theta | x_1, ..., x_k) } as the likelihood function of the first \eqn{k} items and evaluated at \eqn{\theta}. Set finally \eqn{P_{jt}(\theta)} as the probability of answering response category \emph{t} to item \emph{j} for a given ability level \eqn{\theta}. Then, Kullack-Leibler (KL) information is defined as
\deqn{KL_j(\theta || \hat{\theta}_k) = \sum_{t=0}^{g_j} \,P_{jt}(\hat{\theta}_k) \,\log \left( \frac{P_{jt}(\hat{\theta}_k)}{P_{jt}(\theta)}\right).}

In case of dichotomous IRT models, all \eqn{g_l} values reduce to 1, so that item responses \eqn{x_l} equal either 0 or 1. Set simply \eqn{P_j(\theta)} as the probability of answering item \emph{j} correctly for a given ability level \eqn{\theta}. Then, KL information reduces to
\deqn{KL_j(\theta || \hat{\theta}) = P_j(\hat{\theta}) \,\log \left( \frac{P_j(\hat{\theta}_k)}{P_j(\theta)}\right) + [1-P_j(\hat{\theta}_k)] \,\log \left( \frac{1-P_j(\hat{\theta}_k)}{1-P_j(\theta)}\right).}

The quantity that is returned by this \code{KL} function is either: the likelihood function weighted by Kullback-Leibler information (the KL value):
\deqn{KL_j(\hat{\theta}_k) = \int KL_j(\theta || \hat{\theta}_k) \, L(\theta | x_1, ..., x_k) \,d\theta}
or the posterior function weighted by Kullback-Leibler information (the KLP value):
\deqn{KLP_j(\hat{\theta}) = \int KL_j(\theta || \hat{\theta}_k) \, \pi(\theta) \,L(\theta | x_1, ..., x_k) \,d\theta}
 where \eqn{\pi(\theta)} is the prior distribution of the ability level.

 These integrals are approximated by the \code{\link{integrate.catR}} function. The range of integration is set up 
 by the arguments \code{lower}, \code{upper} and \code{nqp}, giving respectively the lower bound, the upper bound and the number of quadrature points. The default range goes from -4 to 4 with length 33 (that is, by steps of 0.25).

To speed up the computation, both the range of integration of values \eqn{\theta} and the values of the likelihood function \eqn{L(\theta)} can be directly provided to the function through the arguments \code{X} and \code{lik}. If \code{X} is set to \code{NULL} (default), the sequence of ability values for integration is determined by the arguments \code{lower}, \code{upper} and \code{nqp} as explained above. If \code{lik} is \code{NULL} (default), it is also internally computed from an implementation of the likelihood function.

The provisional response pattern and the related item parameters are provided by the arguments \code{x} and \code{it.given} respectively. The target item (for which the KL information is computed) is given by its rank number in the item bank, through the \code{item} argument. 

An ability level estimate must be provided to compute KL and KLP information values. Either the value is specified through the \code{theta} argument, or it is left equal to \code{NULL} (default). In this case, ability estimate is computed internally by maximum likelihood, using the \code{\link{thetaEst}} function with arguments \code{it.given} and \code{x}.

Note that the provisional response pattern \code{x} can also be set to \code{NULL} (which is useful when the number \code{nrItems} of starting items is set to zero). In this case, \code{it.given} must be a matrix with zero rows, such as e.g., \code{itemBank[NULL,]}. In this very specific configuration, the likelihood function \eqn{L(\theta | x_1, ..., x_k)} reduces to the constant value 1 on the whole \eqn{\theta} range (that is, a uniform likelihood). 

 The argument \code{type} defines the type of KL information to be computed. The default value, \code{"KL"}, computes the usual Kullback-Leibler information, while the posterior Kullback-Leibler value is obtained with \code{type="KLP"}. For the latter, the \code{priorDist} and \code{priorPar} arguments fix the prior ability distribution. The normal distribution is set up by \code{priorDist="norm"} and then, \code{priorPar} contains the mean and the standard deviation of the normal distribution. If \code{priorDist} is \code{"unif"}, then the uniform distribution is considered, and \code{priorPar} fixes the lower and upper bounds of that uniform distribution. By default, the standard normal prior distribution is assumed. These arguments are ignored whenever \code{method} is \code{"KL"}.
}

\references{
 Barrada, J. R., Olea, J., Ponsoda, V., and Abad, F. J. (2010). A method for the comparison of item selection rules in computerized adaptive testing. \emph{Applied Psychological Measurement}, \emph{20}, 213-229. \doi{10.1177/0146621610370152}

 Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}.
 Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Chang, H.-H., and Ying, Z. (1996). A global information approach to computerized adaptive testing. \emph{Applied Psychological Measurement, 34}, 438-452. \doi{10.1177/014662169602000303}

 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error}. Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 

Magis, D. and Barrada, J. R. (2017). Computerized Adaptive Testing with R: Recent Updates of the Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{Code Snippets}, \emph{76(1)}, 1-18. \doi{10.18637/jss.v076.c01}

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. \doi{10.18637/jss.v048.i08}
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

 # KL for item 1, ML estimate of ability computed
 KL(bank, 1, x, it.given)

 # Current (ML) ability estimate 
 theta <- thetaEst(it.given, x, method = "ML")
 KL(bank, 1, x, it.given, theta = theta)

 # WL ability estimate instead
 theta <- thetaEst(it.given, x, method = "WL")
 KL(bank, 1, x, it.given, theta = theta)

 # KLP for item 1
 KL(bank, 1, x, it.given, theta = theta, type = "KLP")

 # KLP for item 1, different integration range
 KL(bank, 1, x, it.given, theta = theta, type = "KLP", lower = -2, upper = 2, nqp = 20)

 # KL for item 1, uniform prior distribution on the range [-2,2]
 KL(bank, 1, x, it.given, theta = theta, type = "KLP", priorDist = "unif", 
    priorPar = c(-2, 2))

 # Computation of likelihood function beforehand
 L <- function(th, r, param) 
  prod(Pi(th, param)$Pi^r * (1 - Pi(th,param)$Pi)^(1 - r))
 xx <- seq(from = -4, to = 4, length = 33)
 y <- sapply(xx, L, x, it.given) 
 KL(bank, 1, x, it.given, theta = theta, X = xx, lik = y)


## Polytomous models ##

 # Generation of an item bank under GRM with 100 items and at most 4 categories
 m.GRM <- genPolyMatrix(100, 4, "GRM")
 m.GRM <- as.matrix(m.GRM)

 # Selection of two arbitrary items (15 and 20) 
 it.given <- m.GRM[c(15, 20),]

 # Generation of a response pattern (true ability level 0)
 x <- genPattern(0, it.given, model = "GRM")

 # KL for item 1, ML estimate of ability computed
 KL(m.GRM, 1, x, it.given, model = "GRM")

 # Current (ML) ability estimate 
 theta <- thetaEst(it.given, x, method = "ML", model = "GRM")
 KL(m.GRM, 1, x, it.given, theta = theta, model = "GRM")

 # WL ability estimate instead
 theta <- thetaEst(it.given, x, method = "WL", model = "GRM")
 KL(m.GRM, 1, x, it.given, theta = theta, model = "GRM")

 # KLP for item 1
 KL(m.GRM, 1, x, it.given, theta = theta, model = "GRM", type = "KLP")

 # KLP for item 1, different integration range
 KL(m.GRM, 1, x, it.given, theta = theta, model = "GRM", type = "KLP", lower = -2, 
    upper = 2, nqp = 20)

 # KL for item 1, uniform prior distribution on the range [-2,2]
 KL(m.GRM, 1, x, it.given, theta = theta, model = "GRM", type = "KLP", 
    priorDist = "unif", priorPar = c(-2, 2))


 # Loading the cat_pav data
 data(cat_pav)
 cat_pav <- as.matrix(cat_pav)

 # Selection of two arbitrary items (15 and 20) 
 it.given <- cat_pav[c(15, 20),]

 # Generation of a response pattern (true ability level 0)
 x <- genPattern(0, it.given, model = "GPCM")

  # KL for item 1, ML estimate of ability computed
 KL(cat_pav, 1, x, it.given, model = "GPCM")

 # Current (ML) ability estimate 
 theta <- thetaEst(it.given, x, method = "ML", model = "GPCM")
 KL(cat_pav, 1, x, it.given, theta = theta, model = "GPCM")

 # WL ability estimate instead
 theta <- thetaEst(it.given, x, method = "WL", model = "GPCM")
 KL(cat_pav, 1, x, it.given, theta = theta, model = "GPCM")

 # KLP for item 1
 KL(cat_pav, 1, x, it.given, theta = theta, model = "GPCM", type = "KLP")

 # KLP for item 1, different integration range
 KL(cat_pav, 1, x, it.given, theta = theta, model = "GPCM", type = "KLP", lower = -2, 
    upper = 2, nqp = 20)

 # KL for item 1, uniform prior distribution on the range [-2,2]
 KL(cat_pav, 1, x, it.given, theta = theta, model = "GPCM", type = "KLP", 
    priorDist = "unif", priorPar = c(-2, 2))
 }

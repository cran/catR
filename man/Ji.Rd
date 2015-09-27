\name{Ji}
\alias{Ji}

\title{Function \eqn{J(\theta)} for weighted likelihood estimation (dichotomous and polytomous IRT models)}

\description{
 This command returns the \eqn{J(\theta)} function that is necessary to obtain the weighted likelihood estimation of ability with dichotomous and polytomous IRT models, as well as its asymptotic standard error.
 }

\usage{
Ji(th, it, model = NULL, D = 1)
 }

\arguments{
 \item{th}{numeric: the ability value.}
 \item{it}{numeric: a suitable matrix of item parameters. See \bold{Details}.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952). Ignored if \code{model} is not \code{NULL}.} 
}

\value{
 A list with two arguments:
 \item{Ji}{the vector with \eqn{J(\theta)} values (one value per item).}
 \item{dJi}{the vector with first derivatives of the \eqn{J(\theta)} values (one value per item).}
 }
 
\details{
The \eqn{J(\theta)} fucntion is defined by (Samejima, 1998):
\deqn{J(\theta)=\sum_{j=1}^n \,\sum_{k=0}^{g_j} \frac{P_{jk}'(\theta)\,P_{jk}''(\theta)}{P_{jk}(\theta)}}
where \eqn{n} is the number of items; \eqn{g_j} the number of response categories for item \emph{j} (\eqn{j=1, ..., n}); \eqn{P_{jk}(\theta)} the response category probabilities and \eqn{P_{jk}'(\theta)} and \eqn{P_{jk}''(\theta)} the first and second derivatives with respect to \eqn{\theta}. In case of dichotomous IRT models, this reduces to (Warm, 1989):
\deqn{J(\theta)=\sum_{j=1}^n \frac{P_j'(\theta)\,P_j''(\theta)}{P_j(\theta)\,Q_j(\theta)}}
with \eqn{Q_j(\theta)=1-P_j(\theta)}.

This function is useful to compute the weighted likelihood estimates of ability with dichotomous and polytomous IRT models as well as their related asymptotic standard errors. 

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{it} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model
 (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model, \code{"PCM"} for Partical Credit Model, \code{"GPCM"} for Generalized Partial Credit Model, \code{"RSM"} for Rating Scale Model and \code{"NRM"} for Nominal Response Model. The \code{it} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.
}

\references{
Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}. Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error}. Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 

Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. URL \url{http://www.jstatsoft.org/v48/i08/}

Samejima, F. (1998, April). \emph{Expansion of Warm's weighted likelihood estimator of ability for the three-parameter logistic model to generate discrete responses.} PPaper presented at the annual meeting of the National Council on Measurement in Education, San Diego, CA.

Warm, T.A. (1989). Weighted likelihood estimation of ability in item response models. \emph{Psychometrika, 54}, 427-450. doi: 10.1007/BF02294627
 }

\author{
    David Magis \cr
    Department of Education, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{thetaEst}}, \code{\link{semTheta}}, \code{\link{genPolyMatrix}} 
} 

\examples{

## Dichotomous models ##
 
 # Loading the 'tcals' parameters 
 data(tcals)

 # Selecting item parameters only
 tcals <- as.matrix(tcals[,1:4])
 
 # Various J functions and derivatives
 # (various th and D values)
 Ji(th = 0, tcals)
 Ji(th = 0, tcals, D = 1.702)
 Ji(th = 1, tcals)

## Polytomous models ##

 # Generation of an item bank under GRM with 100 items and at most 4 categories
 m.GRM <- genPolyMatrix(100, 4, "GRM")
 m.GRM <- as.matrix(m.GRM)

 # Computation of J function and derivatives for ability level 0
 Ji(0, m.GRM, model = "GRM")

 # Generation of a item bank under PCM with 20 items and at most 3 categories
 m.PCM <- genPolyMatrix(20, 3, "PCM")
 m.PCM <- as.matrix(m.PCM)

 # Computation of J function and derivatives for ability level 1
 Ji(1, m.PCM, model = "PCM")
 }

\name{EPV}
\alias{EPV}

\title{Expected Posterior Variance (EPV)}

\description{
 This command returns the expected posterior variance (EPV) for a given item under dichotomous and polytomous IRT models, as used for Minimum Expected Posterior Variance (MEPV) criterion.
 }

\usage{
EPV(itemBank, item, x, theta, it.given, model = NULL, priorDist = "norm", 
 	priorPar = c(0, 1), D = 1, parInt = c(-4, 4, 33))
 }

\arguments{
 \item{itemBank}{numeric: a suitable matrix of item parameters. See \bold{Details}.}
 \item{item}{numeric: the item (referred to as its rank in the item bank) for which the expected posterior variance must be computed.}
 \item{x}{numeric: a vector of item responses, coded as 0 or 1 only (for dichotomous items) or from 0 to the number of response categories minus one (for polytomous items).}
 \item{theta}{numeric: the provisional ability estimate.}
 \item{it.given}{numeric: a suitable matrix of item parameters for previously administered items. The number of rows of \code{it.given} must be equal to the length of \code{x}.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
 \item{priorDist}{character: specifies the prior distribution. Possible values are \code{"norm"} (default) and \code{"unif"}.}
 \item{priorPar}{numeric: vector of two components specifying the prior parameters (default is \code{c(0,1)}) of the prior ability distribution.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952). Ignored if \code{model} is not \code{NULL}.} 
 \item{parInt}{numeric: vector of three components, defining the sequence of ability values for computing the posterior variance. See \bold{Details}.}
}

\value{
 The expected posterior variance for the selected item.
 }
 
\details{
 The EPV can be used as a rule for selecting the next item in the CAT process (Choi and Swartz, 2009; Owen, 1975; van der Linden, 1998). This command serves as a subroutine for the \code{\link{nextItem}} function.

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{itemBank} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model, \code{"PCM"} for Partical Credit Model, \code{"GPCM"} for Generalized Partial Credit Model, \code{"RSM"} for Rating Scale Model and \code{"NRM"} for Nominal Response Model. The \code{itemBank} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.

Under polytomous IRT models, let \emph{k} be the number of administered items, and set \eqn{x_1, ..., x_k} as the provisional response pattern (where each response \eqn{x_l} takes values in \eqn{\{0, 1, ..., g_l\}}). Set \eqn{\hat{\theta}_k} as the provisional ability estimate (with the first \emph{k} responses) and let \emph{j} be the item of interest (not previously administered). Set also \eqn{P_{jt}(\theta)}
 as the probability of answering response category \emph{t} to item \emph{j} for a given ability level \eqn{\theta} (thus \eqn{t \in \{0, ..., g_j\}}). Finally, set \eqn{Var(\theta | x_1, ..., x_k, t)} as the posterior variance of \eqn{\theta}, given the provisional response pattern (updated by response \eqn{t}). Then, the EPV for item \emph{j} equals
  \deqn{EPV_j = \sum_{t=0}^{g_j} P_{jt}(\hat{\theta}_k)\,Var(\theta | x_1, ..., x_k, t)}.

In case of dichotomous IRT models, all \eqn{g_l} values reduce to 1, so that item responses \eqn{x_l} equal either 0 or 1. Set simply \eqn{P_j(\theta)} as the probability of answering item \emph{j} correctly for a given ability level \eqn{\theta}, and set \eqn{Q_j(\theta)=1-P_j(\theta)}. Finally, set \eqn{Var(\theta | x_1, ..., x_k, 0)} and \eqn{Var(\theta | x_1, ..., x_k, 1)} as the posterior variances of \eqn{\theta}, given the provisional response pattern (updated by response 0 and 1 respectively). Then, the EPV for item \emph{j} reduces to
 \deqn{EPV_j = P_j(\hat{\theta}_k)\,Var(\theta | x_1, ..., x_k, 1) +  Q_j(\hat{\theta}_k)\,Var(\theta | x_1, ..., x_k, 0)}.

 The posterior variances \eqn{Var(\theta | x_1, ..., x_k, x_j)} (where \eqn{x_j} takes value 0 or 1 for dichotomous models, and 0, 1, ..., or \eqn{g_j} for polytomous models) is computed as the squared standard error of the EAP estimate of ability, using the response pattern \eqn{(x_1, ..., x_k, x_j)}. This is done by a joint use of the \code{\link{eapEst}} and \code{\link{eapSem}} functions. 

 The prior distribution is set up by the arguments \code{priorDist} and \code{priorPar}, with the by-default standard normal distribution. The range of integration is defined by the \code{parInt} argument, with by default, the sequence from -4 to 4 and of length 33 (or, by steps of 0.25). See the function
 \code{\link{eapEst}} for further details. 

The provisional response pattern and the related item parameters are provided by the arguments \code{x} and \code{it.given} respectively. The target item (for which the maximum information computed) is given by its number in the item bank, through the
 \code{item} argument. 

Note that the provisional response pattern \code{x} can also be set to \code{NULL} (which is useful when the number \code{nrItems} of starting items is set to zero). In this case, \code{it.given} must be a matrix with zero rows, such as e.g., \code{itemBank[NULL,]}. 
 }

\references{
 Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}. Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Choi, S. W., and Swartz, R. J. (2009). Comparison of CAT item selection criteria for polytomous items. \emph{Applied Psychological Measurement, 32}, 419-440. doi: 10.1177/0146621608327801

 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error}. Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R 
 Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. URL \url{http://www.jstatsoft.org/v48/i08/}

 Owen, R. J. (1975). A Bayesian sequential procedure for quantal response in the context of adaptive mental testing.  \emph{Journal of the American Statistical Association, 70}, 351-356. doi: 10.1080/ 01621459.1975.10479871

 van der Linden, W. J. (1998). Bayesian item selection criteria for adaptive testing. \emph{Psychometrika, 63}, 201-216. doi: 10.1007/BF02294775
 }

\author{
    David Magis \cr
    Department of Education, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{nextItem}}, \code{\link{eapEst}}, \code{\link{eapSem}}, \code{\link{genPolyMatrix}} 
} 

\examples{
\dontrun{

## Dichotomous models ##

 # Loading the 'tcals' parameters 
 data(tcals)

 # Selecting item parameters only
 bank <- as.matrix(tcals[,1:4])
 
 # Selection of two arbitrary items (15 and 20) of the
 # 'tcals' data set
 it.given <- bank[c(15,20),]

 # Creation of a response pattern
 x <- c(0, 1)

 # EPV for item 1, provisional ability level 0
 EPV(bank, 1, x, 0, it.given)

 # With prior standard deviation 2
 EPV(bank, 1, x, 0, it.given, priorPar = c(0,2))


## Polytomous models ##

 # Generation of an item bank under GRM with 100 items and at most 4 categories
 m.GRM <- genPolyMatrix(100, 4, "GRM")
 m.GRM <- as.matrix(m.GRM)

 # Selection of two arbitrary items (15 and 20) 
 it.given <- m.GRM[c(15,20),]

 # Generation of a response pattern (true ability level 0)
 x <- genPattern(0, it.given, model = "GRM")

 # EPV for item 1, provisional ability level 0
 EPV(m.GRM, 1, x, 0, it.given, model = "GRM")

 # With prior standard deviation 2
 EPV(m.GRM, 1, x, 0, it.given, model = "GRM", priorPar = c(0, 2))


 # Generation of an item bank under PCM with 20 items and 4 categories
 m.PCM <- genPolyMatrix(20, 4, "PCM", same.nrCat = TRUE)
 m.PCM <- as.matrix(m.PCM)

 # Selection of two arbitrary items (15 and 20) 
 it.given <- m.PCM[c(15, 20),]

 # Generation of a response pattern (true ability level 0)
 x <- genPattern(0, it.given, model = "PCM")

 # EPV for item 1, provisional ability level 0
 EPV(m.PCM, 1, x, 0, it.given, model = "PCM")

 # With prior standard deviation 2
 EPV(m.PCM, 1, x, 0, it.given, model = "PCM", priorPar = c(0, 2))
 }
 }

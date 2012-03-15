\name{EPV}
\alias{EPV}

\title{Expected Posterior Variance (EPV)}

\description{
 This command returns the expected posterior variance (EPV) for a given item, as used for Minimum Expected Posterior
 Variance (MEPV) criterion.
 }

\usage{
 EPV(itemBank, item, x, theta, it, priorDist="norm", 
 	priorPar=c(0,1), D=1, parInt=c(-4,4,33))
 }

\arguments{
 \item{itemBank}{an item bank, i.e. a list of class \code{itBank}, typically an output of the function \code{createItemBank}.}
 \item{item}{numeric: the item (referred to as its rank in the item bank) for which the maximum information must be computed.}
 \item{x}{binary: a vector of item responses, coded as 0 or 1 only.}
 \item{theta}{numeric: the provisional ability estimate.}
 \item{it}{numeric: a matrix with one row per item and four columns, with the values of the discrimination, the difficulty,
           the pseudo-guessing and the inattention parameters (in this order). The number of rows of \code{it} must be equal to the
           length of \code{x}.}
 \item{priorDist}{character: specifies the prior distribution. Possible values are \code{"norm"} (default) and \code{"unif"}.}
 \item{priorPar}{numeric: vector of two components specifying the prior parameters (default is \code{c(0,1)}) of the prior ability distribution.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952).} 
 \item{parInt}{numeric: vector of three components, defining the sequence of ability values for computing the posterior variance. See \bold{Details}.}
}

\value{
 The expected posterior variance for the selected item.
 }
 
\details{
 The EPV can be used as a rule for selecting the next item in the CAT process (Choi and Swartz, 2009; Owen, 1975; van der Linden, 1998). This command serves
 as a subroutine for the \code{\link{nextItem}} function.

 Let \emph{k} be the number of administered items, and set \eqn{x_1, ..., x_k} as the provisional response pattern. Set \eqn{\hat{\theta}_k} as the 
 provisional ability estimate (with the first \emph{k} responses) and let \emph{j} be the item of interest (not previously administered). Set also \eqn{P_j(\theta)}
 as the probability of answering item \emph{j} correctly for a given ability level \eqn{\theta}, and set \eqn{Q_j(\theta)=1-P_j(\theta)}. Finally, set
 \eqn{Var(\theta | x_1, ..., x_k, 0)} and \eqn{Var(\theta | x_1, ..., x_k, 1)} as the posterior variances of \eqn{\theta}, given the provisional response
 pattern (updated by response 0 and 1 respectively). Then, the EPV for item \emph{j} equals
 \deqn{EPV_j = P_j(\hat{\theta}_k)\,Var(\theta | x_1, ..., x_k, 1) +  Q_j(\hat{\theta}_k)\,Var(\theta | x_1, ..., x_k, 0)}.
 
 The posterior variance \eqn{Var(\theta | x_1, ..., x_k, x_j)} (where \code{x_j} takes value 0 or 1) is computed as the squared standard error of the EAP estimate
 of ability, using the response pattern \eqn{(x_1, ..., x_k, x_j)}. This is done by a joint use of the \code{\link{eapEst}} and \code{\link{eapSem}} functions. 

 The prior distribution is set up by the arguments \code{priorDist} and \code{priorPar}, with the by-default standard normal distribution. The range of
 integration is defined by the \code{parInt} argument, with by default, the sequence from -4 to 4 and of length 33 (or, by steps of 0.25). See the function
 \code{\link{eapEst}} for further details. 

 The item bank is provided through the argument \code{itemBank}. The provisional response pattern and the related item parameters are provided by the arguments
 \code{x} and \code{it} respectively. The target item (for which the maximum information computed) is given by its number in the item bank, through the
 \code{item} argument. 
 }

\references{
 Choi, S. W., and Swartz, R. J. (2009). Comparison of CAT item selection criteria for polytomous items. \emph{Applied PScyhological Measurement, 32}, 419-440.

 Haley, D.C. (1952). Estimation of the dosage mortality relationship when the dose is subject to error. 
 Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 

 Magis, D., and Raiche, G. (in press). Random generation of response patterns under computerized adaptive testing
 with the R package \emph{catR}. \emph{Journal of Statistical Software}.

 Owen, R. J. (1975). A Bayesian sequential procedure for quantal response in the context of adaptive mental testing.  \emph{Journal of the American Statistical
 Association, 70}, 351-356.

 van der Linden, W. J. (1998). Bayesian item selection criteria for adaptive testing. \emph{Psychometrika, 63}, 201-216.
 }

\author{
    David Magis \cr
    Post-doctoral researcher, FNRS (Fonds National de la Recherche Scientifique) \cr
    Department of Mathematics, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{nextItem}}, \code{\link{eapEst}}, \code{\link{eapSem}}
} 

\examples{
\dontrun{

 # Loading the 'tcals' parameters 
 data(tcals)

 # Selecting item parameters only
 tcals <- as.matrix(tcals[,1:4])
 
 # Item bank creation with 'tcals' item parameters
 bank <- createItemBank(tcals)

 # Selection of two arbitrary items (15 and 20) of the
 # 'tcals' data set
 it <- bank$itemPar[c(15,20),]

 # Creation of a response pattern
 x <- c(0,1)

 # MEI for item 1, provisional ability level 0
 EPV(bank, 1, x, 0, it)

 # With prior standard deviation 2
 EPV(bank, 1, x, 0, it, priorPar=c(0,2))
 }
 }

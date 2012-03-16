\name{MEI}
\alias{MEI}

\title{(Maximum) Expected Information (MEI)}

\description{
 This command returns the expected information (EI) for a given item, as used for Maximum Expected Information (MEI) criterion.
 }

\usage{
 MEI(itemBank, item, x, theta, it, method="BM", priorDist="norm",
 	priorPar=c(0,1), D=1, range=c(-4,4), parInt=c(-4,4,33),
 	infoType="observed")
 }

\arguments{
 \item{itemBank}{an item bank, i.e. a list of class \code{itBank}, typically an output of the function \code{createItemBank}.}
 \item{item}{numeric: the item (referred to as its rank in the item bank) for which the maximum information must be computed.}
 \item{x}{binary: a vector of item responses, coded as 0 or 1 only.}
 \item{theta}{numeric: the provisional ability estimate.}
 \item{it}{numeric: a matrix with one row per item and four columns, with the values of the discrimination, the difficulty,
           the pseudo-guessing and the inattention parameters (in this order). The number of rows of \code{it} must be equal to the
           length of \code{x}.}
 \item{method}{character: the ability estimator. Possible values are \code{"BM"} (default), \code{"ML"} and \code{"WL"}. See \bold{Details}.}
 \item{priorDist}{character: specifies the prior distribution. Possible values are \code{"norm"} (default), \code{"unif"} and \code{"Jeffreys"}.
 Ignored if \code{method} is neither \code{"BM"} nor \code{"EAP"}. See \bold{Details}.}
 \item{priorPar}{numeric: vector of two components specifying the prior parameters (default is \code{c(0,1)}) of the prior ability distribution. Ignored if \code{method} is neither \code{"BM"}
 nor \code{"EAP"}, or if \code{priorDist="Jeffreys"}. See \bold{Details}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952).} 
 \item{range}{numeric: vector of two components specifying the range wherein the ability estimate must be looked for (default is \code{c(-4,4)}).
 	Ignored if \code{method=="EAP"}.}
 \item{parInt}{numeric: vector of three components, holding respectively the values of the arguments \code{lower}, \code{upper} and \code{nqp} of the
 \code{\link{eapEst}} command. Default vector is (-4, 4, 33). Ignored if \code{method} is not \code{"EAP"}.}
 \item{infoType}{character: the type of information function to be used. Possible values are \code{"observed"} (default) for observed
                 information function, and \code{"Fisher"} for Fisher information function.}
}

\value{
 The required maximum expected information for the selected item.
 }
 
\details{
 The MEI (van der Linden, 1998; van der Linden and Pashley, 2000) can be used as a rule for selecting the next item in the CAT process
 (see also Choi and Swartz, 2009). This command serves as a subroutine for the \code{\link{nextItem}} function.

 Let \emph{k} be the number of administered items, and set \eqn{x_1, ..., x_k} as the provisional response pattern. Set \eqn{\hat{\theta}_k} as the 
 provisional ability estimate (with the first \emph{k} responses) and let \emph{j} be the item of interest (not previously administered). Set also \eqn{P_j(\theta)}
 as the probability of answering item \emph{j} correctly for a given ability level \eqn{\theta}, and set \eqn{Q_j(\theta)=1-P_j(\theta)}. Finally, set
 \eqn{\hat{\theta}_{k+1}^0} and \eqn{\hat{\theta}_{k+1}^1} as the ability estimates computed under the condition that the response to item \emph{j} is 0 or 1 
 respectively (that is, if the response pattern is updated by 0 or 1 for item \emph{j}). Then, the MEI for item \emph{j} equals
 \deqn{MEI_j = P_j(\hat{\theta}_k)\,I_j(\hat{\theta}_{k+1}^1) +  Q_j(\hat{\theta}_k)\,I_j(\hat{\theta}_{k+1}^0)} 
 where \eqn{I_j(\theta)} is the information function for item \emph{j}.

 Two types of information functions are available. The first one is the observed information function, defined as 
 \deqn{I_j(\theta) = -\frac{\partial^2}{\partial\,\theta^2} \,\log P_j(\theta).}
 (van der Linden, 1998). The second one is Fisher information function:
 \deqn{I_j(\theta) = -E\,\left[\frac{\partial^2}{\partial\,\theta^2} \,\log P_j(\theta)\right].}
 Under the 1PL and the 2PL models, these functions are identical (Veerkamp, 1996).

 The observed and Fisher information functions are specified by the \code{infoType} argument, with respective values \code{"observed"} and \code{"Fisher"}. By
 default, the observed information function is considered (Choi and Swartz, 2009; van der Linden, 1998).

 The estimator of provisional ability is defined by means of the arguments \code{method}, \code{priorDist}, \code{priorPar}, \code{D}, \code{range} and 
 \code{parInt} of the \code{\link{thetaEst}} function. See the corresponding help file for further details.

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

 van der Linden, W. J. (1998). Bayesian item selection criteria for adaptive testing. \emph{Psychometrika, 63}, 201-216.

 van der Linden, W. J., and Pashley, P. J. (2000). Item selection and ability estimlation in adaptive testing. In W. J. van der Linden and C. A. W. Glas (Eds.), 
 \emph{Computerized adaptive testing. Theory and practice} (pp. 1-25). Boston, MA: Kluwer. 

 Veerkamp, W. J. J. (1996). \emph{Statistical inference for adaptive testing}. Internal report. Enschede, The Netherlands: 
 University of Twente. 
 }

\author{
    David Magis \cr
    Post-doctoral researcher, FNRS (Fonds National de la Recherche Scientifique) \cr
    Department of Mathematics, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{Ii}}, \code{\link{OIi}}, \code{\link{nextItem}}, \code{\link{integrate.xy}}, \code{\link{thetaEst}}
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
 MEI(bank, 1, x, 0, it)

 # With Fisher information instead
 MEI(bank, 1, x, 0, it, infoType="Fisher")

 # With WL estimator instead
 MEI(bank, 1, x, 0, it, method="WL")
 }
 }

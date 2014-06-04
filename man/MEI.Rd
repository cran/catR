\name{MEI}
\alias{MEI}

\title{(Maximum) Expected Information (MEI)}

\description{
 This command returns the expected information (EI) for a given item (under dichotomous and polytomous IRT models), as used for Maximum Expected Information (MEI) criterion.
 }

\usage{
MEI(itemBank, item, x, theta, it.given, model = NULL, method = "BM", 
 	priorDist = "norm", priorPar = c(0, 1), D = 1, range = c(-4, 4), 
 	parInt = c(-4, 4, 33), infoType = "observed")
 }

\arguments{
 \item{itemBank}{numeric: a suitable matrix of item parameters. See \bold{Details}.}
 \item{item}{numeric: the item (referred to as its rank in the item bank) for which the expected information must be computed.}
 \item{x}{numeric: a vector of item responses, coded as 0 or 1 only (for dichotomous items) or from 0 to the number of response categories minus one (for polytomous items).}
 \item{theta}{numeric: the provisional ability estimate.}
 \item{it.given}{numeric: a suitable matrix of item parameters for previously administered items. The number of rows of \code{it.given} must be equal to the length of \code{x}.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
 \item{method}{character: the ability estimator. Possible values are \code{"BM"} (default), \code{"ML"} and \code{"WL"}. See \bold{Details}.}
 \item{priorDist}{character: specifies the prior distribution. Possible values are \code{"norm"} (default), \code{"unif"} and \code{"Jeffreys"}.
 Ignored if \code{method} is neither \code{"BM"} nor \code{"EAP"}. See \bold{Details}.}
 \item{priorPar}{numeric: vector of two components specifying the prior parameters (default is \code{c(0,1)}) of the prior ability distribution. Ignored if \code{method} is neither \code{"BM"}
 nor \code{"EAP"}, or if \code{priorDist="Jeffreys"}. See \bold{Details}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952).} 
 \item{range}{numeric: vector of two components specifying the range wherein the ability estimate must be looked for (default is \code{c(-4,4)}). Ignored if \code{method=="EAP"}.}
 \item{parInt}{numeric: vector of three components, holding respectively the values of the arguments \code{lower}, \code{upper} and \code{nqp} of the
 \code{\link{eapEst}} command. Default vector is (-4, 4, 33). Ignored if \code{method} is not \code{"EAP"}.}
 \item{infoType}{character: the type of information function to be used. Possible values are \code{"observed"} (default) for observed information function, and \code{"Fisher"} for Fisher information function.}
}

\value{
 The required expected information for the selected item.
 }
 
\details{
 The MEI (van der Linden, 1998; van der Linden and Pashley, 2000) can be used as a rule for selecting the next item in the CAT process (see also Choi and Swartz, 2009), both with dichotomous and polytomous IRT models. This command serves as a subroutine for the \code{\link{nextItem}} function.

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{itemBank} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model, \code{"PCM"} for Partical Credit Model, \code{"GPCM"} for Generalized Partial Credit Model, \code{"RSM"} for Rating Scale Model and \code{"NRM"} for Nominal Response Model. The \code{itemBank} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.

Under polytomous IRT models, let \emph{k} be the number of administered items, and set \eqn{x_1, ..., x_k} as the provisional response pattern (where each response \eqn{x_l} takes values in \eqn{\{0, 1, ..., g_l\}}). Set \eqn{\hat{\theta}_k} as the provisional ability estimate (with the first \emph{k} responses) and let \emph{j} be the item of interest (not previously administered). Set also \eqn{P_{jt}(\theta)} as the probability of answering response category \emph{t} to item \emph{j} for a given ability level \eqn{\theta}. Finally, set
 \eqn{\hat{\theta}_{k+1}^t} as the ability estimates computed under the condition that the response to item \emph{j} is \emph{t} (with \eqn{t=0, ..., g_j}). Then, the EI for item \emph{j} equals
 \deqn{EI_j = \sum_{t=0}^{g_j} P_{jt}(\hat{\theta}_k)\,I_j(\hat{\theta}_{k+1}^t)} 
 where \eqn{I_j(\theta)} is the information function for item \emph{j}.

In case of dichotomous IRT models, all \eqn{g_l} values reduce to 1, so that item responses \eqn{x_l} equal either 0 or 1. Set simply \eqn{P_j(\theta)} as the probability of answering item \emph{j} correctly for a given ability level \eqn{\theta}, and  set \eqn{Q_j(\theta)=1-P_j(\theta)}. Finally, set \eqn{\hat{\theta}_{k+1}^0} and \eqn{\hat{\theta}_{k+1}^1} as the ability estimates computed under the condition that the response to item \emph{j} is 0 or 1 respectively (that is, if the response pattern is updated by 0 or 1 for item \emph{j}). Then, the EI for item \emph{j} reduces to
 \deqn{EI_j = P_j(\hat{\theta}_k)\,I_j(\hat{\theta}_{k+1}^1) +  Q_j(\hat{\theta}_k)\,I_j(\hat{\theta}_{k+1}^0)}.

 Two types of information functions are available. The first one is the observed information function, defined as 
 \deqn{OI_j(\theta) = -\frac{\partial^2}{\partial \theta^2} \,\log L(\theta | x_j)}
 (van der Linden, 1998), where \eqn{L(\theta | x_j)} is the likelihood related to item \emph{j}. The second one is Fisher information function:
 \deqn{I_j(\theta) = -E\,\left[\frac{\partial^2}{\partial\,\theta^2} \,\log L(\theta | x_j))\right].}
 Under the 1PL and the 2PL models, these functions are identical (Veerkamp, 1996). See also \code{\link{OIi}}. 

 The observed and Fisher information functions are specified by the \code{infoType} argument, with respective values \code{"observed"} and \code{"Fisher"}. By default, the observed information function is considered (Choi and Swartz, 2009; van der Linden, 1998).

 The estimator of provisional ability is defined by means of the arguments \code{method}, \code{priorDist}, \code{priorPar}, \code{D}, \code{range} and \code{parInt} of the \code{\link{thetaEst}} function. See the corresponding help file for further details.

The provisional response pattern and the related item parameters are provided by the arguments \code{x} and \code{it.given} respectively. The target item (for which the maximum information computed) is given by its number in the item bank, through the
 \code{item} argument. 

Note that the provisional response pattern \code{x} can also be set to \code{NULL} (which is useful when the number \code{nrItems} of starting items is set to zero). In this case, \code{it.given} must be a matrix with zero rows, such as e.g., \code{itemBank[NULL,]}. 
 }

\references{
 Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}.
 Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Choi, S. W., and Swartz, R. J. (2009). Comparison of CAT item selection criteria for polytomous items. \emph{Applied Psychological Measurement, 32}, 419-440. doi: 10.1177/0146621608327801

 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error}. Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. URL \url{http://www.jstatsoft.org/v48/i08/}

 van der Linden, W. J. (1998). Bayesian item selection criteria for adaptive testing. \emph{Psychometrika, 63}, 201-216. doi: 10.1007/BF02294775

 van der Linden, W. J., and Pashley, P. J. (2000). Item selection and ability estimation in adaptive testing. In W. J. van der Linden and C. A. W. Glas (Eds.), 
 \emph{Computerized adaptive testing. Theory and practice} (pp. 1-25). Boston, MA: Kluwer. 

 Veerkamp, W. J. J. (1996). \emph{Statistical inference for adaptive testing}. Internal report. Enschede, The Netherlands: 
 University of Twente. 
 }

\author{
    David Magis \cr
    Department of Education, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{Ii}}, \code{\link{OIi}}, \code{\link{nextItem}}, \code{\link{thetaEst}}, \code{\link{genPolyMatrix}} 
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
 it.given <- bank[c(15, 20),]

 # Creation of a response pattern
 x <- c(0, 1)

 # MEI for item 1, provisional ability level 0
 MEI(bank, 1, x, 0, it.given)

 # With Fisher information instead
 MEI(bank, 1, x, 0, it.given, infoType = "Fisher")

 # With WL estimator instead
 MEI(bank, 1, x, 0, it.given, method = "WL")


## Polytomous models ##

 # Generation of an item bank under GRM with 100 items and at most 4 categories
 m.GRM <- genPolyMatrix(100, 4, "GRM")
 m.GRM <- as.matrix(m.GRM)

 # Selection of two arbitrary items (15 and 20) 
it.given <- m.GRM[c(15, 20),]

 # Generation of a response pattern (true ability level 0)
 x <- genPattern(0, it.given, model = "GRM")

 # EPV for item 1, provisional ability level 0
 MEI(m.GRM, 1, x, 0, it.given, model = "GRM")

 # With WL method
 MEI(m.GRM, 1, x, 0, it.given, model = "GRM", method = "WL")

 # With Fisher information
 MEI(m.GRM, 1, x, 0, it.given, model = "GRM", infoType = "Fisher")

 # Generation of an item bank under PCM with 20 items and 4 categories
 m.PCM <- genPolyMatrix(20, 4, "PCM", same.nrCat = TRUE)
 m.PCM <- as.matrix(m.PCM)

 # Selection of two arbitrary items (15 and 20) 
 it.given <- m.PCM[c(15, 20),]

 # Generation of a response pattern (true ability level 0)
 x <- genPattern(0, it.given, model = "PCM")

 # EPV for item 1, provisional ability level 0
 MEI(m.PCM, 1, x, 0, it.given, model = "PCM")

 # With WL method
 MEI(m.PCM, 1, x, 0, it.given, model = "PCM", method = "WL")

 # With Fisher information
 MEI(m.PCM, 1, x, 0, it.given, model = "PCM", infoType = "Fisher") 
}
 }


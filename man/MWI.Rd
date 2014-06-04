\name{MWI}
\alias{MWI}

\title{Maximum likelihood weighted information (MLWI) and maximum posterior weighted information (MPWI)}

\description{
 This command returns the value of the likelihood (MLWI) or the posterior (MPWI) weighted information for a given item and an item bank (both under dichotomous and polytomous IRT models).
 }

\usage{
MWI(itemBank, item, x, it.given, model = NULL, lower = -4, upper = 4, nqp = 33, 
 	type = "MLWI", priorDist = "norm", priorPar = c(0, 1), D = 1)
 }

\arguments{
 \item{itemBank}{numeric: a suitable matrix of item parameters. See \bold{Details}.}
 \item{item}{numeric: the item (referred to as its rank in the item bank) for which the expected information must be computed.}
 \item{x}{numeric: a vector of item responses, coded as 0 or 1 only (for dichotomous items) or from 0 to the number of response categories minus one (for polytomous items). }
 \item{it.given}{numeric: a suitable matrix of item parameters for previously administered items. The number of rows of \code{it.given} must be equal to the length of \code{x}.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
 \item{lower}{numeric: the lower bound for numerical integration (default is -4).}
 \item{upper}{numeric: the upper bound for numerical integration (default is 4).}
 \item{nqp}{numeric: the number of quadrature points (default is 33).}
 \item{type}{character: the type of information to be computed. Possible values are \code{"MLWI"} (default) and \code{"MPWI"}. See \bold{Details}.}
 \item{priorDist}{character: the prior ability distribution. Possible values are \code{"norm"} (default) for the normal distribution, and \code{"unif"} for the uniform distribution. Ignored if \code{type} is not \code{"MPWI"}.}
 \item{priorPar}{numeric: a vector of two components with the prior parameters. If \code{priorDist} is \code{"norm"}, then \code{priorPar} contains the mean and the standard deviation of the normal distribution. If \code{priorDist} is \code{"unif"}, then \code{priorPar} contains the bounds of the uniform distribution. The default values are 0 and 1 respectively. Ignored if \code{type} is not \code{"MPWI"}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952). Ignored if \code{model} is not \code{NULL}.} 
}

\value{
 The required (likelihood or posterior) weighted information for the selected item.
 }
 
\details{
 Both the MLWI (Veerkamp and Berger, 1997) and the MPWI (van der Linden, 1998; van der Linden and Pashley, 2000) can be used as rules for selecting the next item in the CAT process (see also Choi and Swartz, 2009), both under dichotomous and polytomous IRT models. This command serves as a subroutine for the \code{\link{nextItem}} function.

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{itemBank} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model, \code{"PCM"} for Partical Credit Model, \code{"GPCM"} for Generalized Partial Credit Model, \code{"RSM"} for Rating Scale Model and \code{"NRM"} for Nominal Response Model. The \code{itemBank} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.

Under polytomous IRT models, let \emph{k} be the number of administered items, and set \eqn{x_1, ..., x_k} as the provisional response pattern (where each response \eqn{x_l} takes values in \eqn{\{0, 1, ..., g_l\}}). Set also
 \eqn{I_j(\theta)} as the information function of item \emph{j} evaluated at \eqn{\theta}, and set \eqn{L(\theta | x_1, ..., x_k)} as the likelihood function evaluated at \eqn{\theta}, given the provisional response pattern. Then, the LWI for item \emph{j} is given by
 \deqn{LWI_j = \int I_j(\theta) L(\theta | x_1, ..., x_k) d\theta} and the PWI by  
 \deqn{PWI_j = \int I_j(\theta) \pi(\theta) L(\theta | x_1, ..., x_k) d\theta} 
 where \eqn{\pi(\theta)} is the prior distribution of the ability level. 

In case of dichotomous IRT models, all \eqn{g_l} values reduce to 1, so that item responses \eqn{x_l} equal either 0 or 1. But except from this difference, the previous definitions of LWI and PWI remain valid.

 These integrals are approximated by the \code{\link{integrate.catR}} function. The range of integration is set up 
 by the arguments \code{lower}, \code{upper} and \code{nqp}, giving respectively the lower bound, the upper bound and the number of quadrature points. The default range goes from -4 to 4 with length 33 (that is, by steps of 0.25).

 The argument \code{type} defines the type of information to be computed. The default value, \code{"MLWI"}, computes the MLWI value, while the MPWI value is obtained with \code{type="MPWI"}. For the latter, the \code{priorDist} and \code{priorPar} arguments fix the prior ability distribution. 
 The normal distribution is set up by \code{priorDist="norm"} and then, \code{priorPar} contains the mean and the standard deviation of the normal distribution. If \code{priorDist} is \code{"unif"}, then the uniform distribution is considered, and \code{priorPar} fixes the lower and upper bounds of that 
 uniform distribution. By default, the standard normal prior distribution is assumed. This argument is ignored whenever \code{method} is not \code{"MPWI"}.

The provisional response pattern and the related item parameters are provided by the arguments \code{x} and \code{it.given} respectively. The target item (for which the maximum information is computed) is given by its rank number in the item bank, through the \code{item} argument. 

Note that the provisional response pattern \code{x} can also be set to \code{NULL} (which is useful when the number \code{nrItems} of starting items is set to zero). In this case, \code{it.given} must be a matrix with zero rows, such as e.g., \code{itemBank[NULL,]}. In this very specific configuration, the likelihood function \eqn{L(\theta | x_1, ..., x_k)} reduces to the constant value 1 on the whole \eqn{\theta} range (that is, a uniform likelihood). 
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

 Veerkamp, W. J. J., and Berger, M. P. F. (1997). Some new item selection criteria for adaptive testing. \emph{Journal of Educational and Behavioral Statistics, 22}, 203-226. doi: 10.3102/10769986022002203
 }

\author{
    David Magis \cr
    Department of Education, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{Ii}}, \code{\link{nextItem}}, \code{\link{integrate.catR}}, \code{\link{genPolyMatrix}}
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

 # MLWI for item 1
 MWI(bank, 1, x, it.given)

 # MPWI for item 1
 MWI(bank, 1, x, it.given, type = "MPWI")

 # MLWI for item 1, different integration range
 MWI(bank, 1, x, it.given, lower = -2, upper = 2, nqp = 20)

 # MPWI for item 1, uniform prior distribution on the range [-2,2]
 MWI(bank, 1, x, it.given, type = "MPWI", priorDist = "unif", priorPar = c(-2, 2))


## Polytomous models ##

 # Generation of an item bank under GRM with 100 items and at most 4 categories
 m.GRM <- genPolyMatrix(100, 4, "GRM")
 m.GRM <- as.matrix(m.GRM)

 # Selection of two arbitrary items (15 and 20) 
 it.given <- m.GRM[c(15, 20),]

 # Generation of a response pattern (true ability level 0)
 x <- genPattern(0, it.given, model = "GRM")

 # MLWI for item 1
 MWI(m.GRM, 1, x, it.given, model = "GRM")

 # MPWI for item 1
 MWI(m.GRM, 1, x, it.given, model = "GRM", type = "MPWI")

 # MLWI for item 1, different integration range
 MWI(m.GRM, 1, x, it.given, model = "GRM", lower = -2, upper = 2, nqp = 20)

 # MPWI for item 1, uniform prior distribution on the range [-2,2]
 MWI(m.GRM, 1, x, it.given, model = "GRM", type = "MPWI", priorDist = "unif", 
 priorPar = c(-2, 2))


 # Generation of an item bank under PCM with 20 items and 4 categories
 m.PCM <- genPolyMatrix(20, 4, "PCM", same.nrCat = TRUE)
 m.PCM <- as.matrix(m.PCM)

 # Selection of two arbitrary items (15 and 20) 
 it.given <- m.PCM[c(15, 20),]

 # Generation of a response pattern (true ability level 0)
 x <- genPattern(0, it.given, model = "PCM")

  # MLWI for item 1
 MWI(m.PCM, 1, x, it.given, model = "PCM")

 # MPWI for item 1
 MWI(m.PCM, 1, x, it.given, model = "PCM", type = "MPWI")

 # MLWI for item 1, different integration range
 MWI(m.PCM, 1, x, it.given, model = "PCM", lower = -2, upper = 2, nqp = 20)

 # MPWI for item 1, uniform prior distribution on the range [-2,2]
 MWI(m.PCM, 1, x, it.given, model = "PCM", type = "MPWI", priorDist = "unif", 
 priorPar = c(-2, 2))
 }
 }
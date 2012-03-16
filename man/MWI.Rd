\name{MWI}
\alias{MWI}

\title{Maximum likelihood weighted information (MLWI) and maximum posterior weighted information (MPWI)}

\description{
 This command returns the maximum likelihood (MLWI) or the maximum posterior (MPWI) weighted information for a given item and an item bank.
 }

\usage{
 MWI(itemBank, item, x, it, lower=-4, upper=4, nqp=33, 
 	type="MLWI", priorDist="norm", priorPar=c(0,1))
 }

\arguments{
 \item{itemBank}{an item bank, i.e. a list of class \code{itBank}, typically an output of the function \code{createItemBank}.}
 \item{item}{numeric: the item (referred to as its rank in the item bank) for which the maximum information must be computed.}
 \item{x}{binary: a vector of item responses, coded as 0 or 1 only.}
 \item{it}{numeric: a matrix with one row per item and four columns, with the values of the discrimination, the difficulty,
           the pseudo-guessing and the inattention parameters (in this order). The number of rows of \code{it} must be equal to the
          length of \code{x}.}
 \item{lower}{numeric: the lower bound for numercal integration (default is -4).}
 \item{upper}{numeric: the upper bound for numercal integration (default is 4).}
 \item{nqp}{numeric: the number of quadrature points (default is 33).}
 \item{type}{character: the type of information to be computed. Possible values are \code{"MLWI"} (default) and \code{"MPWI"}. See \bold{Details}.}
 \item{priorDist}{character: the prior ability distribution. Possible values are \code{"norm"} (default) for the normal distribution, and \code{"unif"} for
                  the uniform distribution. Ignored if \code{type} is not \code{"MPWI"}.}
 \item{priorPar}{numeric: a vector of two components with the prior parameters. If \code{priorDist} is \code{"norm"}, then \code{priorPar} contains the 
                 mean and the standard deviation of the normal distribution. If \code{priorDist} is \code{"unif"}, then \code{priorPar} contains the bounds
                 of the uniform distribution. The default values are 0 and 1 respectively. Ignored if \code{type} is not \code{"MPWI"}.}
}

\value{
 The required maximum information for the selected item.
 }
 
\details{
 Both the MLWI (Veerkamp and Berger, 1997) and the MPWI (van der Linden, 1998; van der Linden and Pashley, 2000) can be used
 as rules for selecting the next item in the CAT process (see also Choi and Swartz, 2009). This command serves as a subroutine for the
 \code{\link{nextItem}} function.

 Let \emph{k} be the number of administered items, and set \eqn{x_1, ..., x_k} as the binary responses to the first \emph{k} administered items. Set also
 \eqn{I_j(\theta)} as the information function of item \emph{j} evaluated at \eqn{\theta}, and set \eqn{L(\theta | x_1, ..., x_k)} as
 the likelihood function evaluated at \eqn{\theta}, given the provisional response pattern. Then, the MLWI for item \emph{j} is given by
 \deqn{MLWI_j = \int I_j(\theta) L(\theta | x_1, ..., x_k) d\theta} and the MPWI by  
 \deqn{MPWI_j = \int I_j(\theta) \pi(\theta) L(\theta | x_1, ..., x_k) d\theta} 
 where \eqn{\pi(\theta)} is the prior distribution of the ability level. 

 These integrals are approximated by the \code{\link{integrate.xy}} function from the package \code{sfsmisc}. The range of integration is set up 
 by the arguments \code{lower}, \code{upper} and \code{nqp}, giving respectively the lower bound, the upper bound and the number of quadrature points.
 The default range goes from -4 to 4 with length 33 (that is, by steps of 0.25).

 The argument \code{type} defines the type of information to be computed. The default value, \code{"MLWI"}, computes the MLWI value, while the MPWI
 value is obtained with \code{type="MPWI"}. For the latter, the \code{priorDist} and \code{priorPar} arguments fix the prior ability distribution. 
 The normal distribution is set up by \code{priorDist="norm"} and then, \code{priorPar} contains the mean and the standard deviation of the normal distribution. 
 If \code{priorDist} is \code{"unif"}, then the uniform distribution is considered, and \code{priorPar} fixes the lower and upper bounds of that 
 uniform distribution. By default, the standard normal prior distribution is assumed. This argument is ignored whenever \code{method} is not \code{"MPWI"}.

 The item bank is provided through the argument \code{itemBank}. The provisional response pattern and the related item parameters are provided by the arguments
 \code{x} and \code{it} respectively. The target item (for which the maximum information computed) is given by its number in the item bank, through the
 \code{item} argument. 
 }

\references{
 Choi, S. W., and Swartz, R. J. (2009). Comparison of CAT item selection criteria for polytomous items. \emph{Applied PScyhological Measurement, 32}, 419-440.

 Magis, D., and Raiche, G. (in press). Random generation of response patterns under computerized adaptive testing
 with the R package \emph{catR}. \emph{Journal of Statistical Software}.

 van der Linden, W. J. (1998). Bayesian item selection criteria for adaptive testing. \emph{Psychometrika, 63}, 201-216.

 van der Linden, W. J., and Pashley, P. J. (2000). Item selection and ability estimlation in adaptive testing. In W. J. van der Linden and C. A. W. Glas (Eds.), 
 \emph{Computerized adaptive testing. Theory and practice} (pp. 1-25). Boston, MA: Kluwer. 

 Veerkamp, W. J. J., and Berger, M. P. F. (1997). Some new item selection criteria for adaptive testing. \emph{Journal of Educational and Behavioral Statistics, 22}, 203-226.
 }

\author{
    David Magis \cr
    Post-doctoral researcher, FNRS (Fonds National de la Recherche Scientifique) \cr
    Department of Mathematics, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{Ii}}, \code{\link{nextItem}}, \code{\link{integrate.xy}}
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

 # MLWI for item 1
 MWI(bank, 1, x, it)

 # MPWI for item 1
 MWI(bank, 1, x, it, type="MPWI")

 # MLWI for item 1, different integration range
 MWI(bank, 1, x, it, lower=-2, upper=2, nqp=20)

 # MPWI for item 1, uniform prior distribution on the range [-2,2]
 MWI(bank, 1, x, it, type="MPWI", priorDist="unif", priorPar=c(-2,2))
 }
 }

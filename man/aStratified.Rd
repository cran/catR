\name{aStratified}
\alias{aStratified}

\title{Item membership assessment for a-stratified sampling}

\description{
 This command allocates each item from the bank into its specified stratum, defined by increasing discrimination parameters for a-stratified sampling.
 }

\usage{
aStratified(itemBank, K, model = NULL)
 }

\arguments{
 \item{itemBank}{a matrix or data frame of item parameters. Must refer to either nay dichotomous IRT model or to polytomous GRM, MGRM or GPCM models. See \bold{Details}.}
\item{K}{either a single integer value specifying the number of strata, or an integer vector of required numbers of items per strata.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for allowed polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"} and \code{"GPCM"}. See \bold{Details}.}

}

\value{
 A vector with as many components as the number of items in the bank, with values \code{Group1} to \code{GroupG} where \code{G} is the number of strata (i.e. either \code{K}) or the length of vector \code{K}).}

 
\details{
a-stratified sampling (Chang and Ying, 1999) is a content balancing technique that splits the item bank into several strata, the latter being set by increasing discrimination (a) parameters. The first stratum holds the least discriminative items, the second stratum the second set of least discriminative items etc., and the last stratum the most discriminative items. Content balancing item selection is then performed with equiprobable item selection across the strata.

The specification of the strata can be done in two ways, through the \code{K} argument. First, \code{K} can take any positive integer value (larger than one and smaller than the number of items): it specifies then the number of strata to create. Each stratum will hold approximately the same number of items, and in case the ratio between the number of items in the bank and \code{K} is not an integer value, then the last stratum will hold more items than the other ones. Second, \code{K} can be a vector of integer values that specify the number of items per strata. In this second case, the first stratum will hold as many items as the first component of \code{K}, the second startum as many items as the second component of \code{K}, and so on. Note that the sum of \code{K} shoudl coincide with the number of items in the bank.

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{it} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model
 (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model and \code{"GPCM"} for Generalized Partial Credit Model. \emph{Other polytomous IRT models are not allowed for this function and an error message will be returned if the input item bank does not follow one of the three aforementioned models}. The \code{it} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.

The output is a vector with character values \code{"Group1"} for the first startum, \code{"Group2"} for the second stratum etc. It has the same length as the number of items in the \code{itemBank} matrix and permits one-to-one identification of the items in each stratum.
 }

\references{
 Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}. Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

Chang, H.-H., and Ying, Z. (1999). A-stratified multistage computerized adaptive testing. \emph{Applied Psychological Measurement, 23}, 211-222. \doi{10.1177/01466219922031338}

Magis, D. and Barrada, J. R. (2017). Computerized Adaptive Testing with R: Recent Updates of the Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{Code Snippets}, \emph{76(1)}, 1-18. \doi{10.18637/jss.v076.c01}

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. \doi{10.18637/jss.v048.i08}
}

\author{
    David Magis \cr
    Department of Psychology, University of Liege, Belgium \cr
    \email{david.magis@uliege.be}
 }

\seealso{
\code{\link{genPolyMatrix}}
} 

\examples{

## Dichotomous models ##

 # Loading the 'tcals' parameters 
 data(tcals)

 # Selecting item parameters only
 tcals <- as.matrix(tcals[,1:4])

 # Creation of five strata with equal length (17 items each)
 aStratified(tcals, 5)

 # Creation of four strata with prespecified lengths 
 res <- aStratified(tcals, K = c(20, 25, 10, 30))
 table(res) # as expected


## Polytomous models ##

 # GRM with 100 items
 mat <- genPolyMatrix(100, 4, model = "GRM")

 # Creation of four strata with equal length (25 items each)
 aStratified(mat, 5, model = "GRM")
 }

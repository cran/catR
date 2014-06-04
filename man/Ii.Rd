\name{Ii}
\alias{Ii}

\title{Item information functions, first and second derivatives (dichotomous and polytomous models)}

\description{
 This command returns the Fisher information functions for a given ability value and a given matrix of item parameters under either the 4PL model or any suitable polytomous model. Numerical values of the first and second derivatives of the item information functions are also returned.
 }

\usage{
Ii(th, it, model = NULL, D = 1)
 }

\arguments{
 \item{th}{numeric: the ability value.}
 \item{it}{numeric: a suitable matrix of item parameters. See \bold{Details}.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952). Ignored if \code{model} is not \code{NULL}.} 
}

\value{
 A list with three arguments:
 \item{Ii}{the vector with item informations (one value per item).}
 \item{dIi}{the vector with first derivatives of the item information functions (one value per item).}
 \item{d2Ii}{the vector with second derivatives of the item information functions (one value per item).}
 }
 
\details{
 The first and second derivatives are computed algebraically, either from the four-parameter logistic (4PL) model
 (Barton and Lord, 1981) or from the corresponding polytomous model. These derivatives are necessary for both the estimation of ability and the computation of related standard errors. 

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{it} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model
 (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model, \code{"PCM"} for Partical Credit Model, \code{"GPCM"} for Generalized Partial Credit Model, \code{"RSM"} for Rating Scale Model and \code{"NRM"} for Nominal Response Model. The \code{it} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.
 }

\references{
 Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}. Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error}. Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University.

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R 
 Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. URL \url{http://www.jstatsoft.org/v48/i08/}
 }

\author{
    David Magis \cr
    Department of Education, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{Pi}}, \code{\link{thetaEst}}, \code{\link{genPolyMatrix}} 
} 

\examples{
\dontrun{

## Dichotomous models ##
 
 # Loading the 'tcals' parameters 
 data(tcals)

 # Selecting item parameters only
 tcals <- as.matrix(tcals[,1:4])
 
 # Item information functions and derivatives
 # (various th and D values)
 Ii(th = 0, tcals)
 Ii(th = 0, tcals, D = 1.702)
 Ii(th = 1, tcals)

## Polytomous models ##

 # Generation of an item bank under GRM with 100 items and at most 4 categories
 m.GRM <- genPolyMatrix(100, 4, "GRM")
 m.GRM <- as.matrix(m.GRM)

 # Computation of item information and derivatives for ability level 0
 Ii(0, m.GRM, model = "GRM")

 # Generation of a item bank under PCM with 20 items and at most 3 categories
 m.PCM <- genPolyMatrix(20, 3, "PCM")
 m.PCM <- as.matrix(m.PCM)

 # Computation of item information and derivatives for ability level 1
 Ii(1, m.PCM, model = "PCM")
 }
 }

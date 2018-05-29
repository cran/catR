\name{breakBank}
\alias{breakBank}

\title{Breaking the item bank in item parameters and group membership (for content balancing)}

\description{
 This command "breaks" the item bank in two parts, the item parameters as a numeric matrix, and the group membership of the items as a vector of factor levels. These two elements can be used separately for content balancing purposes, among others.
 }

\usage{
breakBank(itemBank)
 }

\arguments{
 \item{itemBank}{a matrix or data frame with item parameters in the first columns and group membership in the last column.}
}

\value{
 A list with two arguments:
 \item{itemPar}{a numeric matrix with the item parameters.}
 \item{cbGroup}{a vector with factor names of item membership, with one name per item in \code{itemPar} element of the output list.}
 }
 
\details{
 The \code{breakBank} function is useful to split the original item bank in two parts, one holding the item parameters only, and one containing the names of item group membership. The former can then directly be plugged in adequate functions such as \code{\link{thetaEst}}.

The function works with both dichotomous and polytomous IRT item banks. In both cases, the group membership must be located as the last column of the matrix or data frame.  

Note that there is no check that the input is correct (that is, the group membership is located in the last column), as \code{breakBank} is mostly devoted to be used at the early stages of the \code{\link{randomCAT}} function.
 }

\references{
Magis, D. and Barrada, J. R. (2017). Computerized Adaptive Testing with R: Recent Updates of the Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{Code Snippets}, \emph{76(1)}, 1-18. \doi{10.18637/jss.v076.c01}

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. \doi{10.18637/jss.v048.i08}
}

\author{
    David Magis \cr
    Department of Psychology, University of Liege, Belgium \cr
    \email{david.magis@uliege.be}
 }

\seealso{
\code{\link{randomCAT}}
} 

\examples{

## Dichotomous models ##

 # Loading the 'tcals' parameters 
 data(tcals)

 # Breaking 'tcals'
 breakBank(tcals)


## Polytomous models ##

 # Creation of the 'cbList' list with arbitrary proportions
 cbList <- list(names =c ("Audio1", "Audio2", "Written1", "Written2",
        "Written3"), props = c(0.1, 0.2, 0.2, 0.2, 0.3))

 # NRM with 100 items
 mat <- genPolyMatrix(100, 4, model = "NRM", cbControl = cbList)

 # Breaking 'mat'
 breakBank(mat)
 }

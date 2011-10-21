\name{Ii}
\alias{Ii}

\title{Item information functions, first and second derivatives (4PL)}

\description{
 This command returns the Fisher information functions for a given matrix of item parameters of the 4PL model 
 and a given ability value. Numerical values of the first and second derivatives of the item information functions
 are also returned.
 }

\usage{
 Ii(th, it, D=1)
 }

\arguments{
 \item{th}{numeric: the ability value.}
 \item{it}{numeric: a matrix with one row per item and four columns, with the values of the discrimination, the difficulty,
           the pseudo-guessing and the inattention parameters (in this order).}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952).} 
}

\value{
 A list with three arguments:
 \item{Ii}{the vector with item informations (one value per item)}
 \item{dIi}{the vector with first derivatives of the item information functions (one value per item)}
 \item{d2Ii}{the vector with second derivatives of the item information functions (one value per item)}
 }
 
\details{
 The first and second derivatives are computed algebraically from the four-parameter logistic (4PL) model
 (Barton and Lord, 1981). These derivatives are necessary for both the estimation of ability and the computation of
 related standard errors. 
 }

\references{
 Barton, M.A., and Lord, F.M. (1981). An upper asymptote for the three-parameter logistic item-response model.
 Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Haley, D.C. (1952). Estimation of the dosage mortality relationship when the dose is subject to error. 
 Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 
 }

\author{
    David Magis \cr
    Post-doctoral researcher, FNRS (Fonds National de la Recherche Scientifique) \cr
    Department of Mathematics, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{Pi}}, \code{\link{thetaEst}}
} 

\examples{
\dontrun{

 # Loading the 'tcals' parameters 
 data(tcals)

 # Selecting item parameters only
 tcals <- as.matrix(tcals[,1:4])
 
 # Item information functions and derivatives
 # (various th and D values)
 Ii(th=0, tcals)
 Ii(th=0, tcals, D=1.702)
 Ii(th=1, tcals)
 }
 }

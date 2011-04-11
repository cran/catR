\name{Pi}
\alias{Pi}

\title{Item response probabilities, first, second and third derivatives (4PL)}

\description{
 This command returns the item reponse probabilities for a given matrix of item parameters of the 4PL model
 and a given ability value. Numerical values of the first, second and third derivatives of the response probabilities
 are also returned.
 }

\usage{
 Pi(th, it, D=1)
 }

\arguments{
 \item{th}{numeric: the ability value.}
 \item{it}{numeric: a matrix with one row per item and four columns, with the values of the discrimination, the difficulty,
           the pseudo-guessing and the inattention parameters (in this order).}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952).} 
}

\value{
 A list with four arguments:
 \item{Pi}{the vector with response probabilities (one value per item)}
 \item{dPi}{the vector with first derivatives of the response probabilities (one value per item)}
 \item{d2Pi}{the vector with second derivatives of the response probabilities (one value per item)}
 \item{d3Pi}{the vector with third derivatives of the response probabilities (one value per item)}
 }
 
\details{
 The first, second and third derivatives are computed algebraically from the four-parameter logistic (4PL) model
 (Barton and Lord, 1981). These derivatives are necessary for both the estimation of ability and the computation of
 related standard errors. 
 }

\note{
 Response probabilites exactly equal to zero are returned as \code{1e-10} values, as well as probabilities exactly equal to one which are
 returned as \code{1-1e-10} values. This is to permit the computation of ability estimates (with the \code{\link{thetaEst}} function) in such extreme cases.
 Many thanks to Pan Tong (University of Texas MD Anderson Cancer Center, USA) who noticed this problem.
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
\code{\link{Ii}}, \code{\link{thetaEst}}
} 

\examples{
\dontrun{

 # Loading the 'tcals' parameters 
 data(tcals)
 tcals <- as.matrix(tcals)
 
 # Response probabilities and derivatives (various th and D values)
 Pi(th=0, tcals)
 Pi(th=0, tcals, D=1.702)
 Pi(th=1, tcals)
 }
 }

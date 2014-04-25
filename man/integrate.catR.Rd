\name{integrate.catR}
\alias{integrate.catR}

\title{Numerical integration by linear interpolation (for catR internal use)}

\description{
 This command computes the integral of function \emph{f(x)} by providing values of \emph{x} and \emph{f(x)}, similarly to the \code{integrate.xy} function of the R package \code{sfsmisc}.
 }

\usage{
integrate.catR(x, y)
 }

\arguments{
 \item{x}{numeric: a vector of \emph{x} values for numerical integration.}
 \item{y}{numeric: a vector of numerical values corresponding to \emph{f(x)} values.}
}

\value{
 The approximated integral.
 }
 
\details{
 This function was written to compute "cheap" numerical integration by providing sequences of \emph{x} values and corresponding computed values \emph{f(x)}. It works similarly as the \code{integrate.xy} function when \code{use.spline=FALSE} is required. It was developed internally to eventually remove dependency of \code{catR} package to package \code{sfsmisc}.
 }


\references{
Maechler, M. et al. (2012). \emph{sfsmisc: Utilities from Seminar fuer Statistik ETH Zurich}. R package version 1.0-23.
  http://CRAN.R-project.org/package=sfsmisc
}

\author{
    David Magis \cr
    Department of Education, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{KL}} and the \code{integrate.xy} function in package \code{sfsmisc}
} 

\examples{
\dontrun{

 # Loading the 'tcals' parameters 
 x <- seq(from = -4, to = 4, length = 33)
 y <- exp(x)
 integrate.catR(x, y) # 54.86381

 # Comparison with integrate.xy
 require(sfsmisc)
 integrate.xy(x, y, use.spline = FALSE) # 54.86381
 integrate.xy(x, y) # 54.58058
 }
 }
 

\name{OIi}
\alias{OIi}

\title{Observed information function (4PL)}

\description{
 This command returns the observed information functions for a given matrix of item parameters of the 4PL model 
 and a given ability value.
 }

\usage{
 OIi(th, it, x, D=1)
 }

\arguments{
 \item{th}{numeric: the ability value.}
 \item{it}{numeric: a matrix with one row per item and four columns, with the values of the discrimination, the difficulty,
           the pseudo-guessing and the inattention parameters (in this order).}
 \item{x}{numeric: the item response (coded as 0 or 1). Can be either a single value or a vector of the same length of the number of items.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952).} 
}

\value{
 A vector with the observed item informations (one per item).
 }
 
\details{
 The observed information function for item \emph{j} is given by
 \deqn{-\frac{\partial^2}{\partial \theta^2} \,\log L(\theta | x_j)}
 where \eqn{\theta} is the ability level, \eqn{L} is the likelihood function and \eqn{x_j} is the item response. For dichotomous item response models with success
 probability \eqn{P_j(\theta)}, it takes the following form:
 \deqn{-\frac{\partial^2}{\partial \theta^2} \,\log L(\theta | x_j) = \frac{P_j \,Q_j \,{P_j'}^2- (x_j-P_j) \,[P_j \,Q_j \,P_j'' + {P_j}^2 \,(P_j-Q_j]}{{P_j}^2 \,{Q_j}^2}}
 where \eqn{P_j=P_j(\theta)}, \eqn{Q_j=1-P_j} and \eqn{P_j'} and \eqn{P_j''} are the first and second derivatives of \eqn{P_j} respectively. 

 Under the 2PL model, the observed information function is exactly equal to Fisher's information function
 \deqn{-E \left[\frac{\partial^2}{\partial \theta^2} \,\log L(\theta | x_j)\right] = \frac{{P_j'}^2}{P_j Q_j}}
 (van der Linden, 1998; Veerkamp, 1996).

 The observed information function is used to compute some item selection criteria, such as the Maximum Expected Information (MEI). See \code{\link{nextItem}}
 for further details.
 }


\references{
 Barton, M.A., and Lord, F.M. (1981). An upper asymptote for the three-parameter logistic item-response model.
 Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Haley, D.C. (1952). Estimation of the dosage mortality relationship when the dose is subject to error. 
 Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University.

 Magis, D., and Raiche, G. (in press). Random generation of response patterns under computerized adaptive testing
 with the R package \emph{catR}. \emph{Journal of Statistical Software}.

 van der Linden, W. (1998). Bayesian item selection criteria for adaptive testing. \emph{Psychometrika, 63}, 201-216.

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
\code{\link{createItemBank}}, \code{\link{nextItem}}
} 

\examples{
\dontrun{

 # Loading the 'tcals' parameters 
 data(tcals)
 
 # Selecting item parameters only
 tcals <- as.matrix(tcals[,1:4])
 
 # Observed information functions 
 # (various th, x and D values)
 OIi(th=0, tcals, x=0)
 OIi(th=0, tcals, x=0, D=1.702)
 OIi(th=0, tcals, x=1)
 OIi(th=1, tcals, x=1)
 }
 }

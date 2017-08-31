\name{genPattern}
\alias{genPattern}

\title{Random generation of item response patterns under dichotomous and polytomous IRT models}

\description{
 This command generates item responses patterns for a given matrix of item parameters of any specified dichotomous or polytomous IRT model and a given (set of) ability value(s). 
 }

\usage{
genPattern(th, it, model = NULL, D = 1, seed = NULL)
 }

\arguments{
 \item{th}{numeric: a vector of true underlynig ability values (can be a singlre value too).}
 \item{it}{numeric: a suitable matrix of item parameters. See \bold{Details}.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952). Ignored if \code{model} is not \code{NULL}.}
\item{seed}{either the random seed value or \code{NULL} (default). See \bold{Details}.} 
}

\value{
 If \code{th} holds a single value, output is a vector with the item responses in the order of appearance of the items in the \code{it} matrix. If \code{th} is a vector of numeric values, output is a response matrix with one row per \code{th} value and one column per item.
 }
 
\details{
 This function permits to randomly generate item responses for a set of given ability levels, specified by the argument \code{thr}, and for a given matrix of item parameters, specified by the argument \code{it}. Both dichotomous and polytomous IRT models can be considered and item responses are generated accordingly.

For dichotomous models, item responses are generated from Bernoulli draws, using the \code{\link{rbinom}} function. For polytomous models they are generated from darws from a multinomial distribution, using the \code{\link{rmultinom}} function. In both cases, success probabilities are obtained from the \code{\link{Pi}} function.

Note that for polytomous models, item responses are coded as 0 (for the first response category), 1 (for the second category), ..., until \eqn{g_j} (for the last category), in agreement with the notations used in the help files of, e.g., the \code{\link{genPolyMatrix}} function.

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{it} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model
 (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model, \code{"PCM"} for Partical Credit Model, \code{"GPCM"} for Generalized Partial Credit Model, \code{"RSM"} for Rating Scale Model and \code{"NRM"} for Nominal Response Model. The \code{it} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.

The random pattern generation can be fixed by setting \code{seed} to some numeric value. By default, \code{seed} is \code{NULL} and the random seed is not fixed.
 }


\references{
 Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}. Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error}. Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University.

Magis, D. and Barrada, J. R. (2017). Computerized Adaptive Testing with R: Recent Updates of the Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{Code Snippets}, \emph{76(1)}, 1-18. \doi{10.18637/jss.v076.c01}

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. \doi{10.18637/jss.v048.i08}
 }

\author{
    David Magis \cr
    Department of Education, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{rbinom}} and \code{\link{rmultinom}} for random draws; \code{\link{genPolyMatrix}}, \code{\link{Pi}}
} 

\examples{

## Dichotomous models ##
 
 # Loading the 'tcals' parameters 
 data(tcals)

 # Selecting item parameters only
 tcals <- as.matrix(tcals[,1:4])
 
 # Generation of a response pattern for ability level 0
 genPattern(th = 0, tcals)

 # Generation of 10 response patterns for ability levels randomly drawn from N(0,1)
 genPattern(th = rnorm(10), tcals)

 # Generation of a single response for the first item only
 genPattern(th = 0, tcals[1,])

## Polytomous models ##

 # Generation of an item bank under GRM with 100 items and at most 4 categories
 m.GRM <- genPolyMatrix(100, 4, "GRM")
 m.GRM <- as.matrix(m.GRM)

 # Generation of a response pattern for ability level 0
 genPattern(0, m.GRM, model = "GRM")

 # Generation of 10 response patterns for ability levels randomly drawn from N(0,1)
 genPattern(rnorm(10), m.GRM, model = "GRM")

 # Generation of a single response for the first item only
 genPattern(0, m.GRM[1,], model = "GRM")

 # Loading the cat_pav data
 data(cat_pav)
 cat_pav <- as.matrix(cat_pav)

 # Generation of a response pattern for ability level 0
 genPattern(0, cat_pav, model = "GPCM")

 # Generation of a single response for the first item only
 genPattern(0, cat_pav[1,], model = "GPCM")
 }

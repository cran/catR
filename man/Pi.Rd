\name{Pi}
\alias{Pi}

\title{Item response probabilities, first, second and third derivatives (dichotomous and polytomous models)}

\description{
 This command returns the item response probabilities for a given ability value and a given matrix of item parameters under either the 4PL model or any suitable polytomous model. Numerical values of the first, second and third derivatives of the response probabilities are also returned.
 }

\usage{
Pi(th, it, model = NULL, D = 1)
 }

\arguments{
 \item{th}{numeric: the ability value.}
 \item{it}{numeric: a suitable matrix of item parameters. See \bold{Details}.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952). Ignored if \code{model} is not \code{NULL}.} 
}

\value{
 Under dichotomous IRT models, a list with four arguments:
 \item{Pi}{the vector with response probabilities (one value per item).}
 \item{dPi}{the vector with first derivatives of the response probabilities (one value per item).}
 \item{d2Pi}{the vector with second derivatives of the response probabilities (one value per item).}
 \item{d3Pi}{the vector with third derivatives of the response probabilities (one value per item).}
Under polytomous IRT models, the aforementioned vectors are replaced by matrices with one row per item (labeled as \code{Item1}, \code{Item2} etc.) and one row per response category. 
 }
 
\details{
Whatever the IRT model, the response probabilities and first, second, and third derivatives are computed algebraically. These derivatives are necessary for both the estimation of ability and the computation of related standard errors. 

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{it} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model
 (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model (Samejima, 1969), \code{"MGRM"} for Modified Graded Response Model (Muraki, 1990), \code{"PCM"} for Partical Credit Model (Masters, 1982), \code{"GPCM"} for Generalized Partial Credit Model (Muraki, 1992), \code{"RSM"} for Rating Scale Model (Andrich, 1978) and \code{"NRM"} for Nominal Response Model (Bock, 1972). The \code{it} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.

The output list contains the response probabilities and the respective derivatives. In case of dichotomous models, only a vector of such values is returned, with one value per item. In case of polytomous models, matrices are returned instead, with one row per item and one column per response category. In case of unequal numbers of response categories (which may happen under GRM, PCM, GPCM and NRM), values for empty response categories are returned as \code{NA} values.
 }

\note{
 For dichotomous IRT models, response probabilites exactly equal to zero are returned as \code{1e-10} values, as well as probabilities exactly equal to one which are returned as \code{1-1e-10} values. This is to permit the computation of ability estimates (with the \code{\link{thetaEst}} function) in such extreme cases.

 Many thanks to Pan Tong (University of Texas MD Anderson Cancer Center, USA) who noticed this problem.
}


\references{
Andrich, D. (1978). A rating formulation for ordered response categories. \emph{Psychometrika, 43}, 561-573. \doi{10.1007/BF02293814}

Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}.
 Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

Bock, R. D. (1972). Estimating item parameters and latent ability when responses are scored in two or more nominal categories. \emph{Psychometrika, 37}, 29-51. \doi{10.1007/BF02291411}

 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error}. Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 

Magis, D. and Barrada, J. R. (2017). Computerized Adaptive Testing with R: Recent Updates of the Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{Code Snippets}, \emph{76(1)}, 1-18. \doi{10.18637/jss.v076.c01}

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. \doi{10.18637/jss.v048.i08}

Masters, G. N. (1982). A Rasch model for partial credit scoring. \emph{Psychometrika, 47}, 149-174. \doi{10.1007/BF02296272}

Muraki, E. (1990). Fitting a polytomous item response model to Likert-type data. \emph{Applied Psychological Measurement, 14}, 59-71. \doi{10.1177/014662169001400106}

Muraki, E. (1992). A generalized partial credit model: Application of an EM algorithm. \emph{Applied Psychological Measurement, 16}, 19-176. \doi{10.1177/014662169201600206}

Samejima, F. (1969). \emph{Estimation of latent ability using a response pattern of graded scores}. Psychometrika Monograph (vol. 17).
 }

\author{
    David Magis \cr
    Department of Education, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{Ii}}, \code{\link{thetaEst}}
} 

\examples{

## Dichotomous models ##

 # Loading the 'tcals' parameters 
 data(tcals)
 
 # Selecting item parameters only
 tcals <- as.matrix(tcals[,1:4])
 
 # Response probabilities and derivatives (various th and D values)
 Pi(th = 0, tcals)
 Pi(th = 0, tcals, D = 1.702)
 Pi(th = 1, tcals)

## Polytomous models ##

 # Generation of an item bank under GRM with 100 items and at most 4 categories
 m.GRM <- genPolyMatrix(100, 4, "GRM")
 m.GRM <- as.matrix(m.GRM)

 # Computation of probabilities and derivatives for ability level 0
 Pi(0, m.GRM, model = "GRM")

 # Loading the cat_pav data
 data(cat_pav)
 cat_pav <- as.matrix(cat_pav)

 # Computation of probabilities and derivatives for ability level 1
 Pi(1, cat_pav, model = "GPCM")
 }

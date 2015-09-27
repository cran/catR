\name{OIi}
\alias{OIi}

\title{Observed information function (dichotomous and polytomous models)}

\description{
 This command returns the observed information functions for a given ability value and a given matrix of item parameters, either under the 4PL model or any suitable polytomous IRT model. 
 }

\usage{
OIi(th, it, x, model = NULL, D = 1)
 }

\arguments{
 \item{th}{numeric: the ability value.}
 \item{it}{numeric: a suitable matrix of item parameters. See \bold{Details}.}
 \item{x}{numeric: a vector of item responses.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952). Ignored if \code{model} is not \code{NULL}.} 
}

\value{
 A vector with the observed item informations (one per item).
 }
 
\details{
 The observed information function for item \emph{j} is given by
 \deqn{OI_j= -\frac{\partial^2}{\partial \theta^2} \,\log L(\theta | x_j)}
 where \eqn{\theta} is the ability level, \eqn{L} is the likelihood function and \eqn{x_j} is the item response. 

For dichotomous IRT models with success probability \eqn{P_j(\theta)}, it takes the following form:
 \deqn{OI_j = \frac{P_j \,Q_j \,{P_j'}^2- (x_j-P_j) \,[P_j \,Q_j \,P_j'' + {P_j}^2 \,(P_j-Q_j]}{{P_j}^2 \,{Q_j}^2}}
 where \eqn{P_j=P_j(\theta)}, \eqn{Q_j=1-P_j} and \eqn{P_j'} and \eqn{P_j''} are the first and second derivatives of \eqn{P_j} respectively. 

For polytomous IRT models, set $X_j$ as the item response, taking values \eqn{k \in \{0, 1, ..., g_j\}}. Set \eqn{P_{jk}(\theta)=Pr(X_j=k | \theta)} as the probability of answering response category \eqn{k} and set \eqn{\tau_{jk}} as the boolean factor equal to 1 if \eqn{X_j=k} and 0 otherwise. Then, the observed information function for item \eqn{j} takes the following form:
\deqn{OI_j = \sum_{k=0}^{g_j} \tau_{jk}\,\left( \frac{{P_{jk}'(\theta)}^2}{{P_{jk}(\theta)}^2}-\frac{P_{jk}''(\theta)}{P_{jk}(\theta)}\right)}
with the same notations for the first and second derivatives as above.

 Under the 2PL model, the observed information function is exactly equal to Fisher's information function
 \deqn{I_j=-E \left[\frac{\partial^2}{\partial \theta^2} \,\log L(\theta | x_j)\right] = \frac{{P_j'}^2}{P_j Q_j}}
 (van der Linden, 1998; Veerkamp, 1996).

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{it} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model
 (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model, \code{"PCM"} for Partical Credit Model, \code{"GPCM"} for Generalized Partial Credit Model, \code{"RSM"} for Rating Scale Model and \code{"NRM"} for Nominal Response Model. The \code{it} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.

 The observed information function is used to compute some item selection criteria, such as the Maximum Expected Information (MEI). See \code{\link{nextItem}} and \code{\link{MEI}} for further details.
 }


\references{
 Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}.
 Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error}. Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University.

 Magis, D., and Raiche, G. (in press). Random generation of response patterns under computerized adaptive testing
 with the R package \emph{catR}. \emph{Journal of Statistical Software}.

 van der Linden, W. (1998). Bayesian item selection criteria for adaptive testing. \emph{Psychometrika, 63}, 201-216. doi: 10.1007/BF02294775

 Veerkamp, W. J. J. (1996). \emph{Statistical inference for adaptive testing}. Internal report. Enschede, The Netherlands: 
 University of Twente. 
 }

\author{
    David Magis \cr
    Department of Education, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{nextItem}}, \code{\link{MEI}}
} 

\examples{

## Dichotomous models ##

 # Loading the 'tcals' parameters 
 data(tcals)
 
 # Selecting item parameters only
 tcals <- as.matrix(tcals[,1:4])
 
 # Observed information functions 
 # (various th, x and D values)
 OIi(th = 0, tcals, x = 0)
 OIi(th = 0, tcals, x = 0, D = 1.702)
 OIi(th = 0, tcals, x = 1)
 OIi(th = 1, tcals, x = 1)


## Polytomous models ##

 # Generation of an item bank under GRM with 100 items and at most 4 categories
 m.GRM <- genPolyMatrix(100, 4, "GRM")
 m.GRM <- as.matrix(m.GRM)

 # Creation of a response pattern (true ability level 0)
 set.seed(1)
 x <- genPattern(0, m.GRM, model = "GRM")

 # Observed information functions (various th values)
 OIi(th = 0, m.GRM, x, model = "GRM")
 OIi(th = 1, m.GRM, x, model = "GRM")

 # Generation of an item bank under PCM with 20 items and 4 categories
 m.PCM <- genPolyMatrix(20, 4, "PCM", same.nrCat = TRUE)
 m.PCM <- as.matrix(m.PCM)

 # Creation of a response pattern (true ability level 0)
 set.seed(1)
 x <- genPattern(0, m.PCM, model = "PCM")

 # Observed information functions (various th values)
 OIi(th = 0, m.PCM, x, model = "PCM")
 OIi(th = 1, m.PCM, x, model = "PCM")
 }

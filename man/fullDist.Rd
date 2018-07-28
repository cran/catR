\name{fullDist}
\alias{fullDist}

\title{Full distribution of ability estimator (dichotomous models only)}

\description{
 This command returns the full distribution of the selected ability estimator, that is, the set of all possible ability estimates and related probabilities, for a given matrix of item parameters and ability estimate (or true value).
 }

\usage{
fullDist(th, it, method = "BM", priorDist = "norm", priorPar = c(0,1), 
 	weight = "Huber", tuCo  = 1, range = c(-4 ,4), parInt = c(-4, 4, 33))
 }

\arguments{
 \item{th}{numeric: the ability estimate of interest (can be a vector of estimates too).}
 \item{it}{numeric: a suitable matrix of item parameters. See \bold{Details}.}
 \item{method}{character: the ability estimator. Possible values are \code{"BM"} (default), \code{"ML"}, \code{"WL"}, \code{"EAP"} and \code{"ROB"}. See \bold{Details}.}
 \item{priorDist}{character: specifies the prior distribution. Possible values are \code{"norm"} (default), \code{"unif"} and \code{"Jeffreys"}. Ignored if \code{method} is neither \code{"BM"} nor \code{"EAP"}. See \bold{Details}.}
 \item{priorPar}{numeric: vector of two components specifying the prior parameters (default is \code{c(0,1)}) of the prior ability distribution. Ignored if \code{method} is neither \code{"BM"} nor \code{"EAP"}, or if \code{priorDist="Jeffreys"}. See \bold{Details}.}
\item{weight}{character: the type of weight function for the robust estimator. Possible values are \code{"Huber"} (default) and \code{"Tukey"}. Ignored if \code{method} is not \code{"ROB"} or if \code{model} is not \code{NULL}. See \bold{Details}.}
\item{tuCo}{numeric: the value of the tuning constant for the weight function (default is 1, suitable with \code{"Huber"} weight). Ignored if \code{method} is not \code{"ROB"} or if \code{model} is not \code{NULL}. See \bold{Details}.}
 \item{range}{numeric: vector of two components specifying the range wherein the ability estimate must be looked for (default is \code{c(-4,4)}). Ignored if \code{method=="EAP"}.}
 \item{parInt}{numeric: vector of three components, holding respectively the values of the arguments \code{lower}, \code{upper} and \code{nqp} of the \code{\link{eapEst}} command. Default vector is (-4, 4, 33). Ignored if \code{method} is not \code{"EAP"}.}
}

\value{
 A matrix with \eqn{(t+1)} columns (where \eqn{t} is the length of vector \code{th}), the first column holding all observable ability estimates and columns 2 to \eqn{(t+1)} with related probabilities of the corresponding components of \code{th}.
 }
 
\details{
The computation of the full distribution of an ability estimator is required to determine the exact standard error of ability, compared to the asympotic SE (ASE) displayed by former versions of the \code{\link{semTheta}} function. This distribution is computed as follows. 

First, all possible response patterns (for the given test length fixed by the number of items specified by \code{it} argument) are generated through the internal \code{dataGen} function. Second, ability estimation is performed with each generated pattern, and the corresponding pattern probability is computed using the set of item parameters \code{it} and the predefined ability level \code{th} of interest. These to components are eventually returned as the full distribution. 

In case of the Rasch (1PL) model, this long process can be shortened as the total test score is a proxy for ability estimation. Hence, with \eqn{n} items, only \eqn{(n+1)} patterns must be created, one for each test score from zero to \eqn{n}. Related probabilities can be derived using the Lord-Wingersky algorithm (1984) that is implemented internally through the \code{LW()} function.

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{it} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model
 (Barton and Lord, 1981). 

Five ability estimators are available: the maximum likelihood (ML) estimator (Lord, 1980), the Bayes modal (BM) estimator (Birnbaum, 1969), the expected a posteriori (EAP) estimator (Bock and Mislevy, 1982), the weighted likelihood (WL) estimator (Warm, 1989) and the robust estimator (Schuster & Yuan, 2011). The selected estimator is specified by the \code{method} argument, with values \code{"ML"}, \code{"BM"}, \code{"EAP"}, \code{"WL"} and \code{"ROB"} respectively.

 For the BM and EAP estimators, three prior distributions are available: the normal distribution, the uniform distribution and the Jeffreys' prior distribution (Jeffreys, 1939, 1946). The prior distribution is specified by the argument \code{priorPar}, with values \code{"norm"}, \code{"unif"} and \code{"Jeffreys"}, respectively. The \code{priorPar} argument is ignored if \code{method="ML"} or \code{method="WL"}.
 
 The argument \code{priorPar} determines either: the prior mean and standard deviation of the normal prior distribution (if \code{priorDist="norm"}), or the range for defining the prior uniform distribution (if \code{priorDist="unif"}). This argument is ignored if \code{priorDist="Jeffreys"}.
 
 The \code{eapPar} argument sets the range and the number of quadrature points for numerical integration in the EAP process. By default, it takes the vector value (-4, 4, 33), that is, 33 quadrature points on the range [-4; 4] (or, by steps of 0.25). See \code{\link{eapEst}} for further details. 

Robust estimation requires an appropriate weight function that depends on an accurate tuning constant. Suggested functions are the Huber weight (Schuester and Yuan, 2011) and the Tukey weight (Mosteller and Tukey, 1977). Both can be set by the \code{weight} argument, with respective values \code{"Huber"} and \code{"Tukey"}. Default function is Huber. Moreover, the \code{tuCo} argument specifies the tuning constant for the weight function. Default value is 1 and suggested for Huber weight (also by default), and value 4 is suggested for Tukey weight (Schuester and Yuan, 2011).

The ability level of interest (that is, for which the probability distribution must be computed) is specified by the\code{th} argument. Note that it can hold a vector of ability levels too; in this case, the probability distribution is computed for each component of \code{th}.
 }


\note{
Under the Rasch (1PL) model, the computation of the full distribution is very efficient even for a large test: due to the Lord-Wingersky algorithm, only \eqn{(n+1)} patterns (where \eqn{n} is the number of items) must be created, one for each possible test score. But for other models, \eqn{2^n} patterns must be generated and used for ability estimation, which can become computationally intensive with long tests. Therefore it is recommended not to make use of this function with more than ten items (except under the Rasch model).
}

\references{
 Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}.
 Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Birnbaum, A. (1969). Statistical theory for logistic mental test models with a prior distribution of ability. \emph{Journal of Mathematical Psychology, 6}, 258-276. \doi{10.1016/0022-2496(69)90005-4}

 Bock, R. D., and Mislevy, R. J. (1982). Adaptive EAP estimation of ability in a microcomputer environment. \emph{Applied Psychological Measurement, 6}, 431-444. \doi{10.1177/014662168200600405}

 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error}. 
 Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 
 
 Jeffreys, H. (1939). \emph{Theory of probability}. Oxford, UK: Oxford University Press.

 Jeffreys, H. (1946). An invariant form for the prior probability in estimation problems. \emph{Proceedings of the Royal Society of London. Series A, Mathematical and Physical Sciences, 186}, 453-461.

 Lord, F.M. (1980). \emph{Applications of item response theory to practical testing problems}. Hillsdale, NJ: Lawrence Erlbaum.

Lord, F. M., and Wingersky, M. S. (1984). Comparison of IRT true-score and equipercentile observed-score equatings. \emph{Applied Psychological Measurement, 8}, 453-461.  \doi{10.1177/014662168400800409}

Magis, D. and Barrada, J. R. (2017). Computerized Adaptive Testing with R: Recent Updates of the Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{Code Snippets}, \emph{76(1)}, 1-18. \doi{10.18637/jss.v076.c01}

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. \doi{10.18637/jss.v048.i08}

Mosteller, F., and Tukey, J. (1977). \emph{Exploratory data analysis and regression}. Reading, MA: Addison-Wesley.

Schuester, C., and Yuan, K.-H. (2011). Robust estimation of latent ability in item response models. \emph{Journal of Educational and Behavioral Statistics, 36}, 720)735.\doi{10.3102/1076998610396890}

 Warm, T.A. (1989). Weighted likelihood estimation of ability in item response models. \emph{Psychometrika, 54}, 427-450. \doi{10.1007/BF02294627}
}

\author{
    David Magis \cr
    \email{david.magis@gmail.com}

    Lianne Ippel \cr
    Maastricht University, The Netherlands \cr
    \email{g.j.e.ippel@gmail.com}
 }

\seealso{
\code{\link{semTheta}}, \code{\link{thetaEst}}, 
} 

\examples{

## Dichotomous models ##

 # Generation of ten items under 1PL model
 it <- genDichoMatrix(10, model = "1PL")
 
 # Full distribution of ML estimator for ability level zero
 fullDist(0, it, method = "ML")

 # Idem with BM estimator (probabilities don't change, only estimated abilities)
 fullDist(0, it, method = "BM")

 # Idem with ability level 1 (only probabilities change)
 fullDist(1, it, method = "BM")

 # Distributions with two ability levels 1 and 0.5
 fullDist(c(1, 0.5), it, method = "BM")

 # Generation of ten items under 2PL model
 it2 <- genDichoMatrix(10, model = "2PL")
 
 # Full distribution of ML estimator for ability level zero
 fullDist(0, it2, method = "ML")
 }

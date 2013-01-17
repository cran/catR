\name{createItemBank}
\alias{createItemBank}

\title{Item bank generation}

\description{
 This command creates an item bank from a matrix of item parameters. Item information functions
 are evaluated for all items and a fine grid of ability levels, to be supplied. Subgroups of items
 can also be specified for content balancing purposes.
 }

\usage{
 createItemBank(items=100, cb=FALSE, model="4PL", 
  aPrior=c("norm",1,0.2), bPrior=c("norm",0,1), 
  cPrior=c("unif",0,0.25), dPrior=c("unif",0.75,1), 
  thMin=-4, thMax=4, step=0.01, seed=1, D=1)
 }

\arguments{
 \item{items}{either an integer value or a matrix (or data frame) of item parameters (and possibly subgroups of items). See \bold{Details}.}
 \item{cb}{logical: should subgroups of items be returned in the item bank for content balancing pruposes? (default is \code{FALSE}).}
 \item{model}{character: the name of the logistic IRT model, with possible values \code{"1PL"}, \code{"2PL"}, \code{"3PL"} or \code{"4PL"}
              (default). Ignored if \code{items} is a matrix or a data frame.}
 \item{aPrior}{vector of three components, specifying the prior distribution and item parameters for generating the item discrimination levels. 
               Ignored if \code{items} is a matrix or a data frame. See \bold{Details}.}
 \item{bPrior}{vector of three components, specifying the prior distribution and item parameters for generating the item difficulty levels. 
               Ignored if \code{items} is a matrix or a data frame. See \bold{Details}.}
 \item{cPrior}{vector of three components, specifying the prior distribution and ite mparameters for generating the item lower asymptote levels. 
               Ignored if \code{items} is a matrix or a data frame. See \bold{Details}.}
 \item{dPrior}{vector of three components, specifying the prior distribution and item parameters for generating the item upper asymptote levels. 
               Ignored if \code{items} is a matrix or a data frame. See \bold{Details}.}
 \item{thMin}{numeric: the lower bound for the fine grid of ability levels (default is -4). See \bold{Details}.} 
 \item{thMax}{numeric: the upper bound for the fine grid of ability levels (default is 4). See \bold{Details}.} 
 \item{step}{numeric: the step value for the fine grid of ability levels (default is 0.01). See \bold{Details}.}
 \item{seed}{numeric: the random seed number for the generation of item parameters (default is 1). See \code{\link{set.seed}} for further details.
             Ignored if \code{items} is a matrix or a data frame.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952).} 
}

\value{
 A list of class "itBank" with four arguments:
 \item{itemPar}{the matrix of item parameters, either provided by \code{items} or generated.}
 \item{theta}{a vector with the ability levels of the fine grid, defined by arguments \code{thMin}, \code{thMax} and \code{step}.}
 \item{infoTab}{a matrix of Fisher information functions, evaluated for each ability level (one row per ability level) 
                and each item (one column per item).}
 \item{cbGroup}{either the fifth column of \code{items} (if \code{cb} was set to \code{TRUE} in addition) or \code{NULL}.}
 }
 
\details{
 If \code{items} is a matrix or a data frame, it has the following format: one row per item and four (or five) columns, with respectively 
 the discrimination \eqn{a_i}, the difficulty \eqn{b_i}, the pseudo-guessing \eqn{c_i} and the inattention \eqn{d_i}
 parameters (Barton and Lord, 1981). A fifth column can be added; in this case it holds the names of the subgroups of items (for content balancing purposes
 for instance). This column is returned in the output list only if \code{items} has five columns and \code{cb} argument is set to \code{TRUE} (default
 value is \code{FALSE}).

 If \code{items} is an integer, it corresponds to the number of items to be included in the item bank. Corresponding item parameters are drawn 
 from distributions to be specified by arguments \code{aPrior}, \code{bPrior}, \code{cPrior} and \code{dPrior} for respective parameters
 \eqn{a_i}, \eqn{b_i}, \eqn{c_i} and \eqn{d_i}. Each of these arguments is of length 3, the first component containing the name of the 
 distribution and the last two components coding the distribution parameters. Possible distributions are:
 \itemize{
 \item{the \emph{normal distribution} \eqn{N(\mu, \sigma^2)}, available for parameters \eqn{a_i} and \eqn{b_i}. It is specified by \code{"norm"} as first argument
       while the latter two arguments contain the values of \eqn{\mu} and \eqn{\sigma} respectively.}
 \item{the \emph{log-normal distribution} \eqn{\log N(\mu, \sigma^2)}, available for parameter \eqn{a_i} only. It is specified by \code{"lnorm"} as first argument
       while the latter two arguments contain the values of \eqn{\mu} and \eqn{\sigma} respectively.}
 \item{the \emph{uniform distribution} \eqn{U([a,b])}, available for all parameters. It is specified by \code{"unif"} as first argument
       while the latter two arguments contain the values of \eqn{a} and \eqn{b} respectively. Note that taking \eqn{a} and \eqn{b} equal to a common value, say
       \eqn{t}, makes all parameters to be equal to \eqn{t}.}
 \item{the \emph{Beta distribution} \eqn{Beta(\alpha, \beta)}, available for parameters \eqn{c_i} and \eqn{d_i}. It is specified by \code{"beta"} as first argument
       while the latter two arguments contain the values of \eqn{\alpha} and \eqn{\beta} respectively.}
}

 Inattention parameters \eqn{d_i} are fixed to 1 if \code{model} is not \code{"4PL"}; pseudo-guessing parameters \eqn{c_i} are 
 fixed to zero if \code{model} is either \code{"1PL"} or \code{"2PL"}; and discrimination parameters \eqn{a_i} are
 fixed to 1 if \code{model="1PL"}. The random generation of item parameters can be controlled by the \code{seed} argument.

 Note that currently it is not possible to randomly draw subgroups of items for generated sets of item parameters. Content balancing can thus be performed
 only with provided item parameters (not generated).
 
 The item bank consists of the(\code{infoTab}) matrix, which holds Fisher information functions (Baker, 1992), evaluated for each
 item in the bank and at each value of a sequence of ability levels. These abilities are ranging from \code{thMin} to \code{thMax}
 by steps of \code{step} units. 

 The returned list contains in addition the sequence of ability levels and the matrix of item parameters. If \code{items} has five 
 columns and \code{cb} argument is set to \code{TRUE}, the returned list has the additional vector \code{cbGroup} that contains the
 elements of the fifth column of \code{items} (i.e. the subgroup names). Otherwise, the \code{cbGroup} element is returned as \code{NULL}.
}

\references{
 Baker, F.B. (1992). \emph{Item response theory: parameter estimation techniques}. New York, NY: Marcel Dekker.

 Barton, M.A., and Lord, F.M. (1981). An upper asymptote for the three-parameter logistic item-response model.
 Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Haley, D.C. (1952). Estimation of the dosage mortality relationship when the dose is subject to error. 
 Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R 
 Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. URL \url{http://www.jstatsoft.org/v48/i08/}
}

\author{
    David Magis \cr
    Post-doctoral researcher, FNRS (Fonds National de la Recherche Scientifique) \cr
    Department of Mathematics, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{Ii}}
} 

\examples{
\dontrun{

 # Loading the 'tcals' parameters 
 data(tcals)
  
 # Item bank creation with 'tcals' item parameters
 createItemBank(tcals)

 # Changing the fine grid of ability levels
 createItemBank(tcals, thMin=-2, thMax=2, step=0.05)

 # Item bank creation with 500 items
 createItemBank(items=500)

 # Item bank creation with 100 items, 2PL model and log-normal 
 # distribution with parameters (0, 0.1225) for discriminations
 createItemBank(items=100, model="2PL", aPrior=c("lnorm",0,0.1225))

 # A completely identical method as for previous example
 createItemBank(items=100, aPrior=c("lnorm",0,0.1225), 
  cPrior=c("unif",0,0), dPrior=c("unif",1,1))

 # Item bank creation with 'tcals' item parameters and keeping
 # the subgroups of items
 createItemBank(tcals, cb=TRUE)
 }
 }

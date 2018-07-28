\name{genDichoMatrix}
\alias{genDichoMatrix}

\title{Item bank generation (dichotomous models)}

\description{
 This command generates an item bank from prespecified parent distributions for use with dichotomous IRT models. Subgroups of items can also be specified for content balancing purposes.
 }

\usage{
genDichoMatrix(items = 100, cbControl = NULL, model = "4PL", 
  aPrior = c("norm", 1, 0.2), bPrior = c("norm", 0, 1), 
  cPrior = c("unif", 0, 0.25), dPrior = c("unif", 0.75, 1), seed = 1)
 }

\arguments{
 \item{items}{integer: the number of items to include in the generated item bank.}
 \item{cbControl}{either a list to define subgroups for content balancing or \code{NULL} (default). See \bold{Details}.}
 \item{model}{character: the name of the logistic IRT model, with possible values \code{"1PL"}, \code{"2PL"}, \code{"3PL"} or \code{"4PL"} (default).}
 \item{aPrior}{vector of three components, specifying the prior distribution and item parameters for generating the item discrimination levels. See \bold{Details}.}
 \item{bPrior}{vector of three components, specifying the prior distribution and item parameters for generating the item difficulty levels. See \bold{Details}.}
 \item{cPrior}{vector of three components, specifying the prior distribution and item parameters for generating the item lower asymptote levels. See \bold{Details}.}
 \item{dPrior}{vector of three components, specifying the prior distribution and item parameters for generating the item upper asymptote levels. See \bold{Details}.}
 \item{seed}{numeric: the random seed number for the generation of item parameters (default is 1). See \code{\link{set.seed}} for further details.}
}

\value{
 A data frame with four or five arguments:
 \item{\code{a}}{the generated item discrimination parameters.}
 \item{\code{b}}{the generated item difficulty parameters.}
 \item{\code{c}}{the generated item lower asymptote parameters.}
 \item{\code{d}}{the generated item upper asymptote parameters.}
 \item{\code{Group}}{(optional) the distribution of subgroup names across items. Ignored if \code{cbControl} is \code{NULL}.}
 }
 
\details{
This function permits to generate an item bank under dichotomous IRT models that is compatible for use with \code{\link{randomCAT}}. 

The number of items to be included in the bank is specified by the \code{items} argument. Corresponding item parameters are drawn from distributions to be specified by arguments \code{aPrior}, \code{bPrior}, \code{cPrior} and \code{dPrior} for respective parameters \eqn{a_i}, \eqn{b_i}, \eqn{c_i} and \eqn{d_i} (Barton and Lord, 1981). Each of these arguments is of length 3, the first component containing the name of the distribution and the last two components coding the distribution parameters. 

Possible distributions are:
 \itemize{
 \item{the \emph{normal distribution} \eqn{N(\mu, \sigma^2)}, available for parameters \eqn{a_i} and \eqn{b_i}. It is specified by \code{"norm"} as first argument while the latter two arguments contain the values of \eqn{\mu} and \eqn{\sigma} respectively.}
 \item{the \emph{log-normal distribution} \eqn{\log N(\mu, \sigma^2)}, available for parameter \eqn{a_i} only. It is specified by \code{"lnorm"} as first argument while the latter two arguments contain the values of \eqn{\mu} and \eqn{\sigma} respectively.}
 \item{the \emph{uniform distribution} \eqn{U([a,b])}, available for all parameters. It is specified by \code{"unif"} as first argument while the latter two arguments contain the values of \eqn{a} and \eqn{b} respectively. Note that taking \eqn{a} and \eqn{b} equal to a common value, say \eqn{t}, makes all parameters to be equal to \eqn{t}.}
 \item{the \emph{Beta distribution} \eqn{Beta(\alpha, \beta)}, available for parameters \eqn{c_i} and \eqn{d_i}. It is specified by \code{"beta"} as first argument while the latter two arguments contain the values of \eqn{\alpha} and \eqn{\beta} respectively.}
}
Inattention parameters \eqn{d_i} are fixed to 1 if \code{model} is not \code{"4PL"}; pseudo-guessing parameters \eqn{c_i} are fixed to zero if \code{model} is either \code{"1PL"} or \code{"2PL"}; and discrimination parameters \eqn{a_i} are
 fixed to 1 if \code{model="1PL"}. The random generation of item parameters can be controlled by the \code{seed} argument.

If required, the distribution of the items across subgroups with specified names can be performed. To do so, the \code{cbControl} argument must be supplied with a list of two arguments: (a) the first argument is called \code{$names} and contains the different names of the subgroups of items; (b) the second argument is called \code{$props} and contains a vector of numeric values, of the same length of \code{names} element, with only positive numbers but not necessarily summing to one. For instance, if \code{props} is set as \code{c(1, 2, 2)} and \code{items} to \code{100}, then the three subgroups will hold respectively 20, 40 and 40 items.

Several constraints apply to the arguments of the \code{cbControl} list. First, both arguments of \code{cbControl} must be of the same length. Second, as already explained, \code{cbControl$props} must either sum to 1 (in case of proportions) or to \code{items} (in case of integer values). Finally, if proportions are provided to \code{cbControl$props}, one can ensure that when multiplied by \code{items} they return integer values (so that they can sum up to \code{items}).

The random generation of item parameters and the random allocation of items to subgroups of items are both under control by the \code{seed} argument.

The output is a data frame with at least four arguments, with names \code{a}, \code{b}, \code{c} and \code{d} for respectively the discrimination \eqn{a_i}, the difficulty \eqn{b_i}, the lower asymptote \eqn{c_i} and the upper asymptote \eqn{d_i} parameters. A fifth argument contains optionally the subgroup names that have been randomly assigned to the generated items, in the proportions specified by the \code{$props} argument of the \code{cbControl} list. 
}


\note{
The current version of \code{genItemBank} is only designed for dichotomous IRT models. Future extensions will hopefully provide the same tool for polytomous IRT models.
}


\references{
 Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}. Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

Magis, D. and Barrada, J. R. (2017). Computerized Adaptive Testing with R: Recent Updates of the Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{Code Snippets}, \emph{76(1)}, 1-18. \doi{10.18637/jss.v076.c01}

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. \doi{10.18637/jss.v048.i08}
}

\author{
    David Magis \cr
    \email{david.magis@gmail.com}
 }

\seealso{
\code{\link{nextItem}}, \code{\link{randomCAT}}
} 

\examples{

 # Item bank generation with 500 items
 genDichoMatrix(items = 500)

 # Item bank generation with 100 items, 2PL model and log-normal distribution with 
 # parameters (0, 0.1225) for discriminations
 genDichoMatrix(items = 100, model = "2PL", aPrior = c("lnorm", 0, 0.1225))

 # A completely identical method as for previous example
 genDichoMatrix(items = 100, aPrior = c("lnorm", 0, 0.1225), 
  cPrior = c("unif", 0, 0), dPrior = c("unif", 1, 1))

 # Item bank generation with prespecified content balancing control options 
 cbList <- list(names = c("Group1", "Group2", "Group3", "Group4"), 
        props = c(0.2,0.4,0.3,0.1))
 genDichoMatrix(items = 100, cbControl = cbList)

 # With proportions that do not sum to one
 cbList <- list(names = c("Group1", "Group2", "Group3", "Group4"), props=c(2, 4, 3, 1))
 genDichoMatrix(items = 100, cbControl = cbList)

 }

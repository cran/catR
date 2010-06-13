\name{createItemBank}
\alias{createItemBank}

\title{Item bank generation}

\description{
 This command creates an item bank from a matrix of item parameters. Item information functions
 are evaluated for all items and a fine grid of ability levels, to be supplied. 
 }

\usage{
 createItemBank(items=100, model="4PL", thMin=-4, thMax=4,
	 step=0.01, seed=1)
 }

\arguments{
 \item{items}{either an integer value or a matrix of item parameters. See \bold{Details}.}
 \item{model}{character: the name of the logistic IRT model, with possible values \code{"1PL"}, \code{"2PL"}, \code{"3PL"} or \code{"4PL"}
              (default). Ignored if \code{items} is a matrix.}
 \item{thMin}{numeric: the lower bound for the fine grid of ability levels (default is -4). See \bold{Details}.} 
 \item{thMax}{numeric: the upper bound for the fine grid of ability levels (default is 4). See \bold{Details}.} 
 \item{step}{numeric: the step value for the fine grid of ability levels (default is 0.01). See \bold{Details}.}
 \item{seed}{numeric: the random seed number for the generation of item parameters (default is 1). See \code{\link{set.seed}} for further details.}
}

\value{
 A list of class "itBank" with three arguments:
 \item{itemPar}{the matrix of item parameters, either provided by \code{items} or generated.}
 \item{theta}{a vector with the ability levels of the fine grid, defined by arguments \code{thMin}, \code{thMax} and \code{step}.}
 \item{infoTab}{a matrix of item information functions, evaluated for each ability level (one row per ability level) 
                and each item (one column per item).}
 }
 
\details{
 If \code{items} is a matrix, it has the following format: one row per item and four columns, with respectively 
 the discrimination \eqn{a_i}, the difficulty \eqn{b_i}, the pseudo-guessing \eqn{c_i} and the inattention \eqn{d_i}
 parameters (Barton and Lord, 1981). If \code{items} is an integer, it corresponds to the number of items to be 
 included in the item bank. Corresponding item parameters are (by default) randomly drawn from the following distributions:
 \eqn{a_i \sim N(1, 0.2^2)}, \eqn{b_i \sim N(0, 1)}, \eqn{c_i \sim U([0, 0.25])} and \eqn{d_i \sim U([0.75, 1])}. Inattention
 parameters \eqn{d_i} are fixed to 1 if \code{model} is not \code{"4PL"}; pseudo-guessing parameters \eqn{c_i} are 
 fixed to zero if \code{model} is either \code{"1PL"} or \code{"2PL"}; and discrimination parameters \eqn{a_i} are
 fixed to 1 if \code{model="1PL"}. The random generation of item parameters can be controlled by the \code{seed} argument.
 
 The item bank consists of a matrix of item information functions (Baker, 1992), evaluated for each itelm in the bank
 and at each value of a sequence of ability levels. These abilities are ranging from \code{thMin} to \code{thMax} by 
 steps of \code{step} units. 

 The returned list contains in addition the sequence of ability levels and the matrix of item parameters.
}

\references{
 Baker, F.B. (1992). \emph{Item response theory: parameter estimation techniques}. New York, NY: Marcel Dekker.

 Barton, M.A., and Lord, F.M. (1981). An upper asymptote for the three-parameter logistic item-response model.
 Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.
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
 tcals <- as.matrix(tcals)
 
 # Item bank creation with 'tcals' item parameters
 createItemBank(tcals)

 # Changing the fine grid of ability levels
 createItemBank(tcals, thMin=-2, thMax=2, step=0.05)

 # Item bank creation with 500 items
 createItemBank(items=500)
 }
 }

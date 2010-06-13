\name{startItems}
\alias{startItems}

\title{Selection of the first items}

\description{
 This command selects the first items of the adaptive test, either randomly or on the basis of their difficulty level.
 }

\usage{
 startItems(it, fixItems=NULL, seed=NULL, nrItems=1,
 	 theta=0, bw=4, range=c(-4,4))
 }

\arguments{
 \item{it}{numeric: a matrix with one row per item and four columns, with the values of the discrimination, the difficulty,
           the pseudo-guessing and the inattention parameters (in this order).}
 \item{fixItems}{either a vector of integer values or \code{NULL} (default). See \bold{Details}.}
 \item{seed}{either a numeric value or \code{NULL} (default). See \bold{Details}.}
 \item{nrItems}{numeric: the number of starting items to be selected (default is 1). Ignored if \code{fixItems} is not \code{NULL}.}
 \item{theta}{numeric: the initial ability level for selecting the first items (default is 0). 
              Ignored if either \code{fixItems} or \code{seed} is not \code{NULL}. See \bold{Details}.}
 \item{bw}{numeric: the length (bandwidth) of the difficulty interval for selecting the first items (default is 4).
          Ignored if either \code{fixItems} or \code{seed} is not \code{NULL}. See \bold{Details}.}
 \item{range}{numeric: the range of diffiulty levels for selecting the first items (default is \code{c(-4,4)}).
              Ignored if either \code{fixItems} or \code{seed} is not \code{NULL}. See \bold{Details}.}
}

\value{
 A list with three arguments:
 \item{items}{the selected items (identified by their number in the item bank).}
 \item{par}{the matrix of item parameters of the selected items (one row per item).}
 \item{bOpt}{the values of the item difficulties chosen for selecting the first items.}
 }
 

\details{
 The first items of the adaptive test can be selected either by fixing these items or by some selection into the item bank. 
 The first items can be fixed by specifying \code{fixItems} as a vector of integer values, each value referring to the corresponding
 item in the item bank. 

 If \code{fixItems} is kep to \code{NULL} (the default choice), then \code{nrItems} will be selected, either randomly or
 according to their difficulty level. In the former, one has to assign a numeric value to the \code{seed} argument. In the latter,
 \code{seed} is kept to \code{NULL} (the default value) and the items are selected as follows. The arguments \code{theta} and
 \code{bw} permit to define a sequence of value, ranging from \code{theta} - \code{bw} to \code{theta} + \code{bw} and of
 length \code{nrItems} (if \code{nrItems} equals 1, the range reduces to \code{theta}). Then, for each value in the sequence,
 the function selects the item in the item bank with difficulty level as close as possible to the sequential value. By default,
 \code{theta} equals 0 and \code{bw} equals 4, meaning that the sequence ranges from -4 to 4 and of length \code{nrItems}.

 The argument \code{range} permits to constraint the range of difficulty levels, by default to the range [-4, 4].
 }



\author{
    David Magis \cr
    Post-doctoral researcher, FNRS (Fonds National de la Recherche Scientifique) \cr
    Department of Mathematics, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }


\examples{
\dontrun{

 # Loading the 'tcals' parameters 
 data(tcals)
 tcals <- as.matrix(tcals)
 
 # Selecting 1 starting item, initial ability estimate is 0
 startItems(tcals) 

 # Selecting 3 starting items, initial ability estimate is 1
 # and bandwidth is 2
 startItems(tcals, nrItems=3, theta=1, bw=2) 

 # Selecting 5 starting items, initial ability estimate is 2
 # and bandwidth is 3
 startItems(tcals, nrItems=5, theta=2, bw=3) 

 # Random selection of 4 starting items
 startItems(tcals, seed=1, nrItems=4)

 # Selection of the first 5 starting items
 startItems(tcals, fixItems=1:5)
 }
 }

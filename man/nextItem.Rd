\name{nextItem}
\alias{nextItem}

\title{Selection of the next item}

\description{
 This command selects the next item to be administered, given the list of previously administered items and the current
 ability estimate, with the 'maximum information' criterion. 
 }

\usage{
 nextItem(itemBank, theta, out=NULL)
 }

\arguments{
 \item{itemBank}{an item bank of class \code{itBank} sa output of the function \code{\link{createItemBank}}.}
 \item{theta}{numeric: the current value of the ability estimate.}
 \item{out}{either a vector of integer values specifying the items previously administered, or \code{NULL} (default).} 
}

\value{
 A list with three arguments:
 \item{item}{the selected item (identified by its number in the item bank).}
 \item{par}{the vector of item parameters of the selected item.}
 \item{info}{the value of the information function for the selected item and the current ability estimate.}
 }
 

\details{
 Currently only one ethod is available for selecting the next item to be administered in the adaptive test. For a given
 current ability estimate, the next item is selected (among the available items) as the one which maximizes the item
 information function (Baker, 1992). The most informative item is selected from the table of item informations provided
 by the bank of items specified with \code{itemBank}. 

 The available items are those that are not specified in the \code{out} argument. By default, \code{out} is \code{NULL},
 which means that all items are available.  
}

\references{
 Baker, F.B. (1992). \emph{Item response theory: parameter estimation techniques}. New York, NY: Marcel Dekker.
}

\author{
    David Magis \cr
    Post-doctoral researcher, FNRS (Fonds National de la Recherche Scientifique) \cr
    Department of Mathematics, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{createItemBank}}
} 

\examples{
\dontrun{

 # Loading the 'tcals' parameters 
 data(tcals)
 tcals <- as.matrix(tcals)
 
 # Item bank creation with 'tcals' item parameters
 bank <- createItemBank(tcals)

 # Selecting the next item, current ability estimate is 0
 nextItem(bank, 0) # item 63 is selected

 # Selecting the next item, current ability estimate is 0
 # and item 63 is removed
 nextItem(bank, 0, out=63) # item 10 is selected

 # Selecting the next item, current ability estimate is 0
 # and items 63 and 10 are removed
 nextItem(bank, 0, out=c(63,10)) # item 62 is selected
 }
 }

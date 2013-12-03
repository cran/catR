\name{test.cbList}
\alias{test.cbList}

\title{Testing the format of the list for content balancing}

\description{
 This command tests whether format of the list to be assigned to \code{cbControl} argument of function \code{\link{nextItem}} 
 is appropriate for content balancing, and returns a warning message otherwise.
 }

\usage{
 test.cbList(list, itemBank)
 }

\arguments{
 \item{list}{a list of arguments to be tested. See \bold{Details}.}
 \item{itemBank}{an item bank of class \code{itBank} as output of the function \code{\link{createItemBank}}, holding subgroups of items (i.e.
                 such that \code{itemBank$cbGroup} is not \code{NULL}.}
}

\value{
 A list with two arguments:
 \item{test}{a logical value indicating whether the format of the list is accurate (\code{TRUE}) or not (\code{FALSE}).}
 \item{message}{either a message to indicate the type of misspecification, or \code{"ok"} if the format is accurate.}
 }
 
\details{
 The \code{test.cbList} function checks whether the list provided in the \code{cbControl} argument of the \code{\link{nextItem}} and
 \code{\link{randomCAT}} functions, is accurate for controlling for ontent balancing. It mainly serves as an initial check for the
 \code{\link{randomCAT}} function. 

 The function returns an \code{"ok"} message if the arguments of \code{list} match the requirement of the list \code{cbControl} for
 content balancing. Otherwise, a message is returned with information about list - type mismatch. This will be the case
 if:
 \itemize{
 \item{\code{list} is not a list or has not exactly two elements,}
 \item{at least one of the argument names is incorrect,}
 \item{the lengths of the arguments are different, or different from the number of subgroups of items,}
 \item{the 'names' element does not match with the subgroups' names,}
 \item{the 'props' element is not numeric or holds negative values.}
 }
 Each mismatch yields a different output message to help in debugging the problem.
 }

\references{
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
\code{\link{nextItem}}, \code{\link{randomCAT}}
} 

\examples{
\dontrun{

 # Loading the 'tcals' parameters 
 data(tcals)
 
 # Item bank creation for content balancing
 bank <- createItemBank(tcals, cb=TRUE)

 # Creation of a correct list with arbitrary proportions
 cbList <- list(names=c("Audio1","Audio2","Written1","Written2",
                "Written3"), props=c(0.1,0.2,0.2,0.2,0.3))

 # Testing 'cbList'
 test.cbList(cbList, bank)

 # Creation of an incorrect list (mismatch in first name)
 cbList <- list(names=c("audio1","Audio2","Written1","Written2",
                "Written3"), props=c(0.1,0.2,0.2,0.2,0.3))
 test.cbList(cbList, bank)

 # Creation of an incorrect list (mismatch in name of second
 # element)
 cbList <- list(names=c("Audio1","Audio2","Written1","Written2",
                "Written3"), prop=c(0.1,0.2,0.2,0.2,0.3))
 test.cbList(cbList, bank)

 # Creation of an incorrect list (second element shorter than
 # first element)
 cbList <- list(names=c("Audio1","Audio2","Written1","Written2",
                "Written3"), props=c(0.1,0.2,0.2,0.2))
 test.cbList(cbList, bank)

 # Creation of an incorrect list (adding a third element)
 cbList <- list(names=c("Audio1","Audio2","Written1","Written2",
                "Written3"), props=c(0.1,0.2,0.2,0.2), third="hi")
 test.cbList(cbList, bank)
 }
 }

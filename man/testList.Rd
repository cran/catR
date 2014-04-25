\name{testList}
\alias{testList}

\title{Testing the format of the input lists}

\description{
 This command tests whether format of the input lists for the random generation of adaptive tests is convenient,
 and returns a warning message otherwise.
 }

\usage{
testList(list, type = "start")
 }

\arguments{
 \item{list}{a list of arguments to be tested. See \bold{Details}.}
 \item{type}{character: the type of list for checking. Possible values are \code{"start"} (default), \code{"test"},
             \code{"stop"} and \code{"final"}. See \bold{Details}.}
}

\value{
 A list with two arguments:
 \item{test}{a logical value indicating whether the format of the list is accurate (\code{TRUE}) or not (\code{FALSE}).}
 \item{message}{either a message to indicate the type of misspecification, or \code{"ok"} if the format is accurate.}
 }
 

\details{
 The \code{testList} function checks whether the list provided in the \code{list} argument is accurate for the selected \code{type}.
 It mainly serves as an initial check for the \code{\link{randomCAT}} function. 

 The four types of lists are: \code{"start"} with the parameters for selecting the first items; \code{"test"} with
 the options of the adaptive test (i.e. method for next item selection, provisional ability estimator and related information);
 \code{"stop"} with the options setting the stopping rule; and \code{"final"} with the options for final ability estimation. 
 See the help file of \code{\link{randomCAT}} for further details about the different lists, their allowed arguments and their contents.

 The function returns an \code{"ok"} message if the arguments of \code{list} match the requirement of the
 corresponding \code{type}. Otherwise, a message is returned with information about list - type mismatch. This will be the case
 if:
 \itemize{
 \item{\code{list} is not a list, or has no argument names,}
 \item{\code{list} has too many arguments for the \code{type} specified,}
 \item{at least one of the argument names is incorrect,}
 \item{the content of at least one argument is not adequate (e.g. character instead of numeric).}
 }
 Each mismatch yields a different output message to help in debugging the problem.
 }

\references{
 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R 
 Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. URL \url{http://www.jstatsoft.org/v48/i08/}
}

\author{
    David Magis \cr
    Department of Education, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{randomCAT}}
} 

\examples{
\dontrun{

 # Creation and test of a 'start' list
 start <- list(nrItems = 3, theta = 0, halfRange = 2)
 testList(start, type = "start")

 # Modification of the list to introduce a mistake
 names(start)[1] <- "nrItem"
 testList(start, type = "start")
 
 # Creation and test of a 'test' list
 test <- list(method = "WL", itemSelect = "bOpt")
 testList(test, type = "test")

# Creation and test of a 'stop' list
 stop <- list(method = "WL")
 testList(stop, type = "test")

# Creation and test of a 'final' list (with mistake)
 final <- list(method = "MAP")
 testList(final, type = "final")
 }
 }

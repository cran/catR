\name{test.cbList}
\alias{test.cbList}

\title{Testing the format of the list for content balancing under dichotomous or polytomous IRT models}

\description{
 This command tests whether format of the list to be assigned to \code{cbControl} argument of function \code{\link{nextItem}} is appropriate for content balancing, and returns a warning message otherwise.
 }

\usage{
test.cbList(list, cbGroup)
 }

\arguments{
 \item{list}{a list of arguments to be tested. See \bold{Details}.}
 \item{cbGroup}{a vector of character or factor names of the subgroups of items in the bank, or \code{NULL}.}
}

\value{
 A list with two arguments:
 \item{test}{a logical value indicating whether the format of the list is accurate (\code{TRUE}) or not (\code{FALSE}).}
 \item{message}{either a message to indicate the type of misspecification, or \code{"ok"} if the format is accurate.}
 }
 
\details{
 The \code{test.cbList} function checks whether the list provided in the \code{cbControl} argument of the \code{\link{nextItem}} and \code{\link{randomCAT}} functions, is accurate for controlling for content balancing. It mainly serves as an initial check for the \code{\link{randomCAT}} function. 

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

## Dichotomous models ##

 # Loading the 'tcals' parameters 
 data(tcals)
 
 # Vector of group membership for 'tcals' data set
 group <- tcals$Group

 # Creation of a correct list with arbitrary proportions
 cbList <- list(names = c("Audio1", "Audio2", "Written1", "Written2", "Written3"), 
        props = c(0.1, 0.2, 0.2, 0.2, 0.3))

 # Testing 'cbList'
 test.cbList(cbList, group)

 # Creation of an incorrect list (mismatch in first name)
 cbList <- list(names = c("audio1", "Audio2", "Written1", "Written2", "Written3"), 
        props=c(0.1, 0.2, 0.2, 0.2, 0.3))
 test.cbList(cbList, group)

 # Creation of an incorrect list (mismatch in name of second
 # element)
 cbList <- list(names = c("Audio1", "Audio2", "Written1", "Written2", "Written3"), 
        prop = c(0.1, 0.2, 0.2, 0.2, 0.3))
 test.cbList(cbList, group)

 # Creation of an incorrect list (second element shorter than
 # first element)
 cbList <- list(names = c("Audio1", "Audio2", "Written1", "Written2", "Written3"), 
         props=c(0.1, 0.2, 0.2, 0.2))
 test.cbList(cbList, group)

 # Creation of an incorrect list (adding a third element)
 cbList <- list(names = c("Audio1", "Audio2", "Written1", "Written2", "Written3"),  
        props = c(0.1, 0.2, 0.2, 0.2), third = "hi")
 test.cbList(cbList, group)


## Polytomous models ##

 # Creation of an appropriate list for content balancing
 # Equal proprotions across subgroups of items
 cbList <- list(names = c("Group1", "Group2", "Group3", "Group4"), props = rep(1, 4))

 # Creation of a "wrong" list
 cbList2 <- list(names=c("group1","group2"),props = c(1, 1))

 # Generation of an item bank under GRM with 100 items and at most 4 categories
 m.GRM <- genPolyMatrix(100, 4, "GRM", cbControl = cbList)

 # Breaking the 'm.GRM' bank
 bank <- breakBank(m.GRM)

 # Testing 'cbList' and 'cbList2'
 test.cbList(cbList, bank$cbGroup)
 test.cbList(cbList2, bank$cbGroup)

 # Generation of an item bank under PCM with 100 items, 4 categories and groups
 m.PCM <- genPolyMatrix(100, 4, "PCM", same.nrCat = TRUE, cbControl = cbList2)
 
 # Breaking the 'm.PCM' bank
 bank2 <- breakBank(m.PCM)

 # Testing 'cbList' and 'cbList2'
 test.cbList(cbList, bank2$cbGroup)
 test.cbList(cbList2, bank2$cbGroup)
 }

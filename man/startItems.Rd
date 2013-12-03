\name{startItems}
\alias{startItems}

\title{Selection of the first items}

\description{
 This command selects the first items of the adaptive test, either randomly or on the basis of their difficulty level.
 }

\usage{
 startItems(itemBank, fixItems=NULL, seed=NULL, nrItems=1,
 	 theta=0, halfRange=2, startSelect="bOpt") 
 }

\arguments{
 \item{itemBank}{an item bank of class \code{itBank} as output of the function \code{\link{createItemBank}}.}
 \item{fixItems}{either a vector of integer values or \code{NULL} (default). See \bold{Details}.}
 \item{seed}{either a numeric value or \code{NULL} (default).Ignored if \code{fixItems} is not \code{NULL}. See \bold{Details}.}
 \item{nrItems}{numeric: the number of starting items to be selected (default is 1). Ignored if \code{fixItems} is not \code{NULL}.}
 \item{theta}{numeric: the initial ability level for selecting the first items (default is 0). 
              Ignored if either \code{fixItems} or \code{seed} is not \code{NULL}. See \bold{Details}.}
 \item{halfRange}{numeric: the half of the range of initial ability values (default is 2). Ignored if either
                 \code{fixItems} or \code{seed} is not \code{NULL}. See \bold{Details}.}
 \item{startSelect}{character: the criterion for selecting the first items. Possible values are \code{"bOpt"}
                  (default) and \code{"MFI"}. See \bold{Details}.}
}

\value{
 A list with four arguments:
 \item{items}{the selected items (identified by their number in the item bank).}
 \item{par}{the matrix of item parameters of the selected items (one row per item).}
 \item{thStart}{the sequence of starting ability values used for selecting the items.}
 \item{startSelect}{the value of the \code{startSelect} argument.}
 }
 

\details{
 This function permits to select the first items of the test. The number of starting items is given by the \code{nrItems} argument,
 with default value 1. 

 The first item(s) of the adaptive test can be selected by one of the following methods.
 \enumerate{
 \item{By specifying the item(s) to be administered. The argument \code{fixItems} then holds the item number(s) as 
       listed in the item bank. Setting \code{fixItems} to \code{NULL} (default value) disables this method.}
 \item{By selecting it (them) randomly into the item bank. The argument \code{seed} permits to fix the random selection by 
       specifying the random seed number. Setting \code{seed} to \code{NULL} (default value) disables this method.}
 \item{By selecting the item(s) according to an initial sequence of ability values (see below). In this case, two criteria
       can be used: either one selects the item(s) whose difficulty level is as close as possible to the inital ability 
       value(s), or one selects the most informative item(s) for the given initial ability value(s). The criterion is
       specified by the \code{startSelect} argument, with values \code{"bOpt"} (default) for the 'difficulty' criterion, 
       and \code{"MFI"} for the 'information' criterion.}
 }
 The third method above will be used if and only if both \code{fixItems} and \code{seed} arguments are fixed to \code{NULL}.
 Otherwise, one of the first two methods will be used (see also \code{\link{testList}} for details about debugging 
 misspecifications of the starting arguments).

 The sequence of initial ability estimates is specified by the triplet of arguments (\code{nrItems}, \code{theta},
 \code{halfRange}). As mentioned above, \code{nrItems} is the number of items to select, and thus the length of the 
 sequence. The \code{theta} value is the central ability value, and \code{halfRange} sets half of the range of the 
 ability values. These three arguments altogether permit to define any type of (equidistant) ability values. For 
 instance, 
 \itemize{
 \item{the set (-1, 1) can be obtained by specifying the triplet to (2, 0, 1);}
 \item{the set (-1, 0, 1) can be obtained by specifying the triplet to (3, 0, 1);}
 \item{the set (-1, 0, 1, 2) can be obtained by specifying the triplet to (4, 0.5, 1.5);}
 \item{etc.}
 }  
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
\code{\link{createItemBank}}, \code{\link{testList}}
} 

\examples{
\dontrun{

 # Loading the 'tcals' parameters 
 data(tcals)

 # Item bank creation with 'tcals' item parameters
 bank <- createItemBank(tcals)
 
 # Random selection of 4 starting items
 startItems(bank, seed=1, nrItems=4)

 # Selection of the first 5 starting items
 startItems(bank, fixItems=1:5)

 # Selecting 1 starting item, initial ability estimate is 0
 startItems(bank) 

 # Selecting 3 starting items, initial ability estimate is 1
 # and half range is 2
 startItems(bank, nrItems=3, theta=1, halfRange=2) 

 # Idem but with 'information' criterion
 startItems(bank, nrItems=3, theta=1, halfRange=2, startSelect="MFI")
 
 # Selecting 5 starting items, initial ability estimate is 2
 # and half range is 3
 startItems(bank, nrItems=5, theta=2, halfRange=3) 
 }
 }

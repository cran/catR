\name{startItems}
\alias{startItems}

\title{Selection of the first items}

\description{
 This command selects the first items of the adaptive test, either randomly or on the basis of their information function.
 }

\usage{
startItems(itemBank, model = NULL, fixItems = NULL, seed = NULL, nrItems = 1, 
 	 theta = 0, D = 1, randomesque = 1, random.seed = NULL, startSelect = "MFI",
 	 nAvailable = NULL, cbControl = NULL, cbGroup = NULL, random.cb=NULL) 
 }

\arguments{
 \item{itemBank}{numeric: a suitable matrix of item parameters. See \bold{Details}.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
 \item{fixItems}{either a vector of integer values or \code{NULL} (default). See \bold{Details}.}
 \item{seed}{either a numeric value, \code{NA} or \code{NULL} (default). Ignored if \code{fixItems} is not \code{NULL}. See \bold{Details}.}
 \item{nrItems}{numeric: the number of starting items to be randomly selected (default is 1). Can be equal to zero to avoid initial selection of items (see \bold{Details}). Used only if \code{fixItems} is \code{NULL} and \code{seed} is not \code{NULL}.}
 \item{theta}{numeric: a vector of the initial ability levels for selecting the first items (default is the single value 0). Ignored if either \code{fixItems} or \code{seed} is not \code{NULL}. See \bold{Details}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952). Ignored if \code{model} is not \code{NULL} and if \code{startSelect} is not \code{"MFI"}.} 
 \item{randomesque}{integer: the number of 'randomesque' items to be picked up optimally for each value of \code{theta} vector, before random selection of a single one. Ignored if either \code{fixItems} or \code{seed} is not \code{NULL}. See \bold{Details}.}
\item{random.seed}{either \code{NULL} (default) or a numeric value to fix the random seed of randomesque selection of the items. Ignored if either \code{fixItems} or \code{seed} is not \code{NULL}.}
 \item{startSelect}{character: the criterion for selecting the first items. Possible values are \code{"bOpt"}, \code{"thOpt"}, \code{"progressive"}, \code{"proportional"}, and \code{"MFI"} (default). See \bold{Details}.}
 \item{nAvailable}{either a boolean vector indicating which items (denoted by 1's) are available at the start of the test and which (denoted by 0's) are not, or \code{NULL} (default). See \bold{Details}.}
 \item{cbControl}{either a list of accurate format to control for content balancing, or \code{NULL}. See \bold{Details}.}
 \item{cbGroup}{either a factor vector of accurate format to control for content balancing, or \code{NULL}. See \bold{Details}.}
\item{random.cb}{either \code{NULL} (default) or a numeric value to fix the selection of subgroups of items for random sampling the starting items. See \bold{Details}.}
}

\value{
 A list with five arguments:
 \item{items}{the selected items (identified by their number in the item bank) or \code{NULL} (if \code{nrItems} is 0).}
 \item{par}{the matrix of item parameters of the selected items (one row per item) or \code{NULL} (if \code{nrItems} is 0).}
 \item{thStart}{the sequence of starting ability values used for selecting the items or \code{NA} (if not applicable) or \code{NULL} (if \code{nrItems} is 0).}
 \item{startSelect}{the value of the \code{startSelect} argument or \code{NA} (if not applicable) or \code{NULL} (if \code{nrItems} is 0).}
\item{names}{either a vector with the names of the selected item(s) or \code{NULL} (if the item bank has no item names).}
 }
 

\details{
 This function permits to select the first item(s) of the test. It works with both dichotomous and polytomous item banks.

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{it} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model
 (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model, \code{"PCM"} for Partical Credit Model, \code{"GPCM"} for Generalized Partial Credit Model, \code{"RSM"} for Rating Scale Model and \code{"NRM"} for Nominal Response Model. The \code{it} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.

The number of starting items is given by the length of \code{fixItems} argument, the \code{nrItems} argument (in case of random selection) or by the length of \code{theta} argument (in case of optimal selection), with default value 1 in all cases. It can be set to zero; in this case, only \code{NULL} values are returned in the output list and the CAT process will start without starting items. 

 The first item(s) of the adaptive test can be selected by one of the following methods.
 \enumerate{
 \item{By specifying the item(s) to be administered. The argument \code{fixItems} then holds the item number(s) as listed in the item bank. Setting \code{fixItems} to \code{NULL} (default value) disables this method.}
 \item{By selecting it (them) randomly into the item bank. The argument \code{seed} permits to fix the random selection by 
specifying the random seed number and the number of selected items is fixed by the \code{nrItems} argument. Setting \code{seed} to \code{NA} disables the random seed (though items are still picked up randomly in the bank); in other words, successive runs of \code{startItems} with \code{seed=NA} may lead to different item(s) selection. Setting \code{seed} to \code{NULL} (default value) disables this selection method.}
 \item{By selecting the item(s) according to an initial sequence of ability values set by the \code{theta} argument. In this case, five criteria can be used, specified through the \code{startSelect} argument:  
\enumerate{
\item{\code{"MFI"} (default): one selects the most informative item(s) for the given initial ability value(s);}
\item{\code{"bOpt"}: one selects the item(s) whose difficulty level is as close as possible to the inital ability value(s);}
\item{\code{"thOpt"}: one selects the item(s) with the ability value where they get their maximum Fisher information is as close as possible to the inital ability value(s) (see Magis, 2013, for further details);}
\item{\code{"progressive"} for the \emph{progressive} method (see \code{\link{nextItem}});}
\item{\code{"proportional"} for the \emph{proportional} method (see \code{\link{nextItem}}).}
}
If the \code{"progressive"} or \code{"proportional"} methods} are selected, this will force the values of \code{fixItems} to \code{NULL}, \code{seed} to \code{NULL}, and \code{nrItems} to 1. Thus, a single item will be selected randomly.}
 
 The third method above will be used if and only if both \code{fixItems} and \code{seed} arguments are fixed to \code{NULL}.
 Otherwise, one of the first two methods will be used (see also \code{\link{testList}} for details about debugging 
 misspecifications of the starting arguments).

 The sequence of initial ability estimates is specified by the argument \code{theta}. For each component of \code{theta} one item will be picked up optimally according to the chosen \code{startSelect} argument. however, it is possible to perform randomesque selection at this stage by setting the \code{randomesque} argument to an integer value larger than one. In this case, the \code{randomesque} most optimal items are chosen per value of \code{theta}, and the final starting item is chosen randomly among the \code{randomesque} items. By default, only one item is picked up per ability level (and is therefore \emph{the} most optimal one).

Only part of the full item bank can be made available for the selection of the first item(s), while others can be dropped out from this first step. This is fixed by the \code{nAvailable} argument, which is a vector with as many components as items in the bank and with zeros and ones only. Values 1 code for available items, values 0 for non-available items. By default, \code{nAvailable} is \code{NULL} and all items are available. Note that \code{nrItems} should never be larger than the number of available items (i.e. \code{sum(nAvailable)}). Otherwise an error message is returned.

Finally, in case of random selection of the first item(s), it is possible to force this selection to approximately match the content balancing control options set by arguments \code{cbControl}  and \code{cbGroup}. Detailed description of these arguments (and their accurate format) are described in \code{\link{nextItem}} function. Note that \code{cbControl} should be tested with the \code{\link{test.cbList}} function prior to using \code{\link{startItems}}.

In practice, when fixing both the number of items to select (through \code{nrItems}) and the proportions of items per subgroups (through \code{cbControl$props}), the selection is made in two steps. First, the number of items per each subgroup is determined in order to reach closest distribution to the theoretical one. Second, within each subgroup the required number of items are randomly drawn (this draw can be fixed with the \code{seed} argument; fixing it to \code{NA} withdraws this option).

In the first step, it is possible that several subgroups are equally proable for selection; in this case the selection of the subgroups can be fixed by the argument \code{random.cb}. Otherwise the subgroups are sampled without any random seed control.
 }

\note{
Currently only the \code{"MFI"} value for \code{startSelect} is allowed for polytomous IRT models (i.e., when \code{model} is not \code{NULL}). Otherwise an error message is returned.
}

\references{
 Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}. Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error.} 
 Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 

Magis, D. (2013). A note on the item information function of the four-parameter logistic model. \emph{Applied Psychological Measurement}, \emph{37}, 304-315. \doi{10.1177/0146621613475471}

Magis, D. and Barrada, J. R. (2017). Computerized Adaptive Testing with R: Recent Updates of the Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{Code Snippets}, \emph{76(1)}, 1-18. \doi{10.18637/jss.v076.c01}

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. \doi{10.18637/jss.v048.i08}
}

\author{
    David Magis \cr
    Department of Psychology, University of Liege, Belgium \cr
    \email{david.magis@uliege.be}

 Juan Ramon Barrada \cr
Department of Psychology and Sociology, Universidad Zaragoza, Spain \cr
\email{barrada@unizar.es}
 }

\seealso{
\code{\link{testList}}, \code{\link{genPolyMatrix}} 
} 

\examples{

## Dichotomous models ##

 # Loading the 'tcals' parameters 
 data(tcals)

 # Item bank creation with 'tcals' item parameters
 bank <- as.matrix(tcals[,1:4])
 
 # Random selection of 4 starting items
 startItems(bank, seed = 1, nrItems = 4)

 # Random selection of 4 starting items without fixing the seed
 startItems(bank, seed = NA, nrItems = 4)
 startItems(bank, seed = NA, nrItems = 4) # may provide a different result!

 ## With content balancing control
 prov <- breakBank(tcals)
 cbGroup <- prov$cbGroup

 # Creation of the 'cbList' list with arbitrary proportions
 cbList <- list(names = c("Audio1", "Audio2", "Written1", "Written2", "Written3"), 
        props = c(0.1, 0.2, 0.2, 0.2, 0.3))

 startItems(bank, seed = 1, nrItems = 3, cbControl = cbList, cbGroup = cbGroup)
 startItems(bank, seed = NA, nrItems = 3, cbControl = cbList, cbGroup = cbGroup,
            random.cb = 1)

 # Selection of the first 5 starting items
 startItems(bank, fixItems = 1:5)

 # Selecting 1 starting item, initial ability estimate is 0
 startItems(bank) 

 # Selecting 3 starting items for ability levels -1, 0 and 2
 startItems(bank, theta = c(-1, 0, 2)) 

 # Same with 5 randomesque items per theta value
 startItems(bank, theta = c(-1, 0, 2), randomesque = 5)

 # 5 randomesque items per theta value, with fixed random seed number
 startItems(bank, theta = c(-1, 0, 2), randomesque = 5, random.seed = 1) 

 # Idem but with 'bOpt' criterion
 startItems(bank, theta = c(-1, 0, 2), startSelect = "bOpt")

 # Selecting only the first 10 items as available items
 avail <- c(rep(1, 10), rep(0, nrow(bank)-10))
 startItems(bank, theta = c(-1, 0, 2), nAvailable = avail) 

\dontrun{
 # Selecting too many items among available ones
 startItems(bank, nrItems = 11, theta = 2, halfRange = 3, nAvailable = avail) 
}

## Polytomous models ##

 # Generation of an item bank under GRM with 100 items and at most 4 categories
 m.GRM <- genPolyMatrix(100, 4, "GRM")
 m.GRM <- as.matrix(m.GRM)

 # Random selection of 4 starting items
 startItems(m.GRM, model = "GRM", seed = 1, nrItems = 4)

 # Selection of the first 5 starting items
 startItems(m.GRM, model = "GRM", fixItems = 1:5)

 # Selecting 3 starting items for theta values -1, 0 and 2
 startItems(m.GRM, model = "GRM", theta = c(-1, 0, 2))

\dontrun{

 # Idem but with 'bOpt' criterion
 startItems(m.GRM, model = "GRM", nrItems = 3, theta = 1, halfRange = 2, 
            startSelect = "bOpt")
 }

 # Selecting only the first 10 items as available items
 avail <- c(rep(1, 10), rep(0, nrow(m.GRM)-10))
 startItems(m.GRM, model = "GRM", theta = c(-1, 0, 2), 
            nAvailable = avail) 

\dontrun{

 # Selecting too many items among available ones
 startItems(m.GRM, model = "GRM", theta = seq(from = -2, to = 2, length = 11), 
            nAvailable = avail) 
}

 # Loading the cat_pav data
 data(cat_pav)
 cat_pav <- as.matrix(cat_pav)

 # Random selection of 4 starting items
 startItems(cat_pav, model = "GPCM", seed = 1, nrItems = 4)

 # Selection of the first 5 starting items
 startItems(cat_pav, model = "GPCM", fixItems = 1:5)

 # Selecting 3 starting items for theta values -1, 0 and 2
 startItems(cat_pav, model = "GPCM", theta = c(-1, 0, 2)) 

\dontrun{

 # Idem but with 'bOpt' criterion
 startItems(cat_pav, model = "GPCM", theta = c(-1, 0, 2), startSelect = "bOpt")
 }

 # Selecting only the first 10 items as available items
 avail <- c(rep(1, 10), rep(0, nrow(cat_pav)-10))
 startItems(cat_pav, model = "GPCM", theta = c(-1, 0, 2), nAvailable = avail) 

\dontrun{

 # Selecting too many items among available ones
 startItems(cat_pav, model = "GPCM", theta = seq(from = -2, to = 2, length = 11), 
            nAvailable = avail) 
 }
 }

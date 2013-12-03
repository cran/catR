\name{nextItem}
\alias{nextItem}

\title{Selection of the next item}

\description{
 This command selects the next item to be administered, given the list of previously administered items and the current
 ability estimate, with several possible criteria. Item exposure and content balancing can also be controlled.
 }

\usage{
 nextItem(itemBank, theta, out=NULL, x=NULL, criterion="MFI",
 	method="BM", priorDist="norm", priorPar=c(0,1), D=1,
 	range=c(-4,4), parInt=c(-4,4,33), infoType="observed",
 	randomesque=1, cbControl=NULL) 
 }

\arguments{
 \item{itemBank}{an item bank of class \code{itBank} as output of the function \code{\link{createItemBank}}.}
 \item{theta}{numeric: the current value of the ability estimate (default is 0).}
 \item{out}{either a vector of integer values specifying the items previously administered, or \code{NULL} (default).} 
 \item{x}{numeric: the provisional response pattern, with the same length as \code{out} (and \code{NULL} by default).
          Ignored if \code{method} is either \code{"MFI"} or \code{"Urry"}. See \bold{Details}.}
 \item{criterion}{character: the method for next item selection. Possible values are \code{"MFI"} (default), \code{"Urry"}
               \code{"MLWI"}, \code{"MPWI"}, \code{"MEI"}, \code{"MEPV"} and \code{random}. See \bold{Details}.}
 \item{method}{character: the ability estimator. Possible values are \code{"BM"} (default), \code{"ML"} and \code{"WL"}. See \bold{Details}.}
 \item{priorDist}{character: the prior ability distribution. Possible values are \code{"norm"} (default) for the normal distribution, and \code{"unif"} for
                  the uniform distribution. Ignored if \code{type} is not \code{"MPWI"}.}
 \item{priorPar}{numeric: a vector of two components with the prior parameters. If \code{priorDist} is \code{"norm"}, then \code{priorPar} contains the 
                 mean and the standard deviation of the normal distribution. If \code{priorDist} is \code{"unif"}, then \code{priorPar} contains the bounds
                 of the uniform distribution. The default values are 0 and 1 respectively. Ignored if \code{type} is not \code{"MPWI"}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952).} 
 \item{range}{numeric: vector of two components specifying the range wherein the ability estimate must be looked for (default is \code{c(-4,4)}).
 	Ignored if \code{method=="EAP"}.}
 \item{parInt}{numeric: a vector of three numeric values, specifying respectively the lower bound, the upper bound and the
               number of quadrature points for numerical integration (default is \code{c(-4,4,33)}). Ignored if \code{method} is either \code{"MFI"}
               or \code{"Urry"}. See \bold{Details}.}
 \item{infoType}{character: the type of information function to be used. Possible values are \code{"observed"} (default) for observed
                 information function, and \code{"Fisher"} for Fisher information function. Ignored if \code{criterion} is not \code{"MEI"}.}
 \item{randomesque}{integer: the number of items to be chosen from the next item selection rule, among those the next item to be administered will be
                   randomly picked up. Default value is 1 and leads to usual selection of the optimal item for the specified \code{criterion}. See \bold{Details}.}
 \item{cbControl}{either a list of accurate format to control for content balancing, or \code{NULL}. See \bold{Details}.}
}

\value{
 A list with eight arguments:
 \item{item}{the selected item (identified by its number in the item bank).}
 \item{par}{the vector of item parameters of the selected item.}
 \item{info}{the value of the MFI, Fisher's information, the MLWI, the MPWI, the MEI, the EPV,or NA (for \code{"random"} criterion) for the selected 
             item and the current ability estimate.}
 \item{criterion}{the value of the \code{criterion} argument.}
 \item{randomesque}{the value of the \code{randomesque} argument.}
 \item{prior.prop}{a vector with empirical proportions of items previously administered for each subgroup of items set by \code{cbControl}.}
 \item{post.prop}{a vector with empirical proportions of items previously administered, together with the one currently selected, for each subgroup
                  of items set by \code{cbControl}.}
 \item{th.prop}{a vector with theoretical proportions given by \code{cbControl$props}.}
 }
 

\details{
 Currently seven methods are available for selecting the next item to be administered in the adaptive test. For a given
 current ability estimate, the next item is selected (among the available items) by using: the maximum Fisher information (MFI) criterion,
 the maximum likelihood weighted information (MLWI) (Veerkamp and Berger, 1997), the maximum posterior weighted information (MPWI)
 (van der Linden, 1998), Urry's procedure (Urry, 1970), the maximum expected information (MEI) criterion (van der Linden, 1998), the minimum 
 expected posterior variance (MEPV) or by selecting the next item completely randomly among the available items. 

 The MFI criterion selects the next item as the one which maximizes the item information function (Baker, 1992). The most informative item
 is selected from the table of item informations provided by the bank of items specified with \code{itemBank}. Urry's procedure consists 
 in selecting as next the item whose difficulty level is closest to the current ability estimate. Under the 1PL model, both Urry and MFI
 methods are equivalent. The MLWI and MPWI criteria select the next item as the one with maximal information, weighted either by the
 likelihood function or the posterior distribution. See the function \code{\link{MWI}} for further details. Finally, the MEI criterion 
 selects the item with maximum expected information, computed with the \code{\link{MEI}} function. 

 The method for next item selection is specified by the \code{criterion} argument. Possible values are \code{"MFI"} for maximum Fisher information
 criterion, \code{"Urry"} for Urry's method, \code{"MLWI"} for maximum likelihood weighted information criterion, \code{"MPWI"} for the
 maximum posterior weighted information criterion, \code{"MEI"} for the maximum expected information criterion, \code{"MEPV"} for minimum expected
 posterior variance, and \code{"random"} for random selection. Other values return an error message.

 For MFI, MEI and Urry criteria, the provisional ability estimate must be supplied throught the \code{theta} argument (by default, it is equal to zero).
 For MLWI and MPWI criteria, this argument is ignored. 

 The available items are those that are not specified in the \code{out} argument. By default, \code{out} is \code{NULL}, which means that all
 items are available.  

 For MEI, MEPV, MLWI and MPWI methods, the provisional response pattern must be provided through the \code{x} argument. It must be of 0/1 entries and of
 the same length as the \code{out} argument. It is ignored with MFI and Urry criteria. Moreover, the range of integration (or posterior variance 
 computation) is specified by the triplet \code{parInt}, where the first, second, and third value correspond to the arguments \code{lower}, \code{upper}
 and \code{nqp} of the \code{\link{MWI}} function, respectively. 

 The \code{method}, \code{priorDist}, \code{priorPar}, \code{D}, \code{range} and \code{intPar} arguments fix the ability estimator. See the 
 \code{\link{thetaEst}} function for further details.

 Finally, for MEI criterion, the type of information function must be supplied through the \code{infoType} argument. It is equal to 
 \code{"observed"} by default, which refers to the observed information function, and the other possible value is \code{"Fisher"}
 for Fisher information function. See the \code{\link{MEI}} funtion for further details. This argumpent is ignored if \code{criterion} is not \code{"MEI"}.

 Item exposure can be controlled by using the so-called \emph{randomesque} approach (Kingsbury an Zara, 1989), which consists in selecting more than one item
 as the best items to be administered (according to the specified \code{criterion}). The final item that is administered is randomly chosen among this set of
 optimal items.  The argument \code{randomesque} controls for the number of optimal items to be selected. The default value is 1, which corresponds to the
 usual framework of selecting the optimal item for next administration. Note that, for compatibility issues, if the number of remaining items is smaller than
 \code{randomesque}, the latter is replaced by this number of remaining items.

 Control for content balancing is also possible, given two conditions: (a) the item bank set by \code{itemBank} holds as element \code{$cbGroup} the names of the
 subgroups of items for content balancing, and (b) the argument \code{cbControl} is a correctly specified list. The correct format for \code{cbControl} is a list
 with two elements. The first one is called \code{names} and holds the names of the subgroups of items (in the order that is prespecified by the user). The second
 element is called \code{props} and contains the (theoretical) proportions of items to be administered from eazch subgroup for content balancing. These proportions
 must be strictly positive but may not sum to one; in this case they are internally normalized to sum to one. Note that if \code{cbControl} is misspecified, then
 the \code{\link{test.cbList}} will return a warning message and the \code{\link{nextItem}} function will stop. 

 Under content balancing, the selection of the next item is done in several steps.
 \enumerate{
 \item{If no item was administered yet, one subgroup is randomly picked up and the optimal item from this subgroup is selected.}
 \item{If at least one subgroup wasn't targeted yet by item selection, one of these subgroups is randomly picked up and the optimal item
       from this subgroup is selected.}
 \item{If at least one item per subgroup was already administered, the empirical relative proportions of items administered per subgroup are computed, and 
       the subgroup(s) whose difference between empirical and theoretical (i.e. given by \code{cbControl$props}) proportions is (are) selected. The optimal
       item is then selected from this subgroup for next administration (in case of several such groups, one group is randomly picked up first).}
 }
 See Kingsbury and Zara (1989) for further details.

 Three vectors of proportions are returned in the output list: \code{$prior.prop} contains the empirical relative proportions for items already administered (i.e.
 passed through the \code{out} argument); \code{$post.prop} contains the same empirical relative proportions but including the optimal item that was just selected;
 and \code{$th.prop} contains the theoretical proportions (i.e. those from \code{cbControl$props} or the normalized values).
}

\note{
 van der Linden and Pashley (2000) also introduced the Maximum Expected Posterior Weighted Information (MEPWI) criterion, as
 a mix of both MEI and MPWI methods. However, Choi and Swartz (2009) established that this method is completely equivalent to
 MPWI. For this reason, MEPWI was not implemented here.
}

\references{
 Baker, F.B. (1992). \emph{Item response theory: parameter estimation techniques}. New York, NY: Marcel Dekker.

 Choi, S. W., and Swartz, R. J. (2009). Comparison of CAT item selection criteria for polytomous items. \emph{Applied PScyhological Measurement, 32}, 419-440.

 Kingsbury, G. G., and Zara, A. R. (1989). Procedures for selecting items for computerized adaptive tests. \emph{Applied Measurement in Education, 2}, 359-375.

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R 
 Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. URL \url{http://www.jstatsoft.org/v48/i08/}

 Urry, V. W. (1970). \emph{A Monte Carlo investigation of logistic test models.} Unpublished doctoral dissertation. West Lafayette, IN: Purdue University.

 van der Linden, W. J. (1998). Bayesian item selection criteria for adaptive testing. \emph{Psychometrika, 63}, 201-216.

 van der Linden, W. J., and Pashley, P. J. (2000). Item selection and ability estimlation in adaptive testing. In W. J. van der Linden and C. A. W. Glas (Eds.), 
 \emph{Computerized adaptive testing. Theory and practice} (pp. 1-25). Boston, MA: Kluwer. 

 Veerkamp, W. J. J., and Berger, M. P. F. (1997). Some new item selection criteria for adaptive testing.
 \emph{Journal of Educational and Behavioral Statistics, 22}, 203-226.
}

\author{
    David Magis \cr
    Post-doctoral researcher, FNRS (Fonds National de la Recherche Scientifique) \cr
    Department of Mathematics, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{createItemBank}}, \code{\link{MWI}}, \code{\link{MEI}}, \code{\link{thetaEst}}, \code{\link{test.cbList}}, \code{\link{randomCAT}}
} 

\examples{
\dontrun{

 # Loading the 'tcals' parameters 
 data(tcals)
 
 # Item bank creation with 'tcals' item parameters
 bank <- createItemBank(tcals)

 ## MFI criterion

 # Selecting the next item, current ability estimate is 0
 nextItem(bank, 0) # item 63 is selected

 # Selecting the next item, current ability estimate is 0
 # and item 63 is removed
 nextItem(bank, 0, out=63) # item 10 is selected

 # Selecting the next item, current ability estimate is 0
 # and items 63 and 10 are removed
 nextItem(bank, 0, out=c(63,10)) # item 62 is selected

 # Item exposure control by selecting three items
 # (selected item will be either 10, 62 or 63)
 nextItem(bank, 0, randomesque = 3)

 ## Urry's method

 # Selecting the next item, current ability estimate is 0
 nextItem(bank, 0, criterion="Urry") # item 24 is selected

 # Selecting the next item, current ability estimate is 0
 # and item 24 is removed
 nextItem(bank, 0, out=24, criterion="Urry")
 
 ## MLWI and MPWI methods

 # Selecting the next item, current response pattern is 0
 # and item 63 was administered first
 nextItem(bank, x=0, out=63, criterion="MLWI") 
 nextItem(bank, x=0, out=63, criterion="MPWI")

 # Selecting the next item, current response pattern is
 # (0,1) and item 19 is removed
 nextItem(bank, x=c(0,1), out=c(63, 19), criterion="MLWI")
 nextItem(bank, x=c(0,1), out=c(63, 19), criterion="MPWI")

 ## MEI method

 # Selecting the next item, current response pattern is 0
 # and item 63 was administered first
 nextItem(bank, x=0, out=63, criterion="MEI") 

 # With Fisher information
 nextItem(bank, x=0, out=63, criterion="MEI", infoType="Fisher")

 ## MEPV method

 # Selecting the next item, current response pattern is 0
 # and item 63 was administered first
 nextItem(bank, x=0, out=63, criterion="MEPV") 

 ## Random method

 # Selecting the next item, item 63 was administered first
 nextItem(bank, out=63, criterion="random") 
 nextItem(bank, out=63, criterion="random")  # may produce a 
                                             # different result

 # Item bank creation for content balancing
 bank2 <- createItemBank(tcals, cb=TRUE)

 # Creation of the 'cbList' list with arbitrary proportions
 cbList <- list(names=c("Audio1","Audio2","Written1","Written2",
                "Written3"), props=c(0.1,0.2,0.2,0.2,0.3))

 # Selecting the next item, MFI criterion, current ability 
 # estimate is 0, items 12,33,46 and 63 previously administered
 nextItem(bank2, 0, out=c(12,33,46,63), cbControl=cbList) 
                                        # item 70 is selected

 }
 }

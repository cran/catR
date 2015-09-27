\name{nextItem}
\alias{nextItem}

\title{Selection of the next item}

\description{
 This command selects the next item to be administered, given the list of previously administered items and the current
 ability estimate, with several possible criteria. Item exposure and content balancing can also be controlled.
 }

\usage{
nextItem(itemBank, model = NULL, theta = 0, out = NULL, x = NULL, 
 	criterion = "MFI", method = "BM", priorDist = "norm", priorPar = c(0, 1), 
 	D = 1, range = c(-4, 4), parInt = c(-4, 4, 33), infoType = "observed", 
 	randomesque = 1, rule = "length", thr = 20, SETH = NULL, AP = 1, 
 	nAvailable = NULL, maxItems = 50, cbControl = NULL, cbGroup = NULL) 
 }

\arguments{
 \item{itemBank}{numeric: a suitable matrix of item parameters. See \bold{Details}.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
 \item{theta}{numeric: the current value of the ability estimate (default is 0). Ignored if \code{criterion} is either \code{"MLWI"}, \code{"MPWI"} or \code{"random"}. See \bold{Details}.}
 \item{out}{either a vector of integer values specifying the items previously administered, or \code{NULL} (default).} 
 \item{x}{numeric: the provisional response pattern, with the same length as \code{out} (and \code{NULL} by default). Ignored if \code{criterion} is either \code{"MFI"}, \code{"bOtp"}, \code{"thOpt"}, \code{"proportional"}, \code{"progressive"} or \code{"random"}. See \bold{Details}.}
 \item{criterion}{character: the method for next item selection. Possible values are \code{"MFI"} (default), \code{"bOpt"}, \code{"thOpt"}, \code{"MLWI"}, \code{"MPWI"}, \code{"MEI"}, \code{"MEPV"}, \code{"progressive"}, \code{"proportional"}, \code{"KL"}, \code{"KLP"}, \code{"GDI"}, \code{"GDIP"} and \code{random}. See \bold{Details}.}
 \item{method}{character: the ability estimator. Possible values are \code{"BM"} (default), \code{"ML"} and \code{"WL"}. Ignored if \code{method} is not \code{"MEI"}. See \bold{Details}.}
 \item{priorDist}{character: the prior ability distribution. Possible values are \code{"norm"} (default) for the normal distribution, and \code{"unif"} for the uniform distribution. Ignored if \code{type} is not \code{"MPWI"}, \code{"KLP"} or \code{"GDIP"}.}
 \item{priorPar}{numeric: a vector of two components with the prior parameters. If \code{priorDist} is \code{"norm"}, then \code{priorPar} contains the mean and the standard deviation of the normal distribution. If \code{priorDist} is \code{"unif"}, then \code{priorPar} contains the bounds of the uniform distribution. The default values are 0 and 1 respectively. Ignored if \code{type} is neither \code{"MPWI"} nor \code{"KLP"}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952).} 
 \item{range}{numeric: vector of two components specifying the range wherein the ability estimate must be looked for (default is \code{c(-4,4)}).}
 \item{parInt}{numeric: a vector of three numeric values, specifying respectively the lower bound, the upper bound and the number of quadrature points for numerical integration (default is \code{c(-4,4,33)}). Ignored if \code{method} is either \code{"MFI"}, \code{"bOpt"}, \code{"thOpt"}, \code{"progressive"}, \code{"proportional"} or \code{"random"}. See \bold{Details}.}
 \item{infoType}{character: the type of information function to be used. Possible values are \code{"observed"} (default) and \code{"Fisher"}. Ignored if \code{criterion} is not \code{"MEI"}. See \bold{Details}.}
 \item{randomesque}{integer: the number of items to be chosen from the next item selection rule, among those the next item to be administered will be randomly picked up. Default value is 1 and leads to usual selection of the optimal item for the specified \code{criterion}. See \bold{Details}.}
\item{rule}{character: the type of stopping rule for next item selection. Possible values are \code{"length"} (default) or \code{"precision"}. Ignored if \code{criterion} is any other value than \code{"progressive"} or \code{"proportional"}. See \bold{Details}.}
\item{thr}{numeric: the threshold related to the stopping rule. Ignored if \code{criterion} is neither \code{"progressive"} nor \code{"proportional"}. See \bold{Details}.}
\item{SETH}{either a numeric value for the provisional standard error\code{NULL} (default). Ignored if \code{criterion} is neither \code{"progressive"} nor \code{"proportional"}, or if \code{rule} is not \code{"precision"}. See \bold{Details}.}
\item{AP}{numeric: the value of the acceleration parameter  (default value is 1). Ignored if \code{criterion} is neither \code{"progressive"} nor \code{"proportional"}. See \bold{Details}.}
\item{nAvailable}{either a numeric vector of zero and one entries to denote respectively which items are not available and are availbale, or \code{NULL} (default). Used for content balancing purposes only. See \bold{Details}.}
\item{maxItems}{integer: the maximum number of items to be administered during the adaptive test (default value is 50). See \bold{Details}.}
 \item{cbControl}{either a list of accurate format to control for content balancing, or \code{NULL}. See \bold{Details}.}
 \item{cbGroup}{either a factor vector of accurate format to control for content balancing, or \code{NULL}. See \bold{Details}.}
}

\details{
 Currently twelve methods are available for selecting the next item to be administered in the adaptive test. All are avilable with dichotomous items and ten out of the twelve are also available for polytomous items. For a given current ability estimate, the next item is selected (among the available items) by using: 
\enumerate{
\item{the maximum Fisher information (MFI) criterion,}
\item{the so-called \code{bOpt} procedure (Urry, 1970) (\bold{not} for polytomous items),}
\item{the so-called \code{thOpt} procedure (see e.g., Barrada, Mazuela and Olea, 2006; Magis, 2013) (\bold{not} for polytomous items),}
\item{the maximum likelihood weighted information (MLWI) (Veerkamp and Berger, 1997),}
\item{the maximum posterior weighted information (MPWI)
 (van der Linden, 1998),}
\item{the maximum expected information (MEI) criterion (van der Linden, 1998),}
\item{the minimum expected posterior variance (MEPV),}
\item{the Kullback-Leibler (KL) divergency criterion (Chang and Ying, 1996),}
\item{the posterior Kullback-Leibler (KLP) criterion (Chang and Ying, 1996),}
\item{the progressive method (Barrada, Olea, Ponsoda, and Abad, 2008, 2010; Revuelta and Ponsoda, 1998),}
\item{the proportional method (Barrada, Olea, Ponsoda, and Abad, 2008, 2010; Segall, 2004),}
\item{the global-discrimination index (GDI) (Kaplan, de la Torre, and Barrada, 2015),}
\item{the posterior global-discrimination index (GDIP) (Kaplan, de la Torre, and Barrada, 2015),}
\item{or by selecting the next item completely \emph{randomly} among the available items.}
}
 
The MFI criterion selects the next item as the one which maximizes the item information function (Baker, 1992). The most informative item is selected from the item informations computed from the bank of items specified with \code{itemBank}. 

The so-called \code{bOpt} method (formerly referred to as \emph{Urry's procedure}) consists in selecting as next the item whose difficulty level is closest to the current ability estimate. Under the 1PL model, both bOpt and MFI methods are equivalent. This method is not available with polytomous items.

The so-called \code{thOpt} method consists in selecting the item for which optimal \eqn{\theta} value (that is, the value \eqn{\theta^*} for which item information is maximal) is closest to the current ability estimate. Under the 1PL and 2PL models, both bOpt and thOpt methods are equivalent. This method is not available with polytomous items.

The MLWI and MPWI criteria select the next item as the one with maximal information, weighted either by the likelihood function or the posterior distribution. See the function \code{\link{MWI}} for further details. 

The MEI criterion selects the item with maximum expected information, computed with the \code{\link{MEI}} function. 

The MEPV criterion selects the item with minimum expected posterior variance, computed with the \code{\link{EPV}} function.

The KL and KLP criteria select the item with maximum Kullback-Leibler (KL) information (or the posterior KL information in case of KLP criterion). This information is computed by the \code{\link{KL}} function.

With the progressive method the item selected is the one that maximizes the sum of two elements, a random part and a part determined by the Fisher information. At the beginning of the test, the importance of the random element is maximum; as the test advances, the information increases its relevance in the item selection. The speed for the transition from purely random selection to purely information based selection is determined by the acceleration parameter, set by the argument \code{AP}, where higher values imply a greater importance of the random element during the test. This method is not available with \code{rule="classification"}.

In the proportional method the items are randomly selected with probabilities of selection determined by their Fisher information  raised to a given power. This power is equal to 0 at the beginning of the test and increases as the test advances. This implies that the test starts with completely random selection and approaches the MFI at the end of the test. Here, the acceleration parameter, set by the argument \code{AP}, plays a similar role than with the progressive method. This method is not available with \code{rule="classification"}.

The GDI and GDIP criteria select the item with maximum global-discrimination index (GDI),  or the posterior GDI  in case of GDIP criterion. This index is computed by the \code{\link{GDI}} function.

The method for next item selection is specified by the \code{criterion} argument. Possible values are \code{"MFI"} for maximum Fisher information criterion, \code{"bOpt"} for bOpt's method, \code{"thOpt"} for the eponym method, \code{"MLWI"} for maximum likelihood weighted information criterion, \code{"MPWI"} for the maximum posterior weighted information criterion, \code{"MEI"} for the maximum expected information criterion, \code{"MEPV"} for minimum expected posterior variance, \code{"KL"} for Kullback-Leibler information method, \code{"KLP"} for posterior Kullback-Leibler information method, \code{"progressive"} for the progressive method, \code{"proportional"} for the proportional method, \code{"GDI"} for global-discrimination index method, \code{"KLP"} for posterior global-discrimination index method, and \code{"random"} for random selection. Other values return an error message. 

For all methods but MLWI, MPWI, GDI, GDIP and random criteria, the provisional ability estimate must be supplied throught the \code{theta} argument (by default, it is equal to zero). For MLWI, MPWI, G and random criteria, this argument is ignored. 

 The available items are those that are not specified in the \code{out} argument. By default, \code{out} is \code{NULL}, which means that all items are available. Typically \code{out} contains the item numbers that have been already administered. Alternatively one can reduce the number of available items by specifying the argument \code{nAvailable} appropriately. The latter is a vector of 0 and 1 values, where 1 corresponds to an available item and 0 to a non-availbale item. It works similarly as for the function \code{\link{startItems}}. If both \code{out} and \code{nAvailable} are supplied, then all items from both vectors are discarded for next item selection.

For MEI, MEPV, MLWI, MPWI, KL, KLP, GDI, and GDIP methods, the provisional response pattern must be provided through the \code{x} argument. It must be of 0/1 entries and of the same length as the \code{out} argument. It is ignored with MFI, bOpt, thOpt, progressive, proportional and random criteria. Moreover, the range of integration (or posterior variance computation) is specified by the triplet \code{parInt}, where the first, second, and third value correspond to the arguments \code{lower}, \code{upper} and \code{nqp} of e.g., the \code{\link{MWI}} function, respectively. 

The \code{method}, \code{priorDist}, \code{priorPar}, \code{D}, \code{range} and \code{intPar} arguments fix the choice of the ability estimator for MEI method. See the \code{\link{thetaEst}} function for further details. Note that \code{priorDist} and \code{priorPar} are also used to select the prior distribution with the MPWI and KLP methods. Finally, \code{parInt} is also used for numerical integration, among others for MLWI, MPWI, MEPV, KL and KLP methods. 

 Finally, for MEI criterion, the type of information function must be supplied through the \code{infoType} argument. It is equal to \code{"observed"} by default, which refers to the observed information function, and the other possible value is \code{"Fisher"} for Fisher information function. See the \code{\link{MEI}} funtion for further details. This argument is ignored if \code{criterion} is not \code{"MEI"}.

 The so-called \emph{randomesque} approach is used to improve item exposure control (Kingsbury and Zara, 1989), which consists in selecting more than one item as the best items to be administered (according to the specified \code{criterion}). The final item that is administered is randomly chosen among this set of optimal items.  The argument \code{randomesque} controls for the number of optimal items to be selected. The default value is 1, which corresponds to the usual framework of selecting the optimal item for next administration. Note that, for compatibility issues, if the number of remaining items is smaller than \code{randomesque}, the latter is replaced by this number of remaining items. The  \emph{randomesque} approach is not considered with the progressive or proportional item selection rules.

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{itemBank} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model, \code{"PCM"} for Partical Credit Model, \code{"GPCM"} for Generalized Partial Credit Model, \code{"RSM"} for Rating Scale Model and \code{"NRM"} for Nominal Response Model. The \code{itemBank} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.

Control for content balancing is also possible, given two conditions: (a) the \code{cbGroup} argument is a vector with the names of the subgroups of items for content balancing (one value per item), and (b) the argument \code{cbControl} is a correctly specified list. The correct format for \code{cbControl} is a list with two elements. The first one is called \code{names} and holds the names of the subgroups of items (in the order that is prespecified by the user). The second element is called \code{props} and contains the (theoretical) proportions of items to be administered from each subgroup for content balancing. These proportions must be strictly positive but may not sum to one; in this case they are internally normalized to sum to one. Note that \code{cbControl} should be tested with the \code{\link{test.cbList}} function prior to using \code{\link{nextItem}}.

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

In case of content balancing control, three vectors of proportions are returned in the output list: \code{$prior.prop} contains the empirical relative proportions for items already administered (i.e.
 passed through the \code{out} argument); \code{$post.prop} contains the same empirical relative proportions but including the optimal item that was just selected;
 and \code{$th.prop} contains the theoretical proportions (i.e. those from \code{cbControl$props} or the normalized values). Note that \code{NA} values are returned when no control for content balancing is specified.
}


\value{
 A list with eight arguments:
 \item{item}{the selected item (identified by its number in the item bank).}
 \item{par}{the vector of item parameters of the selected item.}
 \item{info}{the value of the MFI, Fisher's information, the MLWI, the MPWI, the MEI, the EPV, the unsigned distance between estimated ability and difficulty parameter (or the ability value were maximum Fisher information is provided) or NA (for \code{"random"} criterion) for the selected item and the current ability estimate.}
 \item{criterion}{the value of the \code{criterion} argument.}
 \item{randomesque}{the value of the \code{randomesque} argument.}
 \item{prior.prop}{a vector with empirical proportions of items previously administered for each subgroup of items set by \code{cbControl}.}
 \item{post.prop}{a vector with empirical proportions of items previously administered, together with the one currently selected, for each subgroup
                  of items set by \code{cbControl}.}
 \item{th.prop}{a vector with theoretical proportions given by \code{cbControl$props}.}
 }
 

\note{
 van der linden (1998) also introduced the Maximum Expected Posterior Weighted Information (MEPWI) criterion, as
 a mix of both MEI and MPWI methods (see also van der Linden and Pashley, 2000). However, Choi and Swartz (2009) established that this method is completely equivalent to
 MPWI. For this reason, MEPWI was not implemented here. 
}

\references{
 Baker, F.B. (1992). \emph{Item response theory: parameter estimation techniques}. New York, NY: Marcel Dekker.

Barrada, J. R., Mazuela, P., and Olea, J. (2006). Maximum information stratification method for controlling item exposure in computerized adaptive testing. \emph{Psicothema, 18}, 156-159. 

Barrada, J. R., Olea, J., Ponsoda, V., and Abad, F. J. (2008). Incorporating randomness to the Fisher information for improving item exposure control in CATS. \emph{British Journal of Mathematical and Statistical Psychology, 61}, 493-513. doi: 10.1348/000711007X230937

 Barrada, J. R., Olea, J., Ponsoda, V., and Abad, F. J. (2010). A method for the comparison of item selection rules in computerized adaptive testing. \emph{Applied Psychological Measurement}, \emph{34}, 438-452. doi: 10.1177/0146621610370152

 Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}.
 Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

Chang, H.-H., and Ying, Z. (1996). A global information approach to computerized adaptive testing. \emph{Applied Psychological Measurement, 20}, 213-229.

 Choi, S. W., and Swartz, R. J. (2009). Comparison of CAT item selection criteria for polytomous items. \emph{Applied Psychological Measurement, 32}, 419-440. doi: 10.1177/0146621608327801

 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error}. Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 

 Kaplan, M., de la Torre, J., and Barrada, J. R. (2015). New item selection methods for cognitive diagnosis computerized adaptive testing. \emph{Applied Psychological Measurement, 39}, 167-188. doi: 10.1177/0146621614554650

 Kingsbury, G. G., and Zara, A. R. (1989). Procedures for selecting items for computerized adaptive tests. \emph{Applied Measurement in Education, 2}, 359-375. doi: 10.1207/s15324818ame0204_6

Leroux, A. J., Lopez, M., Hembry, I. and Dodd, B. G. (2013). A comparison of exposure control procedures in CATs using the 3PL model. \emph{Educational and Psychological Measurement, 73}, 857-874. doi: 10.1177/0013164413486802

Magis, D. (2013). A note on the item information function of the four-parameter logistic model. \emph{Applied Psychological Measurement, 37}, 304-315. doi: 10.1177/0146621613475471

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R 
 Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. URL \url{http://www.jstatsoft.org/v48/i08/}

Revuelta, J., and Ponsoda, V. (1998). A comparison of item exposure control methods in computerized adaptive testing. \emph{Journal of Educational Measurement, 35}, 311-327. doi: 10.1111/j.1745-3984.1998.tb00541.x

Segall, D. O. (2004). A sharing item response theory model for computerized adaptive testing. \emph{Journal of Educational and Behavioral Statistics, 29}, 439-460. doi: 10.3102/10769986029004439

 Urry, V. W. (1970). \emph{A Monte Carlo investigation of logistic test models.} Unpublished doctoral dissertation. West Lafayette, IN: Purdue University.

 van der Linden, W. J. (1998). Bayesian item selection criteria for adaptive testing. \emph{Psychometrika, 63}, 201-216. doi: 10.1007/BF02294775

 van der Linden, W. J., and Pashley, P. J. (2000). Item selection and ability estimation in adaptive testing. In W. J. van der Linden and C. A. W. Glas (Eds.), \emph{Computerized adaptive testing. Theory and practice} (pp. 1-25). Boston, MA: Kluwer. 

 Veerkamp, W. J. J., and Berger, M. P. F. (1997). Some new item selection criteria for adaptive testing. \emph{Journal of Educational and Behavioral Statistics, 22}, 203-226. doi: 10.3102/10769986022002203
}

\author{
    David Magis \cr
    Department of Education, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}

 Juan Ramon Barrada \cr
Department of Psychology and Sociology, Universidad Zaragoza, Spain \cr
\email{barrada@unizar.es}

Guido Corradi \cr
Department of Psychology and Sociology, Universidad Zaragoza, Spain \cr
\email{guidocor@gmail.com}
 }

\seealso{
\code{\link{MWI}}, \code{\link{MEI}}, \code{\link{KL}}, \code{\link{thetaEst}}, \code{\link{test.cbList}}, \code{\link{randomCAT}}, \code{\link{genPolyMatrix}} 
} 

\examples{

## Dichotomous models ##

 # Loading the 'tcals' parameters 
 data(tcals)
 
 # Item bank creation with 'tcals' item parameters
 prov <- breakBank(tcals)
 bank <- prov$itemPar
 cbGroup <- prov$cbGroup

 ## MFI criterion

 # Selecting the next item, current ability estimate is 0
 nextItem(bank, theta = 0) # item 63 is selected

 # Selecting the next item, current ability estimate is 0 and item 63 is removed
 nextItem(bank, theta = 0, out = 63) # item 10 is selected

 # Selecting the next item, current ability estimate is 0 and items 63 and 10 are 
 # removed
 nextItem(bank, theta = 0, out = c(63, 10)) # item 62 is selected

 # Item exposure control by selecting three items (selected item will be either 10, 62 
 # or 63)
 nextItem(bank, theta = 0, randomesque = 3)


 ## bOpt method

 # Selecting the next item, current ability estimate is 0
 nextItem(bank, theta = 0, criterion = "bOpt") # item 24 is selected

 # Selecting the next item, current ability estimate is 0 and item 24 is removed
 nextItem(bank, theta = 0, out = 24, criterion = "bOpt")


 ## thOpt method

 # Selecting the next item, current ability estimate is 0
 nextItem(bank, theta = 0, criterion = "thOpt") # item 76 is selected

 # Selecting the next item, current ability estimate is 0 and item 76 is removed
 nextItem(bank, theta = 0, out = 76, criterion = "thOpt") # item 70 is selected


 ## MLWI and MPWI methods

 # Selecting the next item, current response pattern is 0 and item 63 was administered 
 # first
 nextItem(bank, x = 0, out = 63, criterion = "MLWI") 
 nextItem(bank, x = 0, out = 63, criterion = "MPWI")

 # Selecting the next item, current response pattern is (0,1) and item 19 is removed
 nextItem(bank, x = c(0, 1), out = c(63, 19), criterion = "MLWI")
 nextItem(bank, x = c(0, 1), out = c(63, 19), criterion = "MPWI")


 ## MEI method

 # Selecting the next item, current response pattern is 0 and item 63 was administered 
 # first
 # Ability estimation by WL method
 th <- thetaEst(rbind(bank[63,]), 0, method = "WL")
 nextItem(bank, x = 0, out = 63, theta = th, criterion = "MEI") # item 49 is selected

 # With Fisher information
 nextItem(bank, x = 0, out = 63, theta = th, criterion = "MEI", infoType = "Fisher") 
   # item 10 is selected

\dontrun{

 ## MEPV method

 # Selecting the next item, current response pattern is 0 and item 63 was administered 
 # first
 # Ability estimation by WL method
 nextItem(bank, x = 0, out = 63, theta = th, criterion = "MEPV") # item 19 is selected

}

 ## KL and KLP methods

 # Selecting the next item, current response pattern is 0 and item 63 was administered 
 # first
 # Ability estimation by WL method
 nextItem(bank, x = 0, out = 63, theta = th, criterion = "KL") # item 19 is selected
 nextItem(bank, x = 0, out = 63, theta = th, criterion = "KLP") # item 44 is selected


## GDI and GDIP methods

 # Selecting the next item, current response pattern is 0 and item 63 was administered 
 # first
 nextItem(bank, x = 0, out = 63, criterion = "GDI") # item 49 is selected
 nextItem(bank, x = 0, out = 63, criterion = "GDIP") # item 44 is selected


 ## Progressive method

 # Selecting the next item, current ability estimate is 0 and item 63 was administered 
 # first
 # (default options: "length" rule with "thr = 20")
 nextItem(bank, out = 63, theta = 0, criterion = "progressive") 
 nextItem(bank, out = 63, theta = 0, criterion = "progressive")
  # result can be different!


 ## Proportional method

 # Selecting the next item, current ability estimate is 0 and item 63 was administered 
 # first
 # (default options: "length" rule with "thr = 20")
 nextItem(bank, out = 63, theta = 0, criterion = "proportional") 
 nextItem(bank, out = 63, theta = 0, criterion = "proportional")
  # result can be different!


 ## Random method

 # Selecting the next item, item 63 was administered first
 nextItem(bank, out = 63, criterion = "random") 
 nextItem(bank, out = 63, criterion = "random")  # may produce a different result


 ## Content balancing

 # Creation of the 'cbList' list with arbitrary proportions
 cbList <- list(names = c("Audio1", "Audio2", "Written1", "Written2", "Written3"), 
        props = c(0.1, 0.2, 0.2, 0.2, 0.3))

 # Selecting the next item, MFI criterion, current ability estimate is 0, items 12, 33, 
 # 46 and 63 previously administered
 nextItem(bank, theta = 0, out = c(12, 33, 46, 63), cbControl = cbList, 
          cbGroup = cbGroup)  # item 70 is selected


## Polytomous models ##

 # Generation of an item bank under GRM with 100 items and at most 4 categories
 m.GRM <- genPolyMatrix(100, 4, "GRM")
 m.GRM <- as.matrix(m.GRM)

 # Current ability estimate is 0
 # Selecting the next item, current response pattern is 1 and item 70 was administered 
 # first

 ## MFI method
 nextItem(m.GRM, model = "GRM", theta = 0, criterion = "MFI", out = 70)

 ## Progressive method
 nextItem(m.GRM, model = "GRM", theta = 0, criterion = "progressive", out = 70)

 ## KL method
 nextItem(m.GRM, model = "GRM", theta = 0, criterion = "KL", out = 70, x = 1)

 ## MFI with content balancing
 cbList <- list(names = c("Audio1","Audio2","Written1","Written2", "Written3"), 
        props = c(0.1,0.2,0.2,0.2,0.3))
 m.GRM<-genPolyMatrix(100, 4, model = "GRM", cbControl = cbList)
 bank<-breakBank(m.GRM)
 nextItem(bank$itemPar, model = "GRM", theta = 0, criterion = "MFI", out = 70,
  cbControl = cbList, cbGroup = bank$cbGroup)


 # Generation of an item bank under PCM with 20 items and 4 categories
 m.PCM <- genPolyMatrix(20, 4, "PCM", same.nrCat = TRUE)
 m.PCM <- as.matrix(m.PCM)

 # Current ability estimate is 0
 # Selecting the next item, current response pattern is 1 and item 15 was administered 
 # first

 ## MFI method
 nextItem(m.PCM, model = "PCM", theta = 0, criterion = "MFI", out = 15)

 ## Progressive method
 nextItem(m.PCM, model = "PCM", theta = 0, criterion = "progressive", out = 15)

 ## KL method
 nextItem(m.PCM, model = "PCM", theta = 0, criterion = "KL", out = 15, x = 1)

 ## MFI with content balancing
 cbList <- list(names = c("Audio1", "Audio2", "Written1", "Written2", "Written3"), 
        props = c(0.1, 0.2, 0.2, 0.2, 0.3))
 m.PCM<-genPolyMatrix(100, 4, model = "PCM", cbControl = cbList)
 bank<-breakBank(m.PCM)
 nextItem(bank$itemPar, model = "PCM", theta = 0, criterion = "MFI", out = 15, 
          cbControl = cbList, cbGroup = bank$cbGroup)
 }

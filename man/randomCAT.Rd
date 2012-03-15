\name{randomCAT}
\alias{randomCAT}
\alias{print.cat}
\alias{plot.cat}

\title{Random generation of adaptive tests}

\description{
 This command generates a response pattern to an adaptive test, for a given item bank, a true ability level, 
 and several lists of CAT parameters (starting items, stopping rule, provisional and final ability estimators).
 }

\usage{
randomCAT(trueTheta, itemBank, maxItems=50, cbControl=NULL,
 	start=list(fixItems=NULL, seed=NULL, nrItems=1, theta=0, 
 	halfRange=2, startSelect="bOpt"), test=list(method="BM", 
 	priorDist="norm", priorPar=c(0,1), range=c(-4,4), D=1, 
 	parInt=c(-4,4,33), itemSelect="MFI", infoType="observed",
 	randomesque=1), stop=list(rule="length", thr=20,
 	alpha=0.05), final=list(method="BM", priorDist="norm",
 	priorPar=c(0,1), range=c(-4,4),D=1, parInt=c(-4,4,33), 
 	alpha=0.05), save.output=FALSE, output=c("out","default"))
\method{print}{cat}(x, ...)
\method{plot}{cat}(x, ci=FALSE, alpha=0.05, trueTh=TRUE, classThr=NULL, 
 	save.plot=FALSE, save.options=c("plot","default","pdf"),...)
 }

\arguments{
 \item{trueTheta}{numeric: the value of the true ability level.}
 \item{itemBank}{an item bank, i.e. a list of class \code{itBank}, typically an output of the function \code{createItemBank}.}
 \item{maxItems}{numeric: the maximal number of items to be administered (default is 50).}
 \item{cbControl}{either a list of accurate format to control for content balancing, or \code{NULL}. See \bold{Details}.}
 \item{start}{a list with the options for starting the adaptive test. See \bold{Details}.}
 \item{test}{a list with the options for provisional ability estimation and next item selection. See \bold{Details}.}
 \item{stop}{a list with the options of the stopping rule. See \bold{Details}.}
 \item{final}{a list with the options for final ability estimation. See \bold{Details}.}
 \item{save.output}{logical: should the output be saved in an external text file? (default is \code{FALSE}).}
 \item{output}{character: a vector of two components. The first component is the name of the output file, the second component is either
              the file path or \code{"default"} (default). See \bold{Details}.}
 \item{x}{an object of class "cat", typically an output of \code{randomCAT} function.}
 \item{ci}{logical: should the confidence intervals be plotted for each provisional ability estimate? (default is \code{TRUE}).}
 \item{alpha}{numeric: the significance level for provisional confidence intervals (default is 0.05). Ignored if \code{ci} is \code{FALSE}.}
 \item{trueTh}{logical: should the true ability level be drawn by a horizontal line? (default is TRUE).}
 \item{classThr}{either a numeric value giving the classification threshold to be displayed, or \code{NULL}.}
 \item{save.plot}{logical: should the plot be saved in an external figure? (default is \code{FALSE}).}
 \item{save.options}{character: a vector of three components. The first component is the name of the output file, the second component
                    is either the file path or \code{"default"} (default), and the third component is the file extension, either \code{"pdf"} (default)
                    or \code{"jpeg"}. Ignored if \code{save.plot} is \code{FALSE}. See \bold{Details}.}
 \item{...}{other generic arguments to be passed to \code{print} and \code{plot} functions.}
}

\value{
 The function \code{randomCAT} returns a list of class "cat" with the following arguments:
 \item{trueTheta}{the value of the \code{trueTheta} argument.}
 \item{maxItems}{the value of the \code{maxItems} argument.}
 \item{testItems}{a vector with the items that were administered during the test.}
 \item{itemPar}{a matrix with the parameters of the items administered during the test.}
 \item{pattern}{the generated response pattern (as vector of 0 and 1 entries).}
 \item{thetaProv}{a vector with the provisional ability estimates.}
 \item{seprov}{a vector with the standard errors of the provisional ability estimates.}
 \item{thFinal}{the final ability estimate.}
 \item{seFinal}{the standrad error of the final ability estimate.}
 \item{ciFinal}{the confidence interval of the final ability estimate.}
 \item{startFixItems}{the value of the \code{start$fixItems} argument (or its default value if missing).}
 \item{startSeed}{the value of the \code{start$seed} argument (or its default value if missing).}
 \item{startNrItems}{the value of the \code{start$nrItems} argument (or its default value if missing).}
 \item{startTheta}{the value of the \code{start$theta} argument (or its default value if missing).}
 \item{startHalfRange}{the value of the \code{start$halfRange} argument (or its default value if missing).}
 \item{startThStart}{the starting ability values used for selecting the first items of the test.}
 \item{startSelect}{the value of the \code{start$startSelect} argument (or its default value if missing).}
 \item{provMethod}{the value of the \code{test$method} argument (or its default value if missing).}
 \item{provDist}{the value of the \code{test$priorDist} argument (or its default value if missing).}
 \item{provPar}{the value of the \code{test$priorPar} argument (or its default value if missing).}
 \item{provRange}{the value of the \code{test$range} argument (or its default value if missing).}
 \item{provD}{the value of the \code{test$D} argument (or its default value if missing).}
 \item{itemSelect}{the value of the \code{test$itemSelect} argument (or its default value if missing).}
 \item{infoType}{the value of the \code{test$infoType} argument (or its default value if missing).}
 \item{randomesque}{the value of the \code{test$randomesque} argument (or its default value if missing).}
 \item{cbControl}{the value of the \code{cbControl} argument (or its default value if missing).}
 \item{cbGroup}{the value of the \code{itemBank$cbGroup} element of the item bank \code{itemBank}.}
 \item{stopRule}{the value of the \code{stop$rule} argument (or its default value if missing).}
 \item{stopThr}{the value of the \code{stop$thr} argument (or its default value if missing).}
 \item{stopAlpha}{the value of the \code{stop$alpha} argument (or its default value if missing).}
 \item{endWarning}{a logical indactor indicating whether the adaptive test stopped because the stopping rule was satisfied or not.}
 \item{finalMethod}{the value of the \code{final$method} argument (or its default value if missing).}
 \item{finalDist}{the value of the \code{final$priorDist} argument (or its default value if missing).}
 \item{finalPar}{the value of the \code{final$priorPar} argument (or its default value if missing).}
 \item{finalRange}{the value of the \code{final$range} argument (or its default value if missing).}
 \item{finalD}{the value of the \code{final$D} argument (or its default value if missing).}
 \item{finalAlpha}{the value of the \code{final$alpha} argument (or its default value if missing).}
 \item{save.output}{the value of the \code{save.output} argument.}
 \item{output}{the value of the \code{output} argument.}
 
The function \code{print.cat} returns similar (but differently organized) results. 
 }


\details{
 The \code{randomCAT} function generates an adaptive test using an item bank specified by argument \code{itemBank}, and for
 a given true ability level specified by argument \code{trueTheta}. The maximal length of the test can be fixed through
 the \code{maxItems} argument, with a default value of 50 items. 

 Content balancing can be controlled by the \code{cbControl} argument. See the \code{\link{nextItem}} function for further
 details on how to specify \code{cbControl} properly and under which conditions it is operational. By default, content balancing
 is not controlled (see Kingsbury and Zara, 1989, for further details on content balancing).

 The test specification is made by means of four lists of options: one list for the selection of the starting items, 
 one list with the options for provisional ability estimation, one list to define the stopping rule, and one list with
 the options for final ability estimation. These lists are specified respectively by the arguments \code{start}, \code{test},
 \code{stop} and \code{final}. 

 The \code{start} list can contain one or several of the following arguments:
 \itemize{
 \item{\code{fixItems}: either a vector of integer values, setting the items to be administered as first items, or 
	\code{NULL} (default) to let the function select the items.}
 \item{\code{seed}: either a numeric value to fix the random seed for item selection, or \code{NULL} (default) to
	select the items on the basis of their difficulty level. Ignored if \code{fixItems} is not \code{NULL}.}
 \item{\code{nrItems}: the number of first items to be selected (default is 1). Ignored if \code{fixItems} is not \code{NULL}.}
 \item{\code{theta}: the central initial ability value, used to define the range of ability levels for selecting the first items 
	(default is 0). Ignored if either \code{fixItems} or \code{seed} is not \code{NULL}. See \code{\link{startItems}} for further details.}
 \item{\code{halfRange}: the half range of starting ability levels for selecting the first items (default is 2). Ignored if
       either \code{fixItems} or \code{seed} is not \code{NULL}. See \code{\link{startItems}} for further details.}
 \item{\code{startSelect}: the method for selecting the first items of the test, with possible values \code{"bOpt"} (default)
       and \code{"MFI"}. Ignored if either \code{fixItems} or \code{seed} is not \code{NULL}. See \code{\link{startItems}} for further details.}
 }
 These arguments are passed to the function \code{\link{startItems}} to select the first items of the test.

 The \code{test} list can contain one or several of the following arguments:
 \itemize{
 \item{\code{method}: a character string to specify the method for ability estimation. Possible values are: \code{"BM"} (default) 
	for Bayesian modal estimation (Birnbaum, 1969), \code{"ML"} for maximum likelihood estimation (lord, 1980), \code{"EAP"} for
 	expected a posteriori (EAP) estimation (Bock and Mislevy, 1982), and \code{"WL"} for weighted likelihood estimation (Warm, 1989).}
 \item{\code{priorDist}: a character string which sets the prior distribution. Possible values are: \code{"norm"} (default)
	for normal distribution, \code{"unif"} for uniform distribution, and \code{"Jeffreys"} for Jeffreys'
	noninformative prior distribution (Jeffreys, 1939, 1946). ignored if \code{method} is neither \code{"BM"} nor \code{"EAP"}.}
 \item{\code{priorPar}: a vector of two numeric components, which sets the parameters of the prior distribution. If 
	(\code{method="BM"} or \code{method=="EAP"}) and \code{priorDist="norm"}, the components of \code{priorPar} are respectively
	the mean and the standard deviation of the prior normal density. If (\code{method="BM"} or \code{method=="EAP"}) and
	\code{priorDist="unif"}, the components of \code{priorPar} are respectively the lower and upper bound of the prior uniform density.
	Ignored in all other cases. By default, \code{priorPar} takes the parameters of the prior standard normal distribution (i.e., 
	\code{priorPar=c(0,1)}). In addition, \code{priorPar} also provides the prior parameters for the comoutation of MLWI and MPWI
      values for next item selection (see \code{\link{nextItem}} for further details).}
 \item{\code{range}: the maximal range of ability levels, set as a vector of two numeric components. The ability estimate
	will always lie to this interval (set by default to [-4, 4]). Ignored if \code{method=="EAP"}.}
 \item{\code{D}: the value of the metric constant. Default is \code{D=1} for logistic metric. Setting \code{D=1.702}
	yields approximately the normal metric (Haley, 1952).}
 \item{\code{parInt}: a numeric vector of three components, holding respectively the values of the arguments \code{lower}, 
               \code{upper} and \code{nqp} of the \code{\link{eapEst}}, \code{\link{eapSem}} and \code{MWI} commands. It specifies the 
               range of quadrature points for numerical integration, and is used for computing the EAP estimate, its standard error,
               and the MLWI and MPWI values for next item selection. Default vector is (-4, 4, 33), thus setting the range from -4 to 4 by
               steps of 0.25. Ignored if \code{method} is not \code{"EAP"} and if \code{itemSelect} is neither \code{"MLWI"} nor \code{"MPWI"}.}
 \item{\code{itemSelect}: the rule for next item selecion, with possible values \code{"MFI"} (default) for maximum Fisher information criterion;
                   \code{"Urry"} for Urry's procedure; \code{"MLWI"} and \code{"MPWI"} for respectively maximum likelihood and posterior
                   weighted information criterion; \code{"MEPV"} for minimum expected posterior variance; \code{"MEI"} for maximum expected
                   information; and \code{"random"} for random selection. For further details, see \code{\link{nextItem}}.}
 \item{\code{infoType}: character: the type of information function to be used for next item selection. Possible values are \code{"observed"}
                 (default) for observed information function, and \code{"Fisher"} for Fisher information function. Ignored if \code{itemselect}
                 is not \code{"MEI"}.}
 \item{\code{randomesque}: integer: the number of items to be chosen from the next item selection rule, among those the next item to be administered will be
                   randomly picked up. Default value is 1 and leads to usual selection of the optimal item (Kingsbury and Zara, 1989).}
 }
 These arguments are passed to the functions \code{\link{thetaEst}} (or \code{\link{eapEst}}) and \code{\link{semTheta}} (or \code{\link{eapSem}})
 to estimate the ability level and the standard error of this estimate. In addition, some arguments are passed to \code{\link{nextItem}} to select the 
 next item appropriately.

 The \code{stop} list can contain one or several of the following arguments:
 \itemize{
 \item{\code{rule}: a character string specifying the type of stopping rule. Possible values are: \code{"length"} (default),
	to stop the test after a pre-specified number of items administered; \code{"precision"}, to stop the test 
	when the provisional standard error of ability becomes less than or equal to the pre-specified value; and
	\code{"classification"}, for which the test ends whenever the provisional confidence interval (set by the 
	\code{alpha} argument) does not hold the classification threshold anymore.}
 \item{\code{thr}: a numeric value fixing the threshold of the stopping rule. If \code{rule="length"}, \code{thr} is the
	maximal number of items to be administered (in practice, it is replaced by the value of the \code{maxItems} argument
	if the latter is smaller than \code{thr}). If \code{rule="precision"}, \code{thr} is the precision level (i.e. the
	standard error) to be reached before stopping. Finally, if \code{rule="classification"}, \code{thr} corresponds to
	the ability level which serves as a classification rule (i.e. which must not be covered by the provisional 
	confidence interval).}
 \item{\code{alpha}: the significance (or \eqn{\alpha}) level for computing the priovisional confidence interval of ability.
	Ignored if \code{rule} is not \code{"classification"}.}
 }

 Eventually, the \code{final} list can contain one or several arguments of the \code{test} list (with possiblly different values),
 as well as the additional \code{alpha} argument. The latter specifies the \eqn{\alpha} level of the final confidence interval
 of ability, which is computed as
 \deqn{[\hat{\theta}-z_{1-\alpha/2} \; se(\hat{\theta}) ; \hat{\theta}+z_{1-\alpha/2} \; se(\hat{\theta})]}
 where \eqn{\hat{\theta}} and \eqn{se(\hat{\theta})} are respectively the ability estimate and its standard error. Note that the
 argument \code{itemSelect} of the \code{test} list is not used for final estimation of the ability level, and is therefore 
 not allowed into the \code{final} list.

 If some arguments of these lists are missing, they are automatically set to their default value. The contents
 of the lists is checked with the \code{\link{testList}} function, and the adaptive test is generated only if the lists
 are adequately defined. Othgerwise, a message error is printed.

 The output of \code{randomCAT}, as displayed by the \code{print.cat} function, can be stored in a text file provided that
 \code{save.output} is set to \code{TRUE} (the default value \code{FALSE} does not execute the storage). In this case, the name of
 the text file must be given as a character string into the first component of the output argument (default name is \code{"out"}),
 and the path for saving the text file can be given through the second component of output. The default value is \code{"default"},
 meaning that the file will be saved in the current working directory. Any other path can be specified as a character string: see
 the \bold{Examples} section for an illustration. 

 The function \code{plot.cat} represents the set of provisional and final ability estimates througghout the test. Corresponding
 confidence intervals (with confidence level defined by the argument \code{alpha}) are also drawn if \code{ci=TRUE} (which
 is not the default value). The true ability level can be drawn by a horizontal solid line by specifying \code{trueTh=TRUE} (which is 
 the default value); setting it to \code{FALSE} will undo the drawing. Finally, any classification threshold can be additionally 
 displayed by specifying a numeric value to the argument \code{classThr}. The default value \code{NULL} does not display any threshold.

 Finally, the plot can be saved in an external file, either as PDF or JPEG format. First, the argument \code{save.plot} must be set to
 \code{TRUE} (default is \code{FALSE}). Then, the name of the figure, its location and format are specified through the argument 
 \code{save.options}, all as character strings. See the \code{Examples} section for further information and a practical example.
}

\references{
 Birnbaum, A. (1969). Statistical theory for logistic mental test models with a prior distribution of ability.
 \emph{Journal of Mathematical Psychology, 6}, 258-276.

 Bock, R. D., and Mislevy, R. J. (1982). Adaptive EAP estimation of ability in a microcomputer environment. 
 \emph{Applied Psychological Measurement, 6}, 431-444.

 Haley, D.C. (1952). Estimation of the dosage mortality relationship when the dose is subject to error. 
 Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 
 
 Jeffreys, H. (1939). \emph{Theory of probability}. Oxford, UK: Oxford University Press.

 Jeffreys, H. (1946). An invariant form for the prior probability in estimation problems. \emph{Proceedings of the Royal
 Society of London. Series A, Mathematical and Physical Sciences, 186}, 453-461.

 Kingsbury, G. G., and Zara, A. R. (1989). Procedures for selecting items for computerized adaptive tests. \emph{Applied Measurement in Education, 2}, 359-375.

 Lord, F.M. (1980). \emph{Applications of item response theory to practical testing problems}.
 Hillsdale, NJ: Lawrence Erlbaum.

 Magis, D., and Raiche, G. (in press). Random generation of response patterns under computerized adaptive testing
 with the R package \emph{catR}. \emph{Journal of Statistical Software}.

 Urry, V. W. (1970). \emph{A Monte Carlo investigation of logistic test models.} Unpublished doctoral dissertation. West Lafayette, IN: Purdue University.

 van der Linden, W. J. (1998). Bayesian item selection criteria for adaptive testing. \emph{Psychometrika, 63}, 201-216.

 Veerkamp, W. J. J., and Berger, M. P. F. (1997). Some new item selection criteria for adaptive testing.
 \emph{Journal of Educational and Behavioral Statistics, 22}, 203-226.

 Warm, T.A. (1989). Weighted likelihood estimation of ability in item response models. \emph{Psychometrika, 54}, 427-450.
}


\author{
    David Magis \cr
    Post-doctoral researcher, FNRS (Fonds National de la Recherche Scientifique) \cr
    Department of Mathematics, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{testList}}, \code{\link{startItems}}, \code{\link{nextItem}}, \code{\link{thetaEst}}, \code{\link{semTheta}}, \code{\link{eapEst}},
\code{\link{eapSem}}, \code{\link{MWI}}, \code{\link{MEI}}
} 

\examples{
\dontrun{
 # Loading the 'tcals' parameters 
 data(tcals)
 
 # Item bank creation with 'tcals' item parameters
 bank <- createItemBank(tcals)

 # Creation of a starting list: 5 items, initial theta 0, bw 2
 start <- list(nrItems=5, theta=0, halfRange=2)

 # Creation of 'test' list: weighted likelihood
 # estimation of provisional ability, and MEI criterion 
 # for next item selection
 test <- list(method="WL", itemSelect="MEI")

 # Creation of 'final' list: EAP estimation of final
 # ability
 final <- list(method="EAP")

 # Creation of a stopping rule: precision criterion, standard
 # error to be reached 0.3
 stop <- list(rule="precision", thr=0.3)

 # CAT test
 res <- randomCAT(0, bank, start=start, test=test, stop=stop,
 	final=final)

 # New 'test' and 'final' rules (BM and EAP estimation
 # with Jeffreys' prior, randomesque value 5)
 test2 <- list(method="BM", priorDist="Jeffreys", randomesque=5)
 final2 <- list(method="EAP", priorDist="Jeffreys")

 # New stopping rule: classification criterion, with 
 # classification threshold 0 and alpha level 0.05
 stop2 <- list(rule="classification", thr=0, alpha=0.05)

 # CAT test with new 'test', 'stop' and 'final' rules
 res2 <- randomCAT(0, bank, start=start, test=test2, stop=stop2,
 	final=final2)

 # New stopping rule: classification criterion, with 
 # classification threshold 0.5 and alpha level 0.05
 stop3 <- list(rule="classification", thr=0.5, alpha=0.05)

 # CAT test with new 'stop' rule
 res3 <- randomCAT(0, bank, start=start, test=test2, stop=stop3,
 	final=final2)

 # new 'test' and 'stop' rule for next item selection
 test3 <- list(method="WL", itemSelect="MLWI")
 stop4 <- list(rule="length",thr=10)
 res4 <- randomCAT(0, bank, start=start, test=test3, stop=stop4,
 	final=final2)

 # Item bank creation for content balancing
 bank2 <- createItemBank(tcals, cb=TRUE)

 # Creation of the 'cbList' list with arbitrary proportions
 cbList <- list(names=c("Audio1","Audio2","Written1","Written2",
                "Written3"), props=c(0.1,0.2,0.2,0.2,0.3))

 # CAT test with 'start', 'test2', 'stop4' and 'final2' lists
 # and content balancing using 'cbList'
 res5 <- randomCAT(0, bank2, start=start, test=test2, stop=stop4,
 	final=final2, cbControl=cbList)

 # Saving the output in the external 'ou' text file within folder
 # 'Program Files' of hard drive 'C'
 res5 <- randomCAT(0, bank2, start=start, test=test2, stop=stop4,
 	final=final2, cbControl=cbList, save.output = TRUE,
      output=c("out","c:/Program Files/"))

 # Plotting results
 plot(res)
 plot(res, ci=TRUE)
 plot(res, ci=TRUE, trueTh=FALSE)
 plot(res, ci=TRUE, classThr=1)

 # Saving last figure into PDF file 'figure' within folder
 # 'C:/Program Files/'
 plot(res, ci=TRUE, classThr=1, save.plot=TRUE,
      save.options=c("figure","c:/Program Files","pdf"))

 # With mistake
 plot(res, ci=0.05)
 plot(res, classThr=TRUE)
 }
 }

\name{randomCAT}
\alias{randomCAT}
\alias{print.cat}
\alias{plot.cat}

\title{Random generation of adaptive tests}

\description{
 This command generates a response pattern to an adaptive test, for agiven item bank, a true ability level, 
 and several lists of CAT parameters (starting items, stopping rule, provisional and final ability estimators).
 }

\usage{
 randomCAT(trueTheta, itemBank, maxItems=50, 
 	start=list(fixItems=NULL, seed=NULL, nrItems=1, theta=0, 
 	bw=4, range=c(-4,4)), test=list(method="BM", 
 	priorDist="norm", priorPar=c(0,1), range=c(-4,4), D=1, 
 	eapPar=c(-4,4,20), itemSelect="info"), 
	stop=list(rule="length", thr=20, alpha=0.05), 
	final=list(method="BM", priorDist="norm",
	priorPar=c(0,1), range=c(-4,4),D=1, eapPar=c(-4,4,20), 
	alpha=0.05))
 \method{print}{cat}(x, ...)
 \method{plot}{cat}(x, ci=TRUE, alpha=0.05, thr=NULL, ...)
 }

\arguments{
 \item{trueTheta}{numeric: the value of the true ability level.}
 \item{itemBank}{an item bank, i.e. a list of class \code{itBank}, typically an output of the function \code{createItemBank}.}
 \item{maxItems}{numeric: the maximal number of items to be administered (default is 50).}
 \item{start}{a list with the options for starting the adaptive test. See \bold{Details}.}
 \item{test}{a list with the options for provisional ability estimation and next item selection. See \bold{Details}.}
 \item{stop}{a list with the options of the stopping rule. See \bold{Details}.}
 \item{final}{a list with the options for final ability estimation. See \bold{Details}.}
 \item{x}{an object of class "cat", typically an output of \code{randomCAT} function.}
 \item{ci}{logical: should the confidence intervals be plotted for each provisional ability estimate? (default is \code{TRUE}).}
 \item{alpha}{numeric: the significance level for provisional confidence intervals (default is 0.05). Ignored if \code{ci} is \code{FALSE}.}
 \item{thr}{either a numeric value specifying the threshold to be plotted, or \code{NULL} (default) for not displaying any threshold.}
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
 \item{startBw}{the value of the \code{start$bw} argument (or its default value if missing).}
 \item{startbOpt}{the value of the \code{start$bOpt} argument (or its default value if missing).}
 \item{startRange}{the value of the \code{start$range} argument (or its default value if missing).}
 \item{provMethod}{the value of the \code{test$method} argument (or its default value if missing).}
 \item{provDist}{the value of the \code{test$priorDist} argument (or its default value if missing).}
 \item{provPar}{the value of the \code{test$priorPar} argument (or its default value if missing).}
 \item{provRange}{the value of the \code{test$range} argument (or its default value if missing).}
 \item{provD}{the value of the \code{test$D} argument (or its default value if missing).}
 \item{itemSelect}{the value of the \code{test$itemSelect} argument (or its default value if missing).}
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
 
The function \code{print.cat} returns similar (but differently organized) results. 
 }


\details{
 The \code{randomCAT} function generates an adaptive test using an item bank specified by argument \code{itemBank}, and for
 a given true ability level specified by argument \code{trueTheta}. The maximal length of the test can be fixed through
 the \code{maxItems} argument, with a default value of 50 items. 

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
 \item{\code{bw}: the bandwidth value, used to define the range of ability levels for selecting the first items 
	(default is 4). Ignored if either \code{fixItems} or \code{seed} is not \code{NULL}. See \code{\link{startItems}} for further details.}
 \item{\code{range}: the maximal range of difficulty levels, specified as a vector of two numeric values (default is \code{c(-4,4)}). 
	Ignored if either \code{fixItems} or \code{seed} is not \code{NULL}. See \code{\link{startItems}} for further details.}
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
	\code{priorPar=c(0,1)}).}
 \item{\code{range}: the maximal range of ability levels, set as a vector of two numeric components. The ability estimate
	will always lie to this interval (set by default to [-4, 4]). Ignored if \code{method=="EAP"}.}
 \item{\code{D}: the value of the metric constant. Default is \code{D=1} for logistic metric. Setting \code{D=1.702}
	yields approximately the normal metric (Haley, 1952).}
 \item{eapPar}{numeric: vector of three components, holding respectively the values of the arguments \code{lower}, 
               \code{upper} and \code{nqp} of the \code{\link{eapEst}} command. Default vector is (-4, 4, 20). 
               Ignored if \code{method} is not \code{"EAP"}.}
 \item{itemSelect}{character: the rule for next item selecion, with values \code{"info"} (default) for maximum information criterion
			 or \code{"Owen"} for Owen's procedure. See \code{\link{nextItem}} for further details.}
 }
 These arguments are passed to the functions \code{\link{thetaEst}} (or \code{\link{eapEst}}) and \code{\link{semTheta}} (or\code{\link{eapSem}})
 to estimate the ability level and the standard error of this estimate. 

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
 argument \code{itemSelect} of the \code{test} list is not allowed into the \code{final} list.

 If some arguments of these lists are missing, they are automatically set to their default value. The contents
 of the lists is checked with the \code{\link{testList}} function, and the adaptive test is generated only if the lists
 are adequately defined. Othgerwise, a message error is printed.

 The function \code{plot.cat} represents the set of provisional and final ability estimates througghout the test. Corresponding
 confidence intervals (with confidence level defined by the argument \code{alpha}) are also drawn if \code{ci=TRUE} (which
 is the default value). In addition, a classification threshold can be displayed by specifying a numeric value to the
 argument \code{thr}. The default value \code{NULL} does not represent any threshold.
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

 Lord, F.M. (1980). \emph{Applications of item response theory to practical testing problems}.
 Hillsdale, NJ: Lawrence Erlbaum.

 Owen, R. J. (1975). A Bayesian sequential procedure for quantal response in the context of adaptive mental testing. 
 \emph{Journal of the American Statistical Association, 70}, 351-356.

 Warm, T.A. (1989). Weighted likelihood estimation of ability in item response models. \emph{Psychometrika, 54}, 427-450.
}


\author{
    David Magis \cr
    Post-doctoral researcher, FNRS (Fonds National de la Recherche Scientifique) \cr
    Department of Mathematics, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{testList}}, \code{\link{startItems}}, \code{\link{nextItem}}, \code{\link{thetaEst}}, \code{\link{semTheta}}
} 

\examples{
\dontrun{
 # Loading the 'tcals' parameters 
 data(tcals)
 tcals <- as.matrix(tcals)
 
 # Item bank creation with 'tcals' item parameters
 bank <- createItemBank(tcals)

 # Creation of a starting list: 5 items, initial theta 0, bw 2
 start<-list(nrItems=5, theta=0, bw=2)

 # Creation of 'test' and 'final' lists: weighted likelihood
 # estimation of ability (both provisional and final), and
 # maximum information criterion for next item selection
 test<-list(method="WL")
 final<-test

 # Creation of a stopping rule: precision criterion, standard
 # error to be reached 0.3
 stop<-list(rule="precision", thr=0.3)

 # CAT test
 res<-randomCAT(0, bank, start=start, test=test, stop=stop,
 	final=final)

 # New 'test' and 'final' rules (BM and EAP estimation
 # with Jeffreys' prior)
 test2<-list(method="BM", priorDist="Jeffreys")
 final2<-list(method="EAP", priorDist="Jeffreys")

 # New stopping rule: classification criterion, with 
 # classification threshold 0 and alpha level 0.05
 stop2<-list(rule="classification", thr=0, alpha=0.05)

 # CAT test with new 'test', 'stop' and 'final' rules
 res2<-randomCAT(0, bank, start=start, test=test2, stop=stop2,
 	final=final2)

 # New stopping rule: classification criterion, with 
 # classification threshold 0.5 and alpha level 0.05
 stop3<-list(rule="classification", thr=0.5, alpha=0.05)

 # CAT test with new 'stop' rule
 res3<-randomCAT(0, bank, start=start, test=test2, stop=stop3,
 	final=final2)

 # Plotting results
 plot(res)
 plot(res, ci=FALSE)
 plot(res, thr=0)
 plot(res, thr=0.5)

 plot(res2,thr=0)

 plot(res3, thr=0.5)
 }
 }

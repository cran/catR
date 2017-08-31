\name{randomCAT}
\alias{randomCAT}
\alias{print.cat}
\alias{plot.cat}

\title{Random generation of adaptive tests (dichotomous and polytomous models)}

\description{
 This command generates a response pattern to an adaptive test, for a given item bank (with either dichotomous or polytomous models), a true ability level, and several lists of CAT parameters (starting items, stopping rule, provisional and final ability estimators).
 }

\usage{
randomCAT(trueTheta, itemBank, model = NULL, responses = NULL, genSeed = NULL, 
 	cbControl = NULL, nAvailable = NULL, start = list(fixItems = NULL, 
 	seed = NULL, nrItems = 1, theta = 0, D = 1, randomesque = 1, 
 	random.seed = NULL, startSelect = "MFI", cb.control = FALSE, random.cb = NULL), 
 	test = list(method = "BM", 	priorDist = "norm", priorPar = c(0, 1), 
 	range = c(-4, 4), D = 1, parInt = c(-4, 4, 33), itemSelect = "MFI", 
 	infoType = "observed", randomesque = 1, random.seed = NULL, AP = 1, 
 	proRule = "length", proThr = 20,constantPatt = NULL), 
 	stop = list(rule = "length", thr = 20, alpha = 0.05),
 	final = list(method = "BM", priorDist = "norm", 	priorPar = c(0, 1), 
 	range = c(-4, 4), D = 1, parInt = c(-4, 4, 33), 	alpha = 0.05), 
 	allTheta = FALSE, save.output = FALSE, output = c("path", "name", "csv"))
\method{print}{cat}(x, ...)
\method{plot}{cat}(x, ci = FALSE, alpha = 0.05, trueTh = TRUE, classThr = NULL, 
 	save.plot = FALSE, save.options = c("path", "name", "pdf"),...)
 }

\arguments{
 \item{trueTheta}{numeric: the value of the true ability level.}
 \item{itemBank}{numeric: a suitable matrix of item parameters (possibly augmented by group membership for content balancing). See \bold{Details}.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
\item{responses}{either \code{NULL} (default) or a vector of pre-specified item responses with as many components as the rows of \code{itemBank}. See \bold{Details}.}
\item{genSeed}{either a numeric value to fix the random generation of responses pattern or \code{NULL} (default). Ignored if \code{responses} is not \code{NULL}. See \bold{Details}.}
 \item{cbControl}{either a list of accurate format to control for content balancing, or \code{NULL}. See \bold{Details}.}
 \item{nAvailable}{either a boolean vector indicating which items (denoted by 1's) are available at the start of the test and which (denoted by 0's) are not, or \code{NULL} (default). See \bold{Details}.}
 \item{start}{a list with the options for starting the adaptive test. See \bold{Details}.}
 \item{test}{a list with the options for provisional ability estimation and next item selection. See \bold{Details}.}
 \item{stop}{a list with the options of the stopping rule. See \bold{Details}.}
 \item{final}{a list with the options for final ability estimation. See \bold{Details}.}
\item{allTheta}{logical: should all provisional ability estimates and standard errors be computed and returned (including among the starting items)? Default is \code{FALSE}, meaning that provisional ability estimates and standard errors are computed only after the selection of the starting items. Ignored if the \code{$nrItems} of the \code{start} list is equal to zero.}
 \item{save.output}{logical: should the output be saved in an external text file? (default is \code{FALSE}).}
 \item{output}{character: a vector of three components. The first component is either the file path to save the output of \code{"path"} (default), the second component is the name of the output file, and the third component is the file type, either \code{"txt"} or \code{"csv"} (default). See \bold{Details}.}
 \item{x}{an object of class "cat", typically an output of \code{randomCAT} function.}
 \item{ci}{logical: should the confidence intervals be plotted for each provisional ability estimate? (default is \code{TRUE}).}
 \item{alpha}{numeric: the significance level for provisional confidence intervals (default is 0.05). Ignored if \code{ci} is \code{FALSE}.}
 \item{trueTh}{logical: should the true ability level be drawn by a horizontal line? (default is TRUE).}
 \item{classThr}{either a numeric value giving the classification threshold to be displayed, or \code{NULL}.}
 \item{save.plot}{logical: should the plot be saved in an external figure? (default is \code{FALSE}).}
 \item{save.options}{character: a vector of three components. The first component is either the file path or \code{"path"} (default), the second component is the name of the output file or ,\code{"name"} (default), and the third component is the file extension, either \code{"pdf"} (default) or \code{"jpeg"}. Ignored if \code{save.plot} is \code{FALSE}. See \bold{Details}.}
 \item{...}{other generic arguments to be passed to \code{print} and \code{plot} functions.}
}

\value{
 The function \code{randomCAT} returns a list of class "cat" with the following arguments:
 \item{trueTheta}{the value of the \code{trueTheta} argument.}
\item{model}{the value of the \code{model} argument.}
 \item{testItems}{a vector with the items that were administered during the test.}
 \item{itemPar}{a matrix with the parameters of the items administered during the test.}
\item{itemNames}{either a vector wit the names of the selected items during the CAT, or \code{NULL}.}
 \item{pattern}{the generated response pattern (as vector of 0 and 1 entries).}
 \item{thetaProv}{a vector with the provisional ability estimates.}
 \item{seProv}{a vector with the standard errors of the provisional ability estimates.}
 \item{ruleFinal}{either the stopping rule(s) that was (were) satisfied to make the CAT stop, or \code{NULL}.}
 \item{thFinal}{the final ability estimate.}
 \item{seFinal}{the standrad error of the final ability estimate.}
 \item{ciFinal}{the confidence interval of the final ability estimate.}
\item{genSeed}{the value of the \code{genSeed} argument.}
 \item{startFixItems}{the value of the \code{start$fixItems} argument (or its default value if missing).}
 \item{startSeed}{the value of the \code{start$seed} argument (or its default value if missing).}
 \item{startNrItems}{the value of the \code{start$nrItems} argument (or its default value if missing).}
 \item{startTheta}{the value of the \code{start$theta} argument (or its default value if missing).}
 \item{startD}{the value of the \code{start$D} argument (or its default value if missing).}
 \item{startRandomesque}{the value of the \code{start$randomesque} argument (or its default value if missing).}
 \item{startThStart}{the starting ability values used for selecting the first items of the test.}
 \item{startSelect}{the value of the \code{start$startSelect} argument (or its default value if missing).}
\item{startCB}{logical value, being \code{TRUE} if both \code{cbControl} is not \code{NULL} and \code{start$cb.control} is \code{TRUE}, and \code{FALSE} otherwise.}
 \item{provMethod}{the value of the \code{test$method} argument (or its default value if missing).}
 \item{provDist}{the value of the \code{test$priorDist} argument (or its default value if missing).}
 \item{provPar}{the value of the \code{test$priorPar} argument (or its default value if missing).}
 \item{provRange}{the value of the \code{test$range} argument (or its default value if missing).}
 \item{provD}{the value of the \code{test$D} argument (or its default value if missing)or \code{NA} if \code{model} is not \code{NULL}.}
 \item{itemSelect}{the value of the \code{test$itemSelect} argument (or its default value if missing).}
 \item{infoType}{the value of the \code{test$infoType} argument (or its default value if missing). Not returned if \code{model} is not \code{NULL}.}
 \item{randomesque}{the value of the \code{test$randomesque} argument (or its default value if missing).}
 \item{AP}{the value of the \code{test$AP} argument (or its default value if missing).}
\item{constantPattern}{the value of the \code{test$constantPatt} argument (or its default value if missing).}
 \item{cbControl}{the value of the \code{cbControl} argument (or its default value if missing).}
 \item{cbGroup}{the value of the \code{itemBank$cbGroup} element of the item bank \code{itemBank} (for dichotomous IRT models), or the \code{cbGroup} element returned by the \code{\link{breakBank}} function (for polytomous IRT models), or \code{NULL}.}
 \item{stopRule}{the value of the \code{stop$rule} argument (or its default value if missing).}
 \item{stopThr}{the value of the \code{stop$thr} argument (or its default value if missing).}
 \item{stopAlpha}{the value of the \code{stop$alpha} argument (or its default value if missing).}
 \item{endWarning}{a logical indactor indicating whether the adaptive test stopped because the stopping rule(s) was (were) satisfied, or becasue all items in the bank were administered.}
 \item{finalMethod}{the value of the \code{final$method} argument (or its default value if missing).}
 \item{finalDist}{the value of the \code{final$priorDist} argument (or its default value if missing).}
 \item{finalPar}{the value of the \code{final$priorPar} argument (or its default value if missing).}
 \item{finalRange}{the value of the \code{final$range} argument (or its default value if missing).}
 \item{finalD}{the value of the \code{final$D} argument (or its default value if missing), or \code{NA} if \code{model} is not \code{NULL}.}
 \item{finalAlpha}{the value of the \code{final$alpha} argument (or its default value if missing).}
 \item{save.output}{the value of the \code{save.output} argument.}
 \item{output}{the value of the \code{output} argument.}
\item{assigned.responses}{a logical value, being \code{TRUE} if \code{responses} was provided or \code{FALSE} \code{responses} was set to \code{NULL}.}
 
The function \code{print.cat} returns similar (but differently organized) results. 
 }


\details{
 The \code{randomCAT} function generates an adaptive test using an item bank specified by arguments \code{itemBank} and \code{model}, and for a given true ability level specified by argument \code{trueTheta}. 

Dichotomous IRT models are considered whenever \code{model} is set to \code{NULL} (default value). In this case, \code{itemBank} must be a matrix with one row per item and four columns, with the values of the discrimination, the difficulty, the pseudo-guessing and the inattention parameters (in this order). These are the parameters of the four-parameter logistic (4PL) model (Barton and Lord, 1981). 

Polytomous IRT models are specified by their respective acronym: \code{"GRM"} for Graded Response Model, \code{"MGRM"} for Modified Graded Response Model, \code{"PCM"} for Partical Credit Model, \code{"GPCM"} for Generalized Partial Credit Model, \code{"RSM"} for Rating Scale Model and \code{"NRM"} for Nominal Response Model. The \code{itemBank} still holds one row per item, end the number of columns and their content depends on the model. See \code{\link{genPolyMatrix}} for further information and illustrative examples of suitable polytomous item banks.

By default all item responses will be randomly drawn from parent distribution set by the item bank parameters of the \code{itemBank} matrix (using the \code{\link{genPattern}} function for instance). Moreover, the random generation of the item responses can be fixed (for e.g., replication purposes) by assigning some numeric value to the \code{genSeed} argument. By default this argument is equal to \code{NULL} so the random seed is not fixed (and two successive runs of \code{randomCAT} will usually lead to different response patterns).

It is possible, however, to provide a full response pattern of previously recorded responses to each of the item bank, for instance for post-hoc simulations. This is done by providing to the \code{responses} argument a vector of binary entries (without missing values). By default \code{responses} is set to \code{NULL} and item responses will be drawn from the item bank parameters.

With the aforementioned item bank structures, content balancing cannot be controled and \code{cbControl} must be set to \code{NULL} (default value). Otherwise this will most often lead to an error. In order to allow for content balancing control:
\enumerate{
\item{the \code{itemBank} must be updated with an additional column holding the group names;}
\item{the \code{cbControl} argument must be set properly as a list with group names and theoretical proportions for content balancing.}
}
See the \code{\link{nextItem}} function for further details on how to specify \code{cbControl} properly and under which conditions it is operational (see Kingsbury and Zara, 1989, for further details). Separation of the item parameters and the vector of group membership is performed internally through the \code{\link{breakBank}} function (and thus should not be performed prior to CAT generation).

 The test specification is made by means of four lists of options: one list for the selection of the starting items, 
 one list with the options for provisional ability estimation, one list to define the stopping rule, and one list with
 the options for final ability estimation. These lists are specified respectively by the arguments \code{start}, \code{test}, \code{stop} and \code{final}. 

 The \code{start} list can contain one or several of the following arguments:
 \itemize{
 \item{\code{fixItems}: either a vector of integer values, setting the items to be administered as first items, or \code{NULL} (default) to let the function select the items.}
 \item{\code{seed}: either a numeric value to fix the random seed for item selection, \code{NA} to randomly select the items withour fixing the random seed, or \code{NULL} (default) to select the items on the basis of their difficulty level. Ignored if \code{fixItems} is not \code{NULL}.}
 \item{\code{nrItems}: numeric, the number of starting items to be randomly selected (default is 1). Can be equal to zero to avoid initial selection of items (see \bold{Details}). Used only if \code{fixItems} is \code{NULL} and \code{seed} is not \code{NULL}.}
 \item{\code{theta}: numeric, a vector of the initial ability levels for selecting the first items (default is the single value 0). Ignored if either \code{fixItems} or \code{seed} is not \code{NULL}. See \code{\link{startItems}} for further details.}
 \item{\code{D}: numeric, the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952). Ignored if \code{model} is not \code{NULL} and if \code{startSelect} is not \code{"MFI"}.} 
 \item{\code{randomesque}: integer, the number of 'randomesque' items to be picked up optimally for each value of \code{theta} vector, before random selection of a single one. Ignored if either \code{fixItems} or \code{seed} is not \code{NULL}. See \code{\link{startItems}} for further details.}
\item{\code{random.seed}: either \code{NULL} (default) or a numeric value to fix the random seed of randomesque selection of the items. Ignored if either \code{fixItems} or \code{seed} is not \code{NULL}.}
 \item{\code{startSelect}: the method for selecting the first items of the test, with possible values \code{"bOpt"} and \code{"MFI"} (default). Ignored if either \code{fixItems} or \code{seed} is not \code{NULL}. See \code{\link{startItems}} for further details.}
\item{\code{cb.control}: logical value indicating whether control for content balancing should also be done when selecting the starting items. Default is \code{FALSE}. Ignored if argument \code{cbControl} is \code{NULL}.}
\item{\code{random.cb}: either \code{NULL} (default) or a numeric value to fix the selection of subgroups of items in case of content balancing control with starting items. Ignored if either \code{cbControl} is \code{NULL} or if \code{start$cb.control} is \code{FALSE}.}
 }
 These arguments are passed to the function \code{\link{startItems}} to select the first items of the test.

If the argument \code{nrItems} is set to zero, then no starting item is selected and the adaptive process starts with a provisional ability level equal to the value of argument \code{theta} (or its default). Moreover, the likelihood function is then set as a flat, uniform function on the whole ability range. See the \code{\link{nextItem}} function for further details.

 The \code{test} list can contain one or several of the following arguments:
 \itemize{
 \item{\code{method}: a character string to specify the method for ability estimation. Possible values are: \code{"BM"} (default) for Bayesian modal estimation (Birnbaum, 1969), \code{"ML"} for maximum likelihood estimation (Lord, 1980), \code{"EAP"} for expected a posteriori (EAP) estimation (Bock and Mislevy, 1982), and \code{"WL"} for weighted likelihood estimation (Warm, 1989).}
 \item{\code{priorDist}: a character string which sets the prior distribution. Possible values are: \code{"norm"} (default) for normal distribution, \code{"unif"} for uniform distribution, and \code{"Jeffreys"} for Jeffreys' noninformative prior distribution (Jeffreys, 1939, 1946). ignored if \code{method} is neither \code{"BM"} nor \code{"EAP"}.}
 \item{\code{priorPar}: a vector of two numeric components, which sets the parameters of the prior distribution. If (\code{method="BM"} or \code{method=="EAP"}) and \code{priorDist="norm"}, the components of \code{priorPar} are respectively the mean and the standard deviation of the prior normal density. If (\code{method="BM"} or \code{method="EAP"}) and \code{priorDist="unif"}, the components of \code{priorPar} are respectively the lower and upper bound of the prior uniform density. Ignored in all other cases. By default, \code{priorPar} takes the parameters of the prior standard normal distribution (i.e., \code{priorPar=c(0,1)}). In addition, \code{priorPar} also provides the prior parameters for the comoutation of MLWI and MPWI values for next item selection (see \code{\link{nextItem}} for further details).}
 \item{\code{range}: the maximal range of ability levels, set as a vector of two numeric components. The ability estimate will always lie to this interval (set by default to [-4, 4]). Ignored if \code{method=="EAP"}.}
 \item{\code{D}: the value of the metric constant. Default is \code{D=1} for logistic metric. Setting \code{D=1.702} yields approximately the normal metric (Haley, 1952). Ignored if \code{model} is not \code{NULL}.}
 \item{\code{parInt}: a numeric vector of three components, holding respectively the values of the arguments \code{lower}, \code{upper} and \code{nqp} of the \code{\link{eapEst}}, \code{\link{eapSem}} and \code{MWI} commands. It specifies the range of quadrature points for numerical integration, and is used for computing the EAP estimate, its standard error, and the MLWI and MPWI values for next item selection. Default vector is (-4, 4, 33), thus setting the range from -4 to 4 by steps of 0.25. Ignored if \code{method} is not \code{"EAP"} and if \code{itemSelect} is neither \code{"MLWI"} nor \code{"MPWI"}.}
 \item{\code{itemSelect}: the rule for next item selecion, with possible values \code{"MFI"} (default) for maximum Fisher information criterion; \code{"bOpt"} for optimal ability-difficulty match (or Urry's procedure) (not available if \code{model} is not \code{NULL}); \code{"thOpt"} for optimal theta selection (not available if \code{model} is not \code{NULL}); \code{"MLWI"} and \code{"MPWI"} for respectively maximum likelihood and posterior weighted information criterion; \code{"MEPV"} for minimum expected posterior variance; \code{"MEI"} for maximum expected information; \code{"KL"} and \code{"KLP"} for Kullback-Leibler and posterior Kullback-Leibler information methods; \code{"progressive"} and \code{"proportional"} for progressive and proportional methods; ; and \code{"random"} for random selection. For further details, see \code{\link{nextItem}}.}
 \item{\code{infoType}: character: the type of information function to be used for next item selection. Possible values are \code{"observed"} (default) for observed information function, and \code{"Fisher"} for Fisher information function. Ignored if \code{itemselect} is not \code{"MEI"}.}
 \item{\code{randomesque}: integer: the number of items to be chosen from the next item selection rule, among those the next item to be administered will be randomly picked up. Default value is 1 and leads to usual selection of the optimal item (Kingsbury and Zara, 1989).}
\item{\code{random.seed}: either \code{NULL} (default) or a numeric value to fix the random seed of randomesque selection of the items. Ignored if either \code{fixItems} or \code{seed} is not \code{NULL}.}
\item{\code{AP}: the acceleration parameter required for progressive and proportional methods, with default value 1. Ignored with all other selection methods.}
\item{\code{proRule}: the stopping rule considered for progressive and proportional methods, with possible values \code{"length"} (default), \code{"precision"} or both. Ignored with all other selection methods.}
\item{\code{proThr}: the stopping rule threshold considered for progressive and proportional methods. Default value is 20. Ignored with all other selection methods.}
\item{\code{constantPatt}: character: the method to estimate ability in case of constant pattern (i.e. only correct or only incorrect responses). Can be either \code{NULL} (default), \code{"BM"}, \code{"EAP"}, \code{"WL"}, \code{"fixed4"}, \code{"fixed7"} or \code{"var"}. \emph{Currently only implemented for dichotomous IRT models}.}
 }
 These arguments are passed to the functions \code{\link{thetaEst}} and \code{\link{semTheta}} to estimate the ability level and the standard error of this estimate. In addition, some arguments are passed to \code{\link{nextItem}} to select the next item appropriately.

 The \code{stop} list can contain one or several of the following arguments:
 \itemize{
 \item{\code{rule}: a vector of character strings specifying the stopping rules. Possible values are: \code{"length"} (default), to stop the test after a pre-specified number of items administered; \code{"precision"}, to stop the test 
when the provisional standard error of ability becomes less than or equal to the pre-specified value; \code{"classification"}, for which the test ends whenever the provisional confidence interval (set by the \code{alpha} argument) does not hold the classification threshold anymore (this is also called the ACI rule; see e.g. Thomson, 2009); and \code{"minInfo"} to stop the test if the maximum item information of the available items at current ability estimate is smaller than the prespecified threshold. Can take a single value.}
 \item{\code{thr}: a vector of numeric values fixing the threshold(s) of the stopping rule(s). If \code{rule="length"}, \code{thr} is the maximal number of items to be administered. If \code{rule="precision"}, \code{thr} is the precision level (i.e. the standard error) to be reached before stopping. If \code{rule="classification"}, \code{thr} corresponds to the ability level which serves as a classification rule (i.e. which must not be covered by the provisional confidence interval). Finally, If \code{rule="minInfo"}, \code{thr} corresponds to the minimum item information that can be observed in the bank of remaining available items. The \code{"classification"} rule is not available for the progressive and proportional item selection rules.}
 \item{\code{alpha}: the significance (or \eqn{\alpha}) level for computing the priovisional confidence interval of ability.
Ignored if \code{rule} is not \code{"classification"}. \emph{Important:} the \code{thr} value smust be sorted in the same order of appearance as the \code{rule} methods. }
 }

 Eventually, the \code{final} list can contain one or several arguments of the \code{test} list (with possiblly different values), as well as the additional \code{alpha} argument. The latter specifies the \eqn{\alpha} level of the final confidence interval of ability, which is computed as
 \deqn{[\hat{\theta}-z_{1-\alpha/2} \; se(\hat{\theta}) ; \hat{\theta}+z_{1-\alpha/2} \; se(\hat{\theta})]}
 where \eqn{\hat{\theta}} and \eqn{se(\hat{\theta})} are respectively the ability estimate and its standard error. Note that the argument \code{itemSelect} of the \code{test} list is not used for final estimation of the ability level, and is therefore not allowed into the \code{final} list.

 If some arguments of these lists are missing, they are automatically set to their default value. The contents of the lists is checked with the \code{\link{testList}} function, and the adaptive test is generated only if the lists are adequately defined. Otherwise, a message error is printed. Note that the \code{\link{testList}} function works for both dichotomous and polytomous models.

Usually the ability estimates and related standard errors are computed right after the administration of the starting items (that is, if \emph{k} starting items are administered, the first \emph{(k-1)} ability levels and standard errors are missing). This can however be avoided by fixing the argument \code{allTheta} to \code{TRUE} (by default it is \code{FALSE}). In this case, all provisional ability estimates and standard errors are computed and returned, but in the display of th output file, the first \emph{(k-1)} abilities and standard errors are printed in parentheses (otherwise they are returned as \code{NA} values). Note that \code{allTheta} is ignored if no starting item was selected (that is, if argument \code{nrItems} of the \code{start} list is set to zero).

 The output of \code{randomCAT}, as displayed by the \code{print.cat} function, can be stored in a text file provided that  \code{save.output} is set to \code{TRUE} (the default value \code{FALSE} does not execute the storage). In this case, the \code{(output} argument mus hold three character values: the path to where the output file must be stored, the name of  the output file, and the type of output file. If the path is not provided (i.e. left to its default value \code{"path"}), it will be saved in the default working directory. The default name is \code{"name"}, and the default type is \code{"csv"}. Any other value yields a text file. See the \bold{Examples} section for an illustration. 

 The function \code{plot.cat} represents the set of provisional and final ability estimates througghout the test. Corresponding confidence intervals (with confidence level defined by the argument \code{alpha}) are also drawn if \code{ci=TRUE} (which is \emph{not} the default value), except when stepsize adjustment was made for constant patterns (as the standard error cannot be estimated at this stage). The true ability level can be drawn by a horizontal solid line by specifying \code{trueTh=TRUE} (which is the default value); setting it to \code{FALSE} will undo the drawing. Finally, any classification threshold can be additionally displayed by specifying a numeric value to the argument \code{classThr}. The default value \code{NULL} does not display any threshold.

 Finally, the plot can be saved in an external file, either as PDF or JPEG format. First, the argument \code{save.plot} must be set to \code{TRUE} (default is \code{FALSE}). Then, the file path for figure storage, the name of the figure and its format are specified through the argument \code{save.options}, all as character strings. See the \code{Examples} section for further information and a practical example.
}

\references{
 Barrada, J. R., Olea, J., Ponsoda, V., and Abad, F. J. (2010). A method for the comparison of item selection rules in computerized adaptive testing. \emph{Applied Psychological Measurement}, \emph{20}, 213-229. \doi{10.1177/0146621610370152}

 Barton, M.A., and Lord, F.M. (1981). \emph{An upper asymptote for the three-parameter logistic item-response model}.
 Research Bulletin 81-20. Princeton, NJ: Educational Testing Service.

 Birnbaum, A. (1969). Statistical theory for logistic mental test models with a prior distribution of ability.
 \emph{Journal of Mathematical Psychology, 6}, 258-276. \doi{10.1016/0022-2496(69)90005-4}

 Bock, R. D., and Mislevy, R. J. (1982). Adaptive EAP estimation of ability in a microcomputer environment. 
 \emph{Applied Psychological Measurement, 6}, 431-444. \doi{10.1177/014662168200600405}

 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error}. Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 
 
 Jeffreys, H. (1939). \emph{Theory of probability}. Oxford, UK: Oxford University Press.

 Jeffreys, H. (1946). An invariant form for the prior probability in estimation problems. \emph{Proceedings of the Royal Society of London. Series A, Mathematical and Physical Sciences, 186}, 453-461.

 Kingsbury, G. G., and Zara, A. R. (1989). Procedures for selecting items for computerized adaptive tests. \emph{Applied Measurement in Education, 2}, 359-375. \doi{10.1207/s15324818ame0204_6}

 Lord, F. M. (1980). \emph{Applications of item response theory to practical testing problems}. Hillsdale, NJ: Lawrence Erlbaum.

Magis, D. and Barrada, J. R. (2017). Computerized Adaptive Testing with R: Recent Updates of the Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{Code Snippets}, \emph{76(1)}, 1-18. \doi{10.18637/jss.v076.c01}

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. \doi{10.18637/jss.v048.i08}

Thompson, N. A. (2009). Item selection in computerized classification testing. \emph{Educational and Psychological Measurement, 69}, 778-793. \doi{10.1177/0013164408324460}

 Urry, V. W. (1970). \emph{A Monte Carlo investigation of logistic test models.} Unpublished doctoral dissertation. West Lafayette, IN: Purdue University.

 van der Linden, W. J. (1998). Bayesian item selection criteria for adaptive testing. \emph{Psychometrika, 63}, 201-216. \doi{10.1007/BF02294775}

 Veerkamp, W. J. J., and Berger, M. P. F. (1997). Some new item selection criteria for adaptive testing. \emph{Journal of Educational and Behavioral Statistics, 22}, 203-226. \doi{10.3102/10769986022002203}

 Warm, T.A. (1989). Weighted likelihood estimation of ability in item response models. \emph{Psychometrika, 54}, 427-450. \doi{10.1007/BF02294627}
}


\author{
    David Magis \cr
    Department of Education, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}

 Juan Ramon Barrada \cr
Department of Psychology and Sociology, Universidad Zaragoza, Spain \cr
\email{barrada@unizar.es}

 }

\seealso{
\code{\link{testList}}, \code{\link{startItems}}, \code{\link{nextItem}}, \code{\link{thetaEst}}, \code{\link{semTheta}}, \code{\link{eapEst}},
\code{\link{eapSem}}, \code{\link{MWI}}, \code{\link{MEI}}, \code{\link{KL}}, \code{\link{breakBank}}, \code{\link{genPolyMatrix}},
\code{\link{genPattern}}
} 

\examples{

## Dichotomous models ##

 # Loading the 'tcals' parameters 
 data(tcals)

 # Matrix of item parameters (with no content balancing)
 bank <- as.matrix(tcals[,1:4])
 
 # Creation of a starting list: 5 items with evenly spread theta values
 start <- list(theta = -2:2)

 # Creation of 'test' list: weighted likelihood estimation of provisional ability, 
 # and MFI criterion for next item selection
 test <- list(method = "WL", itemSelect = "MFI")

 # Creation of 'final' list: EAP estimation of final ability
 final <- list(method = "EAP")

 # Creation of a stopping rule: precision criterion, standard error to be reached 0.3
 stop <- list(rule = "precision", thr = 0.3)

 # CAT test
 res <- randomCAT(0, bank, start = start, test = test, stop = stop, final = final)

\dontrun{
 # Update of the stopping rule: by adding a length criterion, with threshold of 15 items
 stop <- list(rule = c("precision", "length"), thr = c(0.3, 15))

 # CAT test
 res <- randomCAT(0, bank, start = start, test = test, stop = stop, final = final)


 # Creation of a (purely artificial) response pattern for post-hoc simulation
 resp <- rbinom(nrow(bank), 1, 0.5)
 res.ph <- randomCAT(0, bank, start = start, test = test, stop = stop, final = final,
        responses = resp)

 # New 'test' and 'final' rules (BM and EAP estimation with Jeffreys' prior, 
 # randomesque value 5)
 test2 <- list(method = "BM", priorDist = "Jeffreys", randomesque = 5)
 final2 <- list(method = "EAP", priorDist = "Jeffreys")

 # New stopping rule: classification criterion, with classification threshold 0 and 
 # alpha level 0.05
 stop2 <- list(rule = "classification", thr = 0, alpha = 0.05)

 # CAT test with new 'test', 'stop' and 'final' rules
 res2 <- randomCAT(0, bank, start = start, test = test2, stop = stop2, final = final2)

 # New stopping rule: classification criterion, with classification threshold 0.5
 # and alpha level 0.05
 stop3 <- list(rule = "classification", thr = 0.5, alpha = 0.05)

 # CAT test with new 'stop' rule
 res3 <- randomCAT(0, bank, start = start, test = test2, stop = stop3, final=final2)

 # new 'test' and 'stop' rule for next item selection
 test3 <- list(method = "WL", itemSelect = "MLWI")
 stop4 <- list(rule = "length",thr = 10)
 res4 <- randomCAT(0, bank, start = start, test = test3, stop = stop4, final = final2)

 # Creation of the 'cbList' list with arbitrary proportions
 cbList <- list(names = c("Audio1", "Audio2", "Written1", "Written2",
        "Written3"), props = c(0.1, 0.2, 0.2, 0.2, 0.3))

 # CAT test with 'start', 'test', 'stop4' and 'final2' lists and content balancing 
 # using 'cbList' ('tcals' must be used now for content balancing)
 stop4 <- list(rule = "length",thr = 10)
 res5 <- randomCAT(0, tcals, start = start, test = test, stop = stop4, final = final2, 
      cbControl = cbList)

 # new 'start' list to force content balancing control at the starting step
 start2 <- list(seed = 1, nrItems = 3, cb.control = TRUE)
 res6 <- randomCAT(0, tcals, start = start2, test = test, stop = stop4, final = final2, 
      cbControl = cbList)

 # Using progressive item selection rule and requiring all ability estimates and SEs
 test4 <- list(itemSelect = "progressive")
 res6 <- randomCAT(0, tcals, start = start, test = test4, stop = stop4, final = final, 
      cbControl = cbList, allTheta = TRUE)

 # Saving the output in the external 'out' text file within folder 'Program Files' 
 # of hard drive 'C'
 res5 <- randomCAT(0, tcals, start = start, test = test, stop = stop4, final = final2, 
      cbControl = cbList, save.output = TRUE, 
      output = c("c:/Program Files/", "out", "txt"))

 # Plotting results
 plot(res)
 plot(res, ci = TRUE)
 plot(res, ci = TRUE, trueTh = FALSE)
 plot(res, ci = TRUE, classThr = 1)

 # Saving last figure into PDF file 'figure' within folder 'C:/Program Files/'
 plot(res, ci = TRUE, classThr = 1, save.plot = TRUE,
      save.options = c("c:/Program Files/", "figure", "pdf"))

 # With mistake
 plot(res, ci = 0.05)
 plot(res, classThr = TRUE)


## Polytomous models ##

 # Generation of an item bank under GRM with 100 items and at most 4 categories
 m.GRM <- genPolyMatrix(100, 4, "GRM")
 m.GRM <- as.matrix(m.GRM)

 # CAT options
 start <- list(theta = c(-1, 0),startSelect = "MFI")
 test <- list(method = "BM", itemSelect = "KL")
 stop <- list(rule = "precision", thr = 0.4)
 final <- list(method = "EAP")

 # CAT test
 res <- randomCAT(0, m.GRM,model = "GRM", start = start, test = test, stop = stop,
 	final = final)
 res


 # Creation of an appropriate list for content balancing
 # Equal proprotions across subgroups of items
 cbList <- list(names = c("Group1", "Group2", "Group3", "Group4"), props = rep(1,4))

 # With content balancing, all ability estimates and progressive method
 m.GRM <- genPolyMatrix(100, 4, "GRM", cbControl = cbList)
 test <- list(method = "BM", itemSelect = "progressive")
 res <- randomCAT(0, m.GRM, model = "GRM", start = start, test = test, stop = stop, 
 	final = final, cbControl = cbList, allTheta = TRUE)
 res

 # Loading the cat_pav data
 data(cat_pav)
 cat_pav <- as.matrix(cat_pav)

 stop <- list(rule = "length", thr = 10)
 res <- randomCAT(0, cat_pav, model = "GPCM", start = start, test = test, stop = stop, 
 	final = final, allTheta = TRUE)
 res
 }
}

\name{simulateRespondents}
\alias{simulateRespondents}
\alias{print.catResult}
\alias{plot.catResult}

\title{Simulation of multiple examinees of adaptive tests}

\description{
This command runs a set of adaptive tests, for a given item bank, a set of ability levels, a possible matrix of item responses, and several lists of CAT parameters (starting items, stopping rule, provisional and final ability estimators).}


\usage{
simulateRespondents(thetas, itemBank, responsesMatrix = NULL, model = NULL, 
  	genSeed = NULL, cbControl = NULL, rmax = 1, Mrmax = "restricted", 
  	start = list(fixItems = NULL, seed = NULL, nrItems = 1, theta = 0, 
  	D = 1, randomesque = 1, random.seed = NULL, startSelect = "MFI", 
  	cb.control = FALSE, random.cb = NULL), test = list(method = "BM", 
  	priorDist = "norm", priorPar = c(0,1), weight = "Huber", tuCo = 1, 
  	sem.type = "classic", sem.exact = FALSE, se.ase = 10, range = c(-4, 4), 
  	D = 1, parInt = c(-4, 4, 33), itemSelect = "MFI", infoType = "observed", 
  	randomesque = 1, random.seed = NULL, AP = 1, proRule = "length", 
  	proThr = 20, constantPatt = NULL), stop = list(rule = "length", 
  	thr = 20, alpha = 0.05), final = list(method = "BM", priorDist = "norm", 
  	priorPar = c(0,1), weight = "Huber", tuCo = 1, sem.type = "classic", 
  	sem.exact = FALSE, range = c(-4, 4), D = 1, parInt = c(-4, 4, 33), 
  	alpha = 0.05), save.output = FALSE, output = c("", "catR", "csv"))
\method{print}{catResult}(x, ...)
\method{plot}{catResult}(x, type = "all", deciles = "theta", save.plot = FALSE, 
  	save.options = c("", "plot", "pdf"), res = 300, ...)
 }

\arguments{
 \item{thetas}{numeric: a vector of true ability values for which a CAT must be generated for each component.}
 \item{itemBank}{numeric: a suitable matrix of item parameters (possibly augmented by group membership for content balancing). See \bold{Details}.}
 \item{responsesMatrix}{numeric: either \code{NULL} (default) or a suitable matrix of item responses. See \bold{Details}.}
 \item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. See \bold{Details}.}
\item{genSeed}{either a vector of numeric values to fix the random seed of each generated pattern, or \code{NULL} (default). Ignored if \code{responsesMatrix} is not \code{NULL}. See \bold{Details}.}
 \item{cbControl}{either a list of accurate format to control for content balancing, or \code{NULL}. See \bold{Details}.}
\item{rmax}{numeric: the maximum exposure rate (default is 1).}
\item{Mrmax}{character: the method for controlling maximum exposure rate. Possible values are \code{"restricted"}(default) or \code{"IE"}. See \bold{Details}.}
 \item{start}{a list with the options for starting the adaptive test. See \bold{Details}.}
 \item{test}{a list with the options for provisional ability estimation and next item selection. See \bold{Details}.}
 \item{stop}{a list with the options of the stopping rule. See \bold{Details}.}
 \item{final}{a list with the options for final ability estimation. See \bold{Details}.}
 \item{save.output}{logical: should the output be saved in an external text file? (default is \code{FALSE}).}
 \item{output}{character: a vector of three components. The first component is either the file path to save the output or \code{""} (default), the second component is either the initial part of the name of the output file or \code{"catR"}, and the third component is the file type, either \code{"txt"} or \code{"csv"} (default). See \bold{Details}.}
 \item{x}{an object of class "cat", typically an output of \code{simulateRespondents} function.}
 \item{type}{character: the type of plot to display. Possible values are \code{"all"} (default), \code{"trueEst"}, \code{"expRate"}, \code{"cumExpRate"},
\code{"cumNumberItems"}, \code{"expRatePara"}, \code{"condBias"}, \code{"condRMSE"}, \code{"numberItems"}, \code{"sError"} and \code{="condThr"}. See \bold{Details}.}
 \item{deciles}{whether the deciles number (\code{"deciles"}) or the mean ability level per decile (\code{"theta"}) will be used in the axis of the plots.}
 \item{save.plot}{logical: should the plot be saved in an external figure? (default is \code{FALSE}).}
 \item{save.options}{character: a vector of three components. The first component is either the file path or \code{""} (default), the second component is the name of the output file or ,\code{"plot"} (default), and the third component is the file extension, either \code{"pdf"} (default) or \code{"jpeg"}. Ignored if \code{save.plot} is \code{FALSE}. See \bold{Details}.}
 \item{res}{numeric: the resolution for JPEG figures (default value is 300).}
 \item{...}{other generic arguments to be passed to \code{print} and \code{plot} functions.}
}


\value{
 The function \code{simulateRespondents} returns a list of class "catResult" with the following arguments:
\item{thetas}{the value of the \code{thetas} argument.}
\item{itemBank}{the value of the \code{itemBank} argument.}
\item{responsesMatrix}{the value of the \code{responsesMatrix} argument.}
\item{model}{the value of the \code{model} argument.}
\item{genSeed}{the value of the \code{genSeed} argument.}
\item{cbControl}{the value of the \code{cbControl} argument.}
\item{rmax}{the value of the \code{rmax} argument.}
\item{Mrmax}{the value of the \code{Mrmax} argument.}
\item{start}{the value of the \code{start} argument.}
\item{test}{the value of the \code{test} argument.}
\item{stop}{the value of the \code{stop} argument.}
\item{final}{the value of the \code{final} argument.}
\item{save.output}{the value of the \code{save.output} argument.}
\item{output}{the value of the \code{output} argument.}
\item{estimatedThetas}{a vector with (final) estimated ability levels.} 
\item{correlation}{the correlation between the \code{thetas} vector and estimated ability levels.} 
\item{bias}{the value of the bias between true and estimated ability levels.} 
\item{RMSE}{the value of the RMSE between true and estimated ability levels.}
\item{thrOK}{a vector indicating whether the respondents finished the test satisfying (1) or not (0) the stop criteria.}
\item{exposureRates}{a vector with empirical exposure rates of all item in the bank.} 
\item{testLength}{the mean test length.} 
\item{overlap}{the item overlap rate.}
\item{numberItems}{a vector with the lengths of each adaptive test (i.e. the number of items administered).}
\item{condTheta}{a vector with the mean ability level per decile.}
\item{condBias}{a vector with conditional mean bias per decile.}
\item{condRMSE}{a vector with conditional RMSE per decile.}
\item{condnItems}{a vector with conditional mean test length per decile.}
\item{condSE}{a vector with conditional mean standard error per decile.}
\item{condthrOK}{a vector with conditional proportion of respondents that finish the test satisfying the stop criteria per decile.}
\item{ndecile}{a vector with the number of respondents per decile.}
\item{final.values.df}{a data frame with true ability levels, final ability estimates and standard errors, and test lengths.}
 \item{responses.df}{a data frame with all items administered, all item responses and all provisional ability estimates. -99 is displayed when the actual test length is smaller than the maximum test length for those item positions where no item was administered as the stop criterion was already reached.}
\item{start.time}{the CPU time at the start of the CAT generation.}
\item{finish.time}{the CPU time at the end of all CAT generations.}
}


\details{
 The \code{simulateRespondents} function permits to generate several adaptive tests to a set of respondents defined by their ability levels. It makes a repeated call to an adaptive test using an item bank specified by arguments \code{itemBank} \code{model}, and with the same \code{start}, \code{test}, \code{stop} and \code{final} lists. Content balancing can also be controlled for each respondent with the \code{cbControl} argument. All arguments of \code{simulateRespondents} are used in exactly the same manner as in \code{\link{randomCAT}} (so refer to this function for further information), except the following four.

First, \code{thetas} is now a vector of ability levels, and a CAT will be generated for each component of \code{thetas}. If \code{responsesMatrix} is \code{NULL}, item responses are generated from the IRT model, item bank parameters, and the ability levels. In this case, \code{thetas} can be considered as real ability levels. Otherwise, \code{responsesMatrix} must be provided as a matrix with as many rows as the length of \code{thetas} and as many columns as the number of items in \code{itemBank}. Each row contains the response pattern of one examinee whose ability level is given by the value of the corresponding component of \code{thetas}. Note that only allowable item responses can be included in \code{responsesMatrix} and missing values are not accepted. Fixing the random seed can also be done with the \code{genSeed} argument. The latter must hold as many components as the vector \code{thetas}, otherwise an error message is returned. Each component of \code{genSeed} is used to fix the seed for ach pattern generation.

The option of providing a response matrix through \code{responsesMatrix} is considered for two possible uses: 
\enumerate{
\item{\emph{post-hoc simulations}: examinees provided responses to the full item bank and one wants to test the performance of a CAT with those responses,}
\item{simulations considering \emph{item parameter estimation errors}: responses to the full bank are generated with the real parameters and the CAT is run is the estimated parameters and the responses from the correct model.}
}
Note that if \code{thetas} holds a single value, then the function simply calls \code{\link{randomCAT}} and returns its output instead.

Second, \code{rmax} fixes the desired maximum exposure rate for all items. Default value is 1, allowing thus items to be administered to all respondents without restrictions.

Third, the \code{Mrmax} argument fixes the method to constraint exposure rates to be smaller than the maximum allowed rate. Possible methods are the restricted method (\code{"restricted"}; Revuelta and Ponsoda, 1998) and the item-eligibility method (\code{"IE"}; van der Linden and Veldkamp, 2004). A description of both methods can be found in Barrada, Abad and Veldkamp (2009)

Fourth, if the length of \code{thetas} is greater than 1, \code{save.output} for \code{\link{randomCAT}} is fixed to \code{"FALSE"}. Otherwise, the same file would be overwritten for each new respondent.

The output of \code{simulateRespondents}, as displayed by the \code{print.catResult} function, proposes summary statistics related to overall accuracy (bias, RMSE, etc), conditional accuracy per decile, and item exposure control (minimum and maximum exposure rates, test overlap rate, etc) among all generated CATs. This output can be saved when \code{save.output} is set to \code{TRUE}, and if so, three outpout files are returned: (1) one with the main summary statistics, with overall results and conditional on decile results; (2) the file with respondents' patterns, items administered and provisional ability estimates; and (3) a table with true and estimated ability levels, final standard errors and numbers of items administered per respondent. Specific information provided and saved depends on the \code{rule} used in \code{stop}.

This output can also be graphically displayed with the \code{plot.catResult} function. In addition to the output of the function and the item bank, it takes the argument \code{type} to determine which plot should be returned. Ten different single plots can be displayed:
\enumerate{
\item{\code{"trueEst"}: the scatterplot of true vs. estimated ability levels.}
\item{\code{"expRate"}: the exposure rates of the items, having ranked them according to these exposure rates.}
\item{\code{"cumExpRate"}: the cumulative exposure rates of the items, having ranked them according to these exposure rates.}
\item{\code{"cumNumberItems"}: the test length as a function of cumulative percent of examinees. This plot is not available when \code{rule} is 'length'.}
\item{\code{"expRatePara"}: the scatterplot of item exposure rates vs. item discrimination parameters. This plot is not available when \code{model} is 'PCM' or 'NRM', as in those IRT models there are no discrimination parameters.}
\item{\code{"condBias"}: the conditional bias of ability estimation as a function of the deciles of the true ability levels.}
\item{\code{"condRMSE"}: the conditional RMSE of ability estimation as a function of the deciles of the true ability levels.}
\item{\code{"numberItems"}: the conditional test length as a function of the deciles of the true ability levels. This plot is not available when \code{rule} is 'length'.}
\item{\code{"sError"}: the conditional standard error of ability estimation as a function of the deciles of the true ability levels.}
\item{\code{"condThr"}: the conditional proportions of CATs satisfying the 'precision' or 'classification' stopping rule, as a function of the deciles of the true ability levels. This plot is not available when \code{rule} is 'length'.}
}
In addition, the value \code{"all"} (default value) displays several available plots in a single panel. Displayed plots depend on the \code{rule} used in \code{stop}: (a) with 'length', \code{"trueEst"}, \code{"condBias"}, \code{"condRMSE"}, \code{"expRate"}, \code{"cumExpRate"}, \code{"expRatePara"}; (b) with 'precision' and 'classification', all the plots but \code{"cumExpRate"}.

These plots can be saved as external PDF or JPEG files, by setting \code{save.plot} to \code{TRUE} and defining the arguments of \code{save.options} accurately (see \code{\link{randomCAT}} for further explanations and the \bold{Examples} section below).
}

\references{
 Barrada, J. R., Abad, F. J., and Veldkamp, B. P. (2009). Comparison of methods for controlling maximum exposure rates in computerized adaptive testing. \emph{Psicothema}, \emph{21}, 313-320.

Magis, D. and Barrada, J. R. (2017). Computerized Adaptive Testing with R: Recent Updates of the Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{Code Snippets}, \emph{76(1)}, 1-18. \doi{10.18637/jss.v076.c01}

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. \doi{10.18637/jss.v048.i08}
 
Revuelta, J., and Ponsoda, V. (1998). A comparison of item exposure control methods in computerized adaptive testing. \emph{Journal of Educational Measurement, 35}, 311-327. \doi{10.1111/j.1745-3984.1998.tb00541.x}

van der Linden, W.J., and Veldkamp, B.P. (2004). Constraining item exposure in computerized adaptive testing with shadow tests. \emph{Journal of Educational and Behavioral Statistics, 29}, 273-291. \doi{10.3102/10769986029003273}
}


\author{
Juan Ramon Barrada \cr
Department of Psychology and Sociology, Universidad Zaragoza, Spain \cr
\email{barrada@unizar.es}

    David Magis \cr
    Department of Psychology, University of Liege, Belgium \cr
    \email{david.magis@uliege.be}

Guido Corradi \cr
Department of Psychology and Sociology, Universidad Zaragoza, Spain \cr
\email{guidocor@gmail.com}

 }

\seealso{
\code{\link{randomCAT}}
} 

\examples{

## Dichotomous IRT model ##

 # Loading the 'tcals' parameters 
 data(tcals)
 bank <- as.matrix(tcals[,1:4])
 
 # Creation of a starting list with three theta values
 start <- list(theta = -1:1, randomesque = 5)

 # Creation of 'test' list: maximum likelihood estimation and
 # progressive method
 test <- list(method = "ML", itemSelect = "progressive")

 # Creation of a stopping rule: precision criterion, standard
 # error to be reached 0.3
 stop <- list(rule = "precision", thr = 0.3)

 # Creation of 'final' list: ML estimation of final ability
 final <- list(method = "ML")

 # Generation of ten respondents
 set.seed(1)
 thetas <- rnorm(10)

 # Default CAT generations, output not saved
 res <- simulateRespondents(thetas, bank, start = start, test = test, stop = stop, 
  final = final)

 # Maximum exposure restricted to 0.8
 res2 <- simulateRespondents(thetas, bank, start = start, test = test, stop = stop, 
  final = final, rmax = 0.8)

\dontrun{

 # Output saved
 res3 <- simulateRespondents(thetas, bank, start = start, test = test, stop = stop, 
  final = final, save.output = TRUE, output = c("C:/Program Files/", "out", "txt"))
}

 # With content balancing #

 # Creation of an appropriate list for content balancing
 # Equal proportions across subgroups of items
 cbList <- list(names = c("Audio1", "Audio2", "Written1", "Written2", "Written3"), 
        props = c(0.1, 0.2, 0.2, 0.2, 0.3))

 # CAT test (same options as above)
 res4 <- simulateRespondents(thetas, tcals, start = start, test = test, stop = stop, 
  final = final, cbControl = cbList)

 # Plotting and saving output #

 # Plotting all possible panels
 plot(res)
 plot(res, deciles = "deciles")

\dontrun{

 # Saving the plot in the "fig" pdf file in "c:/Program Files/"
 plot(res, save.plot = TRUE, save.options = c("c:/Program Files/", "fig", "pdf"))
}

 # Plotting the 'trueEst' type of plot
 plot(res, type = "trueEst")
 }

\name{checkStopRule}
\alias{checkStopRule}

\title{Checking whether the stopping rule is satisfied}

\description{
 This command tests whether one of the specified stopping rules is satisfied in order to stop the CAT.
 }

\usage{
checkStopRule(th, se, N, it = NULL, model = NULL, D = 1, stop)
 }

\arguments{
 \item{th}{numeric: the current ability estimate.}
 \item{se}{numeric: the current standard error estimate.}
 \item{N}{numeric: the current number of administered items.}
 \item{it}{the matrix of parameters of currently available items (default is \code{NULL}). Ignored if element \code{$rule} of \code{stop} list does not hold the value \code{"minInfo"}. See \bold{Details}.}
\item{model}{either \code{NULL} (default) for dichotomous models, or any suitable acronym for polytomous models. Possible values are \code{"GRM"}, \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"}, \code{"RSM"} and \code{"NRM"}. Ignored if element \code{$rule} of \code{stop} list does not hold the value \code{"minInfo"}. See \bold{Details}.}
 \item{D}{numeric: the metric constant. Default is \code{D=1} (for logistic metric); \code{D=1.702} yields approximately the normal metric (Haley, 1952). Ignored if \code{model} is not \code{NULL} or if element \code{$rule} of \code{stop} list does not hold the value \code{"minInfo"}. } 
 \item{stop}{a list with valuable element names and contents to set the stopping rule, as an input of the \code{randomCAT} or \code{simulateRespondents} functions.}
}

\value{
 A list with two arguments:
 \item{decision}{a logical value indicating whether at least one of the stopping rules was satisfied (\code{TRUE}) or not (\code{FALSE}).}
 \item{rule}{either a vector with the names of the stopping rules that were satisfied, or \code{NULL} if \code{decision} is \code{FALSE}.}
 }
 

\details{
 The \code{checkStopRule} function checks whether at least one of the stopping rules was satisfied at the current step of the CAT test process. It mainly serves as an internal application for\code{\link{randomCAT}} function. 

 The \code{stop} list must be supplied accordingly to the requested list in the \code{randomCAT()} and in the \code{simulateRespondents()} functions.

Three input values must be supplied: \code{th} as the current ability estimate; \code{se} as the current standard error value related to \code{th} estimate; and \code{N} as the current test length (i.e., number of administered items). In addition, if the \code{stop$rule} vector holds the option \code{"minInfo"}, three additional input value smust be supplied: \code{it} with the item parameters or all available items in the bank (i.e., previously administered items should not be set as input); \code{model} to specify the type of IRT model, either dichotomous or polytomous (see \code{\link{Pi}} fir further details); and possibly the scaling constant \code{D} set to one by default.

All stopping rules are being assessed and if at least one of them is satisfied, the output list will hold the vector of rules's nmaes that were satisfied through the \code{rule} argument). If none of the stopping rules were satisfied, this \code{rule} output argument is simply \code{NULL}. 
 }


\references{
 Haley, D.C. (1952). \emph{Estimation of the dosage mortality relationship when the dose is subject to error}. Technical report no 15. Palo Alto, CA: Applied Mathematics and Statistics Laboratory, Stanford University. 

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. URL \url{http://www.jstatsoft.org/v48/i08/}
}

\author{
    David Magis \cr
    Department of Education, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{randomCAT}}, \code{\link{simulateRespondents}}, \code{\link{Pi}}
} 

\examples{

# Creation of a 'stop' list with two possible rules
 stop <- list(rule = c("length", "precision"), thr = c(20, 0.3))

# Example of successful 'length' rule
 checkStopRule(th = 0.35, se = 0.41, N = 20, stop = stop)

# Example of successful 'precision' rule
 checkStopRule(th = 0.35, se = 0.29, N = 15, stop = stop)

# Example of jointly successful 'length' and 'precision' rules
 checkStopRule(th = 0.35, se = 0.29, N = 20, stop = stop)

# Example without sucessfull rule
 checkStopRule(th = 0.35, se = 0.31, N = 18, stop = stop) 

# Creation of a short bank of available items under 2PL
 it <- genDichoMatrix(items = 5, model = "2PL", seed = 1)

# Computation of maximum information at ability level 0.35
maxI <- max(Ii(0.35, it)$Ii)

# Creation of a 'stop' list with four possible rules and too large threshold for 'minInfo'
 stop <- list(rule = c("length", "precision", "classification", "minInfo"), 
              thr = c(20, 0.3, 1, maxI-0.01), alpha = 0.05)

# Example with sucessfull 'classification' rule only
 checkStopRule(th = 0.35, se = 0.31, N = 18, it = it, stop = stop) 

# Creation of a 'stop' list with four possible rules and too large threshold for 'minInfo'
 stop <- list(rule = c("length", "precision", "classification", "minInfo"), 
              thr = c(20, 0.3, 1, maxI+0.01), alpha = 0.05)

# Example with sucessfull 'minInfo' rule only
 checkStopRule(th = 0.35, se = 0.55, N = 18, it = it, stop = stop) 

}

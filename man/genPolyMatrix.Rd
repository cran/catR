\name{genPolyMatrix}
\alias{genPolyMatrix}

\title{Item bank generation (polytomous models)}

\description{
This command generates an item bank from prespecified parent distributions for use with polytomous IRT models. Subgroups of items can also be specified for content balancing purposes.
 }

\usage{
genPolyMatrix(items = 100, nrCat = 3, model = "GRM", seed = 1, same.nrCat = FALSE,
 	cbControl = NULL)
 }

\arguments{
 \item{items}{integer: the number of items to generate (default is 100).}
\item{nrCat}{integer: the (maximum) number of response categories to generate (default is 3).}
\item{model}{character: the type of polytomous IRT model. Possible values are \code{"GRM"} (default), \code{"MGRM"}, \code{"PCM"}, \code{"GPCM"} and \code{"NRM"}. See \bold{Details}.}
\item{seed}{numeric: the random seed for item parameter generation (default is 1).}
\item{same.nrCat}{logical: should all items have the same number of response categories? (default is \code{FALSE}. Ignored if \code{model} is either \code{"MGRM"} or \code{"RSM"}. See \bold{Details}.}
 \item{cbControl}{either a list of accurate format to control for content balancing, or \code{NULL}. See \bold{Details}.}
}

\value{
 A matrix with \code{items} rows and as many columns as required for the considered IRT model: 
\itemize{
\item{\eqn{\max_j \,g_j+1} columns, holding parameters \eqn{(\alpha_j, \beta_{j1}, ..., \beta_{j,g_j})} if \code{model} is \code{"GRM"};}
\item{\eqn{g+2} columns, holding parameters \eqn{(\alpha_j, b_j, c_1, ..., c_g)} if \code{model} is \code{"MGRM"};}
\item{\eqn{\max_j \,g_j} columns, holding parameters \eqn{(\delta_{j1}, ..., \delta_{j,g_j})}  if \code{model} is \code{"PCM"};}
\item{\eqn{\max_j \,g_j+1} columns, holding parameters \eqn{(\alpha_j, \delta_{j1}, ..., \delta_{j,g_j})} if \code{model} is \code{"GPCM"};}
\item{\eqn{g+1} columns, holding parameters \eqn{(\lambda_j, \delta_1, ..., \delta_g)} if \code{model} is \code{"RSM"};}
\item{\eqn{2\,\max_j\, g_j} columns, holding parameters \eqn{(\alpha_{j1}, c_{j1}, \alpha_{j2}, c_{j2}, ..., \alpha_{j,g_j}, c_{j, g_j})} if \code{model} is \code{"NRM"}.}
}
If \code{cbControl} is not \code{NULL}, the output matrix contains an additional colum for item membership is included.
 }
 
\details{
The \code{genPolyMatrix} permits to quickly generate a polytomous item bank in suitable format for further use in e.g. computing item response probabilities with the \code{\link{Pi}}. 

The six polytomous IRT models that are supported are: 
\enumerate{
\item{the \emph{Graded Response Model} (GRM; Samejima, 1969);}
\item{the \emph{Modified Graded Response Model} (MGRM; Muraki, 1990);}
\item{the \emph{Partial Credit Model} (PCM; Masters, 1982);}
\item{the \emph{Generalized Partial Credit Model} (GPCM; Muraki, 1992);}
\item{the \emph{Rating Scale Model} (RSM; Andrich, 1978);}
\item{the \emph{Nominal Response Model} (NRM; Bock, 1972).}
}
Each model is specified through the \code{model} argument, with its accronym surrounded by double quotes (i.e. \code{"GRM"} for GRM, \code{"PCM"} for PCM, etc.). The default value is \code{"GRM"}.

For any item \eqn{j}, set \eqn{(0, ..., g_j)} as the \eqn{g_j+1} possible response categories. The maximum number of response categories can differ across items under the GRM, PCM, GPCM and NRM, but they are obviously equal across items under the MGRM and RSM. In the latter, set \eqn{g} as the (same) number of response categories for all items. It is possible however to require all items to have the same number of response categories, by fixing the \code{same.nrCat} argument to \code{TRUE}.

In case of GRM, PCM, GPCM or NRM with \code{same.nrCat} being \code{FALSE}, the number of response categories \eqn{g_j+1} per item is drawn from a Poisson distribution with parameter \code{nrCat}, and this number is restricted to the interval [2; \code{nrCat}]. This ensure at least two response categories and at most \code{nrCat} categories. In all other cases, each \eqn{g_j+1} is trivially fixed to \eqn{g+1 = } \code{nrCat}.

Denote further \eqn{P_{jk}(\theta)} as the probability of answering response category \eqn{k \in \{0, ..., g_j\}} of item \eqn{j}. For GRM and MGRM, response probabilities \eqn{P_{jk}(\theta)} are defined through cumulative probabilities, while for PCM, GPCM, RSM and NRM they are directly computed.

For GRM and MGRM, set \eqn{P_{jk}^*(\theta)} as the (cumulative) probability of asnwering response category \eqn{k} or "above", that is \eqn{P_{jk}^*(\theta) = Pr(X_j \geq k | \theta)} where \eqn{X_j} is the item response. It follows obviously that for any \eqn{\theta}, \eqn{P_{j0}^*(\theta) = 1} and \eqn{P_{jk}^*(\theta) = 0} when \eqn{k>g_j}. Furthermore, response category probabilities are found back by the relationship \eqn{P_{jk}(\theta)= P_{jk}^*(\theta)-P_{j,k+1}^*(\theta)}. Then, the GRM is defined by (Samejima, 1969)
\deqn{P_{jk}^*(\theta)=\frac{\exp\,[\alpha_j\,(\theta-\beta_{jk})]}{1+\exp\,[\alpha_j\,(\theta-\beta_{jk})]}} 
and the MGRM by (Muraki, 1990)
\deqn{P_{jk}^*(\theta)=\frac{\exp\,[\alpha_j\,(\theta-b_j+c_k)]}{1+\exp\,[\alpha_j\,(\theta-b_j+c_k)]}.}

The PCM, GPCM, RSM and NRM are defined as "divide-by-total" models (Embretson and Reise, 2000). The PCM has following response category probability (Masters, 1982):
\deqn{P_{jk}(\theta)=\frac{\exp\,\sum_{t=0}^k (\theta-\delta_{jt})}{\sum_{r=0}^{g_j}\,\exp\, \sum_{t=0}^r (\theta-\delta_{jt})}\quad \mbox{with} \quad  \sum_{t=0}^0 (\theta-\delta_{jt})=0.}

The GPCM has following response category probability (Muraki, 1992):
\deqn{P_{jk}(\theta)=\frac{\exp\,\sum_{t=0}^k \alpha_j\,(\theta-\delta_{jt})}{\sum_{r=0}^{g_j}\,\exp\, \sum_{t=0}^r \alpha_j\,(\theta-\delta_{jt})}\quad \mbox{with} \quad  \sum_{t=0}^0 \alpha_j\,(\theta-\delta_{jt})=0.}

The RSM has following response category probability (Andrich, 1978):
\deqn{P_{jk}(\theta)=\frac{\exp\,\sum_{t=0}^k [\theta-(\lambda_j+\delta_t)]}{\sum_{r=0}^{g_j}\,\exp\, \sum_{t=0}^r [\theta-(\lambda_j+\delta_t)]}\quad \mbox{with} \quad  \sum_{t=0}^0 [\theta-(\lambda_j+\delta_t)]=0.}

Finally, the NRM has following response category probability (Bock, 1972):
\deqn{P_{jk}(\theta)=\frac{\exp (\alpha_{jk}\,\theta+c_{jk})}{\sum_{r=0}^{g_j} \exp (\alpha_{jr}\,\theta+c_{jr})}\quad \mbox{with} \quad  \alpha_{j0}\,\theta+c_{j0}=0.}

The following parent distributions are considered to generate the different item parameters. The \eqn{\alpha_j} parameters of GRM, MGRM and GPCM, as well as the \eqn{\alpha_{jk}} parameters of the NRM, are drawn from a log-normal distribution with mean 0 and standard deviation 0.1225. All other parameters are drawn from a standard normal distribution. Moreover, the \eqn{\beta_{jk}} parameters of the GRM and the \eqn{c_k} parameters of the MGRM are sorted respectively in increasing and decreasing order of \eqn{k}, to ensure decreasing trend in the cumulative \eqn{P_{jk}^*(\theta)} probabilities.

The output is a matrix with one row per item and as many columns as required to hold all item parameters. In case of missing response categories, the corresponding parameters are replaced by \code{NA} values. Column names refer to the corresponding model parameters. See \bold{Details} for further explanations and \bold{Examples} for illustrative examples.

Finally, the output matrix can contain an additional vector with the names of the subgroups to be used for content balancing purposes. To do so, the argument \code{cbControl} (with default value is \code{NULL}) must contain a list of two elements: (a) the \code{names} element with the names of the subgroups, and (b) the \code{props} elements with proportions of items per subgroup (of the same length of \code{names} element, with only positive numbers but not necessarily summing to one). The \code{cbControl} argument is similar to the one in \code{\link{nextItem}} and \code{\link{randomCAT}} functions to control for content balancing. The output matrix contains then an additional column, with the names of the subgroups randomly allocated to each item by using random multinomial draws with the probabilities given by \code{cbControl$props}. 
 }


\references{
Andrich, D. (1978). A rating formulation for ordered response categories. \emph{Psychometrika, 43}, 561-573. doi: 10.1007/BF02293814

Bock, R. D. (1972). Estimating item parameters and latent ability when responses are scored in two or more nominal categories. \emph{Psychometrika, 37}, 29-51. doi: 10.1007/BF02291411

Embretson, S. E., and Reise, S. P. (2000). \emph{Item response theory for psychologists}. Mahwah, NJ: Lawrence Erlbaum Associates.

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. URL \url{http://www.jstatsoft.org/v48/i08/}

Masters, G. N. (1982). A Rasch model for partial credit scoring. \emph{Psychometrika, 47}, 149-174. doi: 10.1007/BF02296272

Muraki, E. (1990). Fitting a polytomous item response model to Likert-type data. \emph{Applied Psychological Measurement, 14}, 59-71. doi: 10.1177/014662169001400106

Muraki, E. (1992). A generalized partial credit model: Application of an EM algorithm. \emph{Applied Psychological Measurement, 16}, 19-176. doi: 10.1177/014662169201600206

Samejima, F. (1969). \emph{Estimation of latent ability using a response pattern of graded scores}. Psychometrika Monograph (vol. 17).
}



\author{
    David Magis \cr
    Department of Education, University of Liege, Belgium \cr
    \email{david.magis@ulg.ac.be}
 }

\seealso{
\code{\link{Pi}}
} 

\examples{

# All generated item banks have 10 items and at most four response categories

 # GRM 
 genPolyMatrix(10, 4, model = "GRM")

 # GRM with same number of response categories
 genPolyMatrix(10, 4, model = "GRM", same.nrCat = TRUE)

 # MGRM 
 genPolyMatrix(10, 4, model = "MGRM")

 # MGRM with same number of response categories
 genPolyMatrix(10, 4, model = "MGRM", same.nrCat = TRUE) # same result

 # PCM 
 genPolyMatrix(10, 4, model = "PCM")

 # PCM with same number of response categories
 genPolyMatrix(10, 4, model = "PCM", same.nrCat = TRUE) 

 # GPCM 
 genPolyMatrix(10, 4, model = "GPCM")

 # GPCM with same number of response categories
 genPolyMatrix(10, 4, model = "GPCM", same.nrCat = TRUE) 

 # RSM 
 genPolyMatrix(10, 4, model = "RSM")

 # RSM with same number of response categories
 genPolyMatrix(10, 4, model = "RSM", same.nrCat = TRUE) # same result

 # NRM 
 genPolyMatrix(10, 4, model = "NRM")

 # NRM with same number of response categories
 genPolyMatrix(10, 4, model = "NRM", same.nrCat = TRUE)  

## Content balancing

 # Creation of the 'cbList' list with arbitrary proportions
 cbList <- list(names = c("Audio1", "Audio2", "Written1", "Written2", "Written3"), 
        props = c(0.1, 0.2, 0.2, 0.2, 0.3))

 # NRM with 100 items
 genPolyMatrix(100, 4, model = "NRM", cbControl = cbList)
 }

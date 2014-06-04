\name{tcals}
\alias{tcals}

\docType{data}

\title{Items parameters of the TCALS 1998 data set and subgroups of items}

\description{
 The TCALS (\emph{Test de Connaissance en Anglais Langue Seconde}) is an aptitude test of English language as a second language
 in the French speaking college of Outaouais (Gatineau, QC, Canada). The test consists of 85 items and is administered 
 every year to newly incoming students. The item parameters of the year 1998 have been estimated under the 3PL model.
 Inattention parameters are therefore fixed to one. Subgroups of items are also included for content balancing purposes.
}

\format{
 A matrix with 85 rows and five columns, respectively holding the discrimination, difficulty, pseudo-guessing and inattention parameters 
 as calibrated on the results of the 1998 application of the TCALS questionnaire. The fifth column holds the name of the subgroups of items:
 \itemize{
 \item{\emph{Audio1}: listening comprehension of sentences (items \bold{1} to \bold{12}).}
 \item{\emph{Audio2}: listening comprehension of dialogs and short texts (items \bold{13} to \bold{33}).}
 \item{\emph{Written1}: written vocabulary exercises (items \bold{34} to \bold{46}).}
 \item{\emph{Written2}: written grammar exercises (items \bold{47} to \bold{63}).}
 \item{\emph{Written3}: written exercises of other types: reading and mistake detection (items \bold{64} to \bold{85}).}
 }
}

\source{
 The TCALS test was originally developed by Laurier, Froio, Pearo and Fournier (1998) and item parameters were obtained from Raiche (2002).
}


\references{
 Laurier, M., Froio, L., Pearo C., and Fournier, M. (1998). \emph{Test de classement d'anglais langue seconde au collegial}. 
 Montreal, Canada: College de Maisonneuve.

 Magis, D., and Raiche, G. (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the R 
 Package \emph{catR}. \emph{Journal of Statistical Software}, \emph{48 (8)}, 1-31. URL \url{http://www.jstatsoft.org/v48/i08/}

 Raiche, G. (2002). \emph{Le depistage du sous-classement aux tests de classement en anglais, langue seconde, au collegial
 [The detection of under classification to the collegial English as a second language placement tests]}. Gatineau, QC: College de l'Outaouais.
}


\keyword{datasets}
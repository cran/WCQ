\name{wcq}
\alias{wcq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{WCQ QTL mapping}
\description{Performs QTL mapping using Chen-Qin two sample mean test.}

\usage{
wcq(marker_data, trait_data, alleles)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{marker_data}{m x n matrix/array containing marker score data; m is the number of markers, n is the number of samples}
  \item{trait_data}{t x n matrix/array containing trait value data; t is the number of traits, n is the number of samples}
  \item{alleles}{vector containing the marker scores representing the alleles; currently supports diploids only; third entry assumed to be score for heterozygous marker}
}
\value{
\item{pval_matrix}{m x t array/matrix containing p-values of each marker as being a potential QTL for a trait}
}

\author{Jan Michael Yap}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(sample.markers)
data(sample.traits)
alleles <- c(1,2,3)
wcq(sample.markers, sample.traits, alleles)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{QTL mapping}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

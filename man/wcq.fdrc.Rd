\name{wcq.fdrc}
\alias{wcq.fdrc}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{WCQ QTL detection with false discovery rate control}

\description{Performs WCQ with optional false discovery rate control.}

\usage{
wcq.fdrc(marker_data, trait_data, alleles, fdrc.method = "none")
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{marker_data}{m x n matrix/array containing marker score data; m is the number of markers, n is the number of samples}
  \item{trait_data}{t x n matrix/array containing trait value data; t is the number of traits, n is the number of samples}
  \item{alleles}{vector containing the marker scores representing the alleles; currently supports diploids only; third entry assumed to be score for heterozygous marker}
  \item{fdrc.method}{string value representing method to be used for false discovery rate control; uses the p.adjust method in the R default library; default value is "none"}
}

\value{
\item{adjusted_pval_matrix}{m x t array/matrix containing adjusted p-values of each marker as being a potential QTL for a trait; if fdrc.method is set to "none", returns the pval_matrix computed using the wcq function}
}

\author{Jan Michael Yap}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{p.adjust}
}

\examples{
data(sample.markers)
data(sample.traits)
alleles <- c(1,2,3)
wcq.fdrc(sample.markers, sample.traits, alleles, fdrc.method="bonferroni")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QTL detection }
\keyword{ False discovery rate control }

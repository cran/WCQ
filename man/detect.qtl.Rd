\name{detect.qtl}
\alias{detect.qtl}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Detect markers with QTL effects}
\description{Performs detection of QTL effects on list of markers on a number of traits. Optional false discovery rate control can also be performed.}

\usage{
detect.qtl(marker_data, trait_data, alleles, fdrc.method = "none", threshold = 0.05)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{marker_data}{m x n matrix/array containing marker score data; m is the number of markers, n is the number of samples}
  \item{trait_data}{t x n matrix/array containing trait value data; t is the number of traits, n is the number of samples}
  \item{alleles}{vector containing the marker scores representing the alleles; currently supports diploids only; third entry assumed to be score for heterozygous marker}
  \item{fdrc.method}{string value representing method to be used for false discovery rate control; uses the p.adjust method in the R default library; default value is "none"}
  \item{threshold}{p-value threshold to be used for determining significant QTL effect; default value is 0.05}

}

\value{
  \item{qtl.list}{list containing markers with detected significant QTL effect; has length equal to the number of traits}
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
qlist <- detect.qtl(sample.markers, sample.traits, alleles, fdrc.method="bonferroni")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QTL detection }
\keyword{ False discovery rate control }

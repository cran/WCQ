\name{onesamplemeantest}
\alias{onesamplemeantest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Chen-Qin Multidimensional One-Sample Mean Test}
\description{
Multidimensional one-sample mean test used by WCQ.
}
\usage{
onesamplemeantest(X1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X1}{matrix/array containing population data to be tested}
}
\value{
 \item{pval}{p-value of the test}
}
\references{Chen,Y.X., Qin,Y.L. (2010) A Two-Sample Test for High-Dimensional Data with Applications to Gene-Set Testing, The Annals of Statistics, 38(2), 808-835.}
\author{Jan Michael Yap}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
sample.data <- array(rnorm(100),c(10,10))
onesamplemeantest(sample.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ One-sample Mean Test }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

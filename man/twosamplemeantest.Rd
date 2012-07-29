\name{twosamplemeantest}
\alias{twosamplemeantest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Chen-Qin Multidimensional Two-Sample Mean Test}
\description{Multidimensional two-sample test of means used by WCQ}
\usage{
twosamplemeantest(X1, X2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X1,X2}{matrices/arrays containing population values for the test}
}
\value{
 \item{pval}{p-value of the test}
}
\references{Chen,Y.X., Qin,Y.L. (2010) A Two-Sample Test for High-Dimensional Data with Applications to Gene-Set Testing, The Annals of Statistics, 38(2), 808-835.}
\author{Jan Michael Yap}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
sample.data1 <- array(rnorm(100),c(10,10))
sample.data2 <- array(rnorm(100),c(10,10))
twosamplemeantest(sample.data1, sample.data2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Two-sample mean test }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

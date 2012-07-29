\name{sample.markers}
\alias{sample.markers}
\docType{data}
\title{Dummy Marker Data}
\description{
Dummy marker data generated using R/qtl package. 1 and 2 correspond to the marker scores of homozygous alleles, while 3 is for the heterozygous allele.
}
\usage{data(sample.markers)}
\format{
  The format is:
 int [1:200, 1:20] 3 2 2 3 3 3 2 3 3 3 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:200] "D1M1" "D1M2" "D1M3" "D1M4" ...
  ..$ : NULL
}
\source{Broman,K., et al. (2003) R/qtl: QTL mapping in experimental crosses, Bioinformatics, 19, 889-890.}
\examples{
data(sample.markers)
}
\keyword{datasets}

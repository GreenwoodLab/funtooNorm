% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SampleSet.R
\name{getGRanges}
\alias{getGRanges}
\title{Return a list}
\usage{
getGRanges(object)
}
\arguments{
\item{object}{}
}
\value{
a GRange object of all the methylated positions
}
\description{
Return a list
}
\examples{
require(minfiData)
mySampleSet=fromRGChannelSet(RGsetEx)
getGRanges(mySampleSet)
}

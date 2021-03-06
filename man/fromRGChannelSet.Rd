% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SampleSet.R
\name{fromRGChannelSet}
\alias{fromRGChannelSet}
\title{Creates an object of class SampleSet from a RGChannelSet {minfi}}
\usage{
fromRGChannelSet(myRGChannelSet)
}
\arguments{
\item{myRGChannelSet}{: RGChannelSet, from minfi package, should contain a
cell_type vector in pData}
}
\value{
An object of class 'SampleSet'
}
\description{
Creates a object of class SampleSet from the raw unprocessed
data in RGChannelSet
}
\examples{
require(minfiData)
pData(RGsetEx)$cell_type <- rep(c("type1","type2"),3)
mySampleSet=fromRGChannelSet(RGsetEx)

}

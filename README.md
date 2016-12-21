# funtooNorm 
[![Build Status]
(https://travis-ci.org/GreenwoodLab/funtooNorm.svg?branch=master)]
(https://travis-ci.org/GreenwoodLab/funtooNorm)

The R package ```funtooNorm```  provides a function for normalization of 
Illumina Infinium Human Methylation 450 BeadChip (Illumina 450K) data 
when there are samples from multiple tissues or cell types.

## Installation options
Install from bioconductor  
```shell
source("http://www.bioconductor.org/biocLite.R")  
biocLite("funtooNorm")
```  
Download the current build <a href="https://github.com/GreenwoodLab/funtooNorm/releases" ><b>here</b></a> and install it with
``` shell
$ R CMD INSTALL funtooNorm_0.99.9.tar.gz
```

Or, if you want to build from source, you can also install from GitHub using the [devtools](https://cran.r-project.org/package=devtools)
package in `R`: 
```r
library(devtools)
install_github('GreenwoodLab/funtooNorm', local=TRUE, build_vignettes = TRUE)
```

If you have difficulties with the function ```install_github``` you can download
and install from source using the commands ```R CMD build``` and ```install```:
``` shell
$ git clone https://github.com/GreenwoodLab/funtooNorm.git
$ R CMD build ./funtooNorm
$ R CMD INSTALL funtooNorm_0.99.9.tar.gz
=======
install_github('greenwoodLab/funtooNorm')
```

## Usage

There are two functions in the package, ```funtoonorm``` and ```agreement```. 

The output of ```funtoonorm``` is two matrices: both a normalized methylation data set of beta values and also raw, non-normalized beta values. As an input, the user needs to provide signal A and signal B matrices of data extracted from Illumina IDAT files, as well as control probes signals. The function also requires a list of cell types or tissues. Besides its main purpose of normalization, ```funtoonorm``` can be run in a validation mode and its graphical output is used to choose optimal parameters for normalization.

The function agreement accesses the performance of normalization measuring intra-replicate differences before and after normalization. It takes the output of funtoonorm as an input.

For more details, see the vignette provided with the package or download the pdf file from
<a href="https://github.com/GreenwoodLab/funtooNorm/blob/master/vignettes/funtooNorm.pdf">
<b>here</b></a>.


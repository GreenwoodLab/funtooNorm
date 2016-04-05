# funtooNorm 
[![Build Status]
(https://travis-ci.org/GreenwoodLab/funtooNorm.svg?branch=master)]
(https://travis-ci.org/GreenwoodLab/funtooNorm)

The R package ```funtooNorm```  provides a function for normalization of 
Illumina Infinium Human Methylation 450 BeadChip (Illumina 450K) data 
when there are samples from multiple tissues or cell types.

## Installation options
Download the current build <a href="https://github.com/GreenwoodLab/funtooNorm/releases" ><b>here</b></a> and install it with
``` shell
$ R CMD INSTALL funtooNorm_N.NN.N.tar.gz
```

Or, if you want to build from source, you can also install from GitHub using the [devtools](http://cran.r-project.org/web/packages/devtools/index.html)
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
$ R CMD INSTALL funtooNorm_N.NN.N.tar.gz
=======
install_github('greenwoodLab/funtooNorm')
```

## Usage

There are two functions in the package, ```funtoonorm``` and ```agreement```. 

Please see the vignette provided with the package or download the pdf file from
<a href="https://github.com/RaphaelRaphael/funtooNorm/blob/march2016/vignettes/funtooNorm.pdf">
<b>here</b></a>.


# funtooNorm 


The R package <b>funtooNorm</b>  provides a function for normalization of Illumina Infinium Human Methylation 450
BeadChip (Illumina 450K) data when there are samples from multiple tissues or cell types.

## Installation

<pre>
>library(devtools)
>install_github('adminGreenwoodLab/funtooNorm', local=TRUE, build_vignettes = TRUE)
</pre>
If you have difficulties with direct istallation from github you can download and install source with R CMD build and install commands:
<pre>
>git clone https://github.com/adminGreenwoodLab/funtooNorm.git
>R CMD build ./funtooNorm
>R CMD INSTALL funtooNorm_N.NN.N.tar.gz
</pre>

## Usage

There are two functions in the package, <i>funtoonorm</i> and <i>agreement</i>. 

The output of <i>funtoonorm</i> function is a normalized methylation data set of beta values and raw, not normalized beta values. As an input, user needs to have matrices with signal A and signal B extracted from Illumina IDAT files, as well as control probes signals. The function also requires a list of cell types or tissues. Besides its main purpose of normalization, <i>funtoonorm</i> can be run in validation regime and its graphical output is used to choose optimal parameters for normalization.  

Function <i>agreement</i> accesses the performance of normalization measuring intra-replicate differences before and after normalization. It takes outpit of <i>funtoonorm</i> as an input.

For more details, see the vignette provided in package.


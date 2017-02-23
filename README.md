# RcisTarget
RcisTarget is an R-package to identify transcription factor (TF) binding motifs over-represented on a gene list. 







*Note:* This is a version in development. The main version of the package will be submitted Bioconductor. 

To install this BETA version, you can run the following commands from R:
```
library(devtools)
install_github("aertslab/RcisTarget")

# You might need to install these packages first:
install.packages("devtools", "data.table", "zoo")
source("https://bioconductor.org/biocLite.R"); biocLite("BiocInstaller")
```

A **[tutorial](http://scenic.aertslab.org/tutorials/RcisTarget_tutorial.html)** (vignette) is included in the package.
An HTML version of the tutorial, and the **databases** required to use RcisTarget are available at http://scenic.aertslab.org.

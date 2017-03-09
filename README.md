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

____
Note: From version 0.6+, RcisTarget uses a new format for the ranking databases. As temporary solution, you can convert the current databases with the following code (adapt according to the organism/distance around TSS):

```
library(RcisTarget.hg19.motifDatabases)
data(hg19_10kbpAroundTss_motifRanking) # Addapt database name for 

motifRankings <- rankingWrapper(rankings=hg19_10kbpAroundTss_motifRanking, 
          rowType="gene", colType="motif", org="human", genome="hg19", description="")
```

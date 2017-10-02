# RcisTarget
RcisTarget is an R-package to identify transcription factor (TF) binding motifs over-represented on a gene list. 







*Note:* This is a version in development. The main version of the package will be submitted Bioconductor. 
Current version of **SCENIC** requires [RcisTarget 0.99.0](http://scenic.aertslab.org/downloads/Rpackages/RcisTarget_0.99.0.tar.gz)

To install this BETA version, you can run the following commands from R:
```
devtools::install_github("aertslab/RcisTarget")

# You might need to install these packages first:
install.packages(c("devtools", "data.table", "zoo", "BiocGenerics", "AUCell"))
source("https://bioconductor.org/biocLite.R"); biocLite("BiocInstaller")
```

A **[tutorial](http://scenic.aertslab.org/tutorials/RcisTarget.html)** (vignette) is included in the package.
An HTML version of the tutorial, and the **databases** required to use RcisTarget are available at http://scenic.aertslab.org.

____
Note: From version 0.7+, RcisTarget uses a new format for the ranking databases. As temporary solution, you can convert the previous databases with the following code (adapt according to the organism/distance around TSS):

```
library(RcisTarget.mm9.motifDatabases)
data(mm9_10kbpAroundTss_motifRanking) # Adapt database name

motifRankings <- rankingWrapper(rankings=mm9_10kbpAroundTss_motifRanking, 
          rowType="gene", colType="motif", org="mouse", genome="mm9", maxRank = 5000, description="")
```


# Run the enrichment (or load previous results)
load(paste(file.path(system.file('examples', package='RcisTarget')),
           "motifEnrichmentTable_wGenes.RData", sep="/"))

# Add link to logo
newMotifErnTable <- addLogo(motifEnrichmentTable_wGenes)

# Show table
library(DT)
datatable(newMotifErnTable[,-c("enrichedGenes"), with=FALSE],
          escape = FALSE,
          filter="top",
          options=list(pageLength=5))


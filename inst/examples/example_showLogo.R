# Run the enrichment (or load previous results)
load(paste(file.path(system.file('examples', package='RcisTarget')),
           "motifEnrichmentTable_wGenes.RData", sep="/"))

# Show table as HTML
showLogo(motifEnrichmentTable_wGenes)

#this R file written by /home/con/DNA_Methylation/Public_data/Metilene_public_data/Rnbeads/make_R_files.pl
library(RnBeads)

dataDir <- file.path(getwd())
unlink('1.1_results', recursive = TRUE)
resultDir <- file.path(getwd(),'1.1_results')
warnings_file <- file("pathfinder_warnings.Rout", open="wt")
sink(warnings_file, type="message")

# dataset and file locations
datasetDir <- file.path(getwd())
bedDir <- file.path(datasetDir)
sampleSheet <- file.path(datasetDir, 'sample_annotation1_1.csv')
reportDir <- file.path(resultDir)
################################################################################
# (1) Set analysis options
################################################################################
rnb.options(
	identifiers.column                = "BED",
	import.bed.style                  = "EPP",
	assembly                          = "hg19",
	filtering.low.coverage.masking    = TRUE,
	filtering.greedycut               = FALSE,
	filtering.missing.value.quantile  = 0.5,
	filtering.high.coverage.outliers  = TRUE
)
# optionally disable some parts of the analysis to reduce runtime
rnb.options(
	exploratory.intersample           = FALSE,
	exploratory.region.profiles       = character(0),
	differential.comparison.columns   = c("case"),
	differential.report.sites         = TRUE
)
################################################################################
# (2) Run the analysis
################################################################################
rnb.run.analysis(
	dir.reports=reportDir,
	sample.sheet=sampleSheet,
	data.dir=getwd(),
	data.type="bs.bed.dir"
)

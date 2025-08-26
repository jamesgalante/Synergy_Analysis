
# Combine the hic contact values per chromosome
# Load libraries
suppressPackageStartupMessages(
	library(data.table)
)

combined_data <- rbindlist(lapply(snakemake@input$paired_predictions_with_hic_per_chrom, fread))

fwrite(combined_data, snakemake@output$paired_predictions_with_hic, sep = "\t")

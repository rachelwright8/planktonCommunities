# Data and scripts to analyze eukaryotic plankton communities differences from eight reefs in Bocas del Toro

# Data files
1) licor_data.csv == Light logger data
2) minmeanmax_tempdata.csv == Temperature logger data
3) orig_raw_data_total.mapping.txt == ASV sample information

Get the sequencing files at: https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA507270

# Scripts
1) licorData_code.R == Light data analysis
2) tempData_code.R == Temperature analysis
3) Complete_analysis.R == Plankton community genetic analysis, including analysis using dada2 (analysis through taxonomy assignment), MCMC.OTU (ASV filtering), vegan (PCoA), and DESeq2 (differential abundance testing)

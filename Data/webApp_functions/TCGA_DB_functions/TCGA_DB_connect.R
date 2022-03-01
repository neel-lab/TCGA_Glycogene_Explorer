# This script connects to the TCGA Glycogene Metaanalysis results database
# Tables:
# - patient metadata table: TCGA_patient_metadata
# - TCGA VST count data: TCGA_count_data
# - TCGA differential expression data : TCGA_diffExp_data
# - TCGA GlycoPathway enrichment : TCGA_glycoPathway_enrichments
# - TCGA abbreviations to full names table: TCGA_Diseases

library(RSQLite)
TCGA_DB<-dbConnect(RSQLite::SQLite(),'./Data/TCGA_Glycogene_Metaanalysis.sqlite')

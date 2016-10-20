#!/usr/bin/Rscript
################################################################################
#Toolbox 16S
#
#Parse NCBI taxdump data and store in R-readable format
#
#=> load raw NCBI taxdump files
#=> filter, prune and beautify
#=> store in R-readable format
#
#2016-10-19
#sebastian.schmidt@embl.de
################################################################################

################################################################################
################################################################################
#Preallocate global data structures
PARAM <- list();
PARAM$folder.ncbi_tax <- "/Users/sebastian/Desktop/Arbeit/Projects/Toolbox_16S/reference_databases/NCBI_Taxonomy/"
################################################################################
################################################################################

################################################################################
################################################################################
#To prepare names.dmp file to hold scientific names exclusively, run:
#[cat names.dmp | grep "scientific name" | perl -p -e "s/\t\|//g" > scientific_names.tsv]

#Parse scientific names
tmp.dat <- read.table(paste0(PARAM$folder.ncbi_tax, "scientific_names.tsv"), sep="\t", quote="$$$", comment.char="@")
ncbi.names <- data.frame(Tax_ID=tmp.dat$V1, Scientific_Name=tmp.dat$V2, Unique_Name=tmp.dat$V3)
ncbi.names$Scientific_Name <- as.character(ncbi.names$Scientific_Name)
ncbi.names$Unique_Name <- as.character(ncbi.names$Unique_Name)

#Store
save(ncbi.names, file=paste0(PARAM$folder.ncbi_tax, "ncbi_tax.names.RData"))
################################################################################
################################################################################



q()


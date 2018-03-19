# This file shows how to work with the L1000 data in R

# You need to download 3 L1000 data files from: 
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742
# Count data for 978 landmark genes:
#   GSE92742_Broad_LINCS_Level2_GEX_epsilon_n1269922x978.gctx.gz
# Metadata about each gene:
#   GSE92742_Broad_LINCS_gene_info.txt.gz
# Metadata about each experiment:
#   GSE92742_Broad_LINCS_inst_info.txt.gz

# You need to install the R packages plyr and rhdf5
# install.packages('plyr')
# 
# source('https://bioconductor.org/biocLite.R')
# biocLite('rhdf5')

# Path where downloaded L1000 data is stored and where derived files will be
# written
datapath <- '~/data/L1000'

# Read expression data stored in the L1000 gctx file format
source('utils/cmapR-io.R')
source('utils/cmapR-utils.R')
L1000 <- list(expression=parse.gctx(paste(datapath,
  'GSE92742_Broad_LINCS_Level2_GEX_epsilon_n1269922x978.gctx',
  sep='/'))@mat)

# Read metadata about genes and experiments
load_info <- function(path) {
  read.delim(paste(datapath, path, sep='/'), sep='\t', header=TRUE,
    stringsAsFactors=FALSE)
}
L1000$gene_info <- load_info('GSE92742_Broad_LINCS_gene_info.txt.gz')
L1000$inst_info <- load_info('GSE92742_Broad_LINCS_inst_info.txt.gz')
L1000$inst_info[L1000$inst_info == -666] <- NA

# Reorder inst_info/gene_info to column/row order of L1000
L1000$inst_info <- plyr::join(
  data.frame(inst_id=colnames(L1000$expression)), L1000$inst_info)
L1000$gene_info <- plyr::join(
  data.frame(pr_gene_id=rownames(L1000$expression)), L1000$gene_info)

# Save the L1000 R data object to disk
dir.create(paste(datapath, 'derived', sep='/'),  showWarnings=FALSE)
outfile <- gzfile(paste(datapath, 'derived', 'L1000.rds.gz', sep='/'))
saveRDS(L1000, outfile)
close(outfile)

# Example of how to load and use the L1000 R data object
#
# # load L1000 data
# L1000 <- readRDS('~/data/L1000/derived/L1000.rds.gz')
# 
# # find the indices of experiments with named compounds
# compound_inst_ix <- L1000$inst_info$pert_type == 'trt_cp' &
#   substr(L1000$inst_info$pert_iname, 1, 3) != 'BRD'
# 
# # select the expressions from experiments with named compounds and set easier
# # to read row and column names
# compound_expression <- t(L1000$expression[, compound_inst_ix])
# colnames(compound_expression) <- L1000$gene_info$pr_gene_symbol
# rownames(compound_expression) <- L1000$inst_info$pert_iname[compound_inst_ix]
#
# # info about all columns in gene_info and inst_info can be found at
# # https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU/edit#heading=h.l6bq0r1aih50

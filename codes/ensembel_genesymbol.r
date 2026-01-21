# Load required libraries
library(biomaRt)
library(dplyr)

# Step 1: Load the CSV file containing ENSEMBL IDs and read counts
# Assuming the first column is ENSEMBL IDs and the rest are sample counts
file_path <- "count_matrix.csv"
count_matrix <- read.csv(file_path, row.names = 1)

# Step 2: Use biomaRt to map ENSEMBL IDs to gene symbols (for Homo sapiens)
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Get gene symbols for the ENSEMBL IDs
ensembl_ids <- rownames(count_matrix)
gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   filters = "ensembl_gene_id",
                   values = ensembl_ids,
                   mart = mart)

# Step 3: Replace ENSEMBL IDs with gene symbols in the count matrix
# Merge gene symbols with the original count matrix
count_matrix$ensembl_gene_id <- rownames(count_matrix)
count_matrix_mapped <- merge(gene_info, count_matrix, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")

# Remove rows with missing gene symbols
count_matrix_mapped <- count_matrix_mapped %>% filter(hgnc_symbol != "")

# Step 4: Save the result as a new CSV file with gene symbols instead of ENSEMBL IDs
# Set gene symbols as row names and drop ENSEMBL IDs column
rownames(count_matrix_mapped) <- count_matrix_mapped$hgnc_symbol
count_matrix_mapped <- count_matrix_mapped[, -c(1, 2)]  # Remove ENSEMBL ID and gene symbol columns

# Write the new file to disk
output_file <- "gene_symbol_read_count_matrix.csv"
write.csv(count_matrix_mapped, output_file, row.names = TRUE)

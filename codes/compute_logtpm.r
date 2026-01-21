# Load required libraries
library(biomaRt)
library(DESeq2)

# Step 1: Load the count matrix
# Assuming the CSV file has the format where rows are ENSEMBL IDs and columns are samples
count_matrix <- read.csv("gene_symbol_read_count_matrix.csv", row.names = 1)

# Step 2: Use biomaRt to get gene lengths for your ENSEMBL IDs
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Fetch the gene lengths (in base pairs)
gene_lengths <- getBM(attributes = c("ensembl_gene_id", "transcript_length"),
                      filters = "ensembl_gene_id",
                      values = rownames(count_matrix),
                      mart = mart)

# Ensure that the ENSEMBL IDs match with the rows of your count matrix
gene_lengths <- gene_lengths[match(rownames(count_matrix), gene_lengths$ensembl_gene_id), ]

# Convert gene length to kilobases (as TPM uses kilobases)
gene_lengths$length_kb <- gene_lengths$transcript_length / 1000

# Step 3: Prepare DESeq2 dataset
# For simplicity, assuming you only care about normalization and not differential expression
col_data <- data.frame(row.names = colnames(count_matrix), condition = rep("condition", ncol(count_matrix)))

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~ condition)

# Perform DESeq normalization (for sequencing depth)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# Step 4: Calculate TPM
# Calculate RPK (reads per kilobase)
rpk <- sweep(normalized_counts, 1, gene_lengths$length_kb, FUN = "/")

# Calculate scaling factors (sum of RPKs for each sample)
scaling_factors <- colSums(rpk) / 1e6

# Calculate TPM (RPK divided by scaling factor)
tpm <- sweep(rpk, 2, scaling_factors, FUN = "/")

# Step 5: Apply log(TPM + 1) transformation
log_tpm <- log2(tpm + 1)

# Step 6: Save the result
write.csv(log_tpm, "log_TPM_normalized_data.csv", row.names = TRUE)

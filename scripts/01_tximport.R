# --- START OF USER CONFIGURATION ---
# Please modify the following paths according to your system environment.

# Root directory containing the Kallisto quantification results.
dir_sample <- "/path/to/your/kallistoQuantResultsGENCODE"

# Directory containing the reference genome annotation file (GTF).
dir_reference <- "/path/to/your/Reference"

# Root directory for all analysis outputs from this project.
dir_output <- "/path/to/your/project/output"
# --- END OF USER CONFIGURATION ---

# --- Analysis Script ---

# Load required packages
library(tximport)

# Ensure the output directory exists
if (!dir.exists(dir_output)) {
  dir.create(dir_output, recursive = TRUE)
}
 
# Transcripts need to be associated with gene IDs for gene-level summarization.
# Import the GTF File
gencode.ann <- read.table(file.path(dir_reference, "gencode.v43.annotation.gtf"),
                          header = FALSE,
                          sep = "\t",
                          col.names = c("chromosome name",
                                        "annotation source",
                                        "feature type",
                                        "genomic start location",
                                        "genomic end location",
                                        "score(not used)",
                                        "genomic strand",
                                        "genomic phase (for CDS features)",
                                        "additional information as key-value pairs"))
 
# Find transcripts
gencode.ann.tx <- as.data.frame(gencode.ann[gencode.ann$feature.type == "transcript", ncol(gencode.ann)])
colnames(gencode.ann.tx) <- "additional"
gencode.ann.tx[, 2:3] <- do.call("rbind", strsplit(gencode.ann.tx$additional, "; "))[, 1:2]
colnames(gencode.ann.tx)[2:3] <- c("gene_id", "transcript_id")
gencode.ann.tx$gene_id <- do.call("rbind", strsplit(gencode.ann.tx$gene_id, " "))[, 2]
gencode.ann.tx$transcript_id <- do.call("rbind", strsplit(gencode.ann.tx$transcript_id, " "))[, 2]
 
# Make a data.frame called tx2gene with two columns: 1) transcript ID and 2) gene ID.
tx2gene <- data.frame(transcript_id = gencode.ann.tx$transcript_id,
                      gene_id = gencode.ann.tx$gene_id)
 
# kallisto
files <- file.path(dir_sample, filenames, "abundance.h5")
names(files) <- filenames
all(file.exists(files))
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
txi.kallisto$counts[1:6, 1:6]
 
# kallisto with TSV files
files.tsv <- file.path(dir_sample, filenames, "abundance.tsv")
names(files.tsv) <- filenames
all(file.exists(files.tsv[1]))
txi.kallisto.tsv <- tximport(files.tsv, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
txi.kallisto.tsv$counts[1:6, 1:6]
 
# Save count file
setwd(dir_output)
write.table(txi.kallisto.tsv$counts, file = "txi_kallisto_tsv_count.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
 
# > txi.kallisto.tsv$counts[1:6, 1:6] 
#												A_001    A_002 A_003     A_004    A_005    A_006
# ENSG00000000003.15 20.00000 11.00000    62  7.182721 45.00000 42.00001
# ENSG00000000005.6   0.00000  0.00000     0  0.000000  0.00000  0.00000
# ENSG00000000419.14 48.99945 54.99992    13 19.000000 34.01435 25.00000
# ENSG00000000457.14 25.99996  0.00000     0  0.000000  8.00000  0.00000
# ENSG00000000460.17  0.00000  0.00000     0  5.999998  0.00000  0.00000
# ENSG00000000938.13  0.00000  0.00000     0  0.000000  0.00000  0.00000

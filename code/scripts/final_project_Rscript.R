# Load required libraries
library(dada2)
library(phyloseq)
library(ggplot2)

# Set path and list files
path <- "C:/Users/joaqu/Documents/ee282_project_analysis/16s_seq"
fnFs <- sort(list.files(path, pattern = "READ1-Sequences.*\\.fastq$", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "-READ"), `[`, 1)

# Quality profile and filtering
plotQualityProfile(fnFs[1:4])
filtFs <- file.path(path, "filtered", paste0(sample.names, "-READ1_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, truncLen=210, maxN=0, maxEE=2, truncQ=2, 
                     rm.phix=TRUE, compress=TRUE, multithread=FALSE) 
head(out)

# Error learning, dereplication, and denoising
errF <- learnErrors(filtFs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
names(derepFs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

# Create sequence table and remove chimeras
seqtab <- makeSequenceTable(dadaFs)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim) / sum(seqtab)

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, denoisedF = sapply(dadaFs, getN), nonchim = rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)

# Taxonomy assignment
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/joaqu/Documents/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "C:/Users/joaqu/Documents/silva_species_assignment_v128.fa")
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

# Phyloseq object creation
samples.out <- rownames(seqtab.nochim)
samdf <- data.frame(SampleID = samples.out)
rownames(samdf) <- samples.out
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

# Remove mock samples
ps <- prune_samples(sample_names(ps) != "Mock", ps)

# Rename sample names
new_names <- c(
  "mR406-L1-P184-ACCAGAAATGTC" = "LF-1",
  "mR411-L1-P123-GCGTGCCCGGCC" = "IF",
  "mR411-L1-P171-GTCACGGACATT" = "LF-2",
  "mR411-L1-P306-GAACAAAGAGCG" = "EF-1"
)
sample_names(ps) <- ifelse(sample_names(ps) %in% names(new_names), 
                           new_names[sample_names(ps)], 
                           sample_names(ps))
sample_data(ps)$SampleID <- sample_names(ps)

# Visualizations
plot_richness(ps, x="SampleID", measures=c("Shannon", "Simpson")) +
  labs(title = "Alpha Diversity by Sample")

# Extract and save richness table
richness_data <- estimate_richness(ps, measures = c("Shannon", "Simpson"))
write.csv(richness_data, "C:/Users/joaqu/Documents/ee282_project_analysis/CSV_data/richness_data.csv")
View(richness_data) # View richness data in a table format

# Top 4 families relative abundance (by total family counts)
ps.family <- tax_glom(ps, taxrank="Family")
top4_family_total_counts <- names(sort(taxa_sums(ps.family), decreasing=TRUE))[1:4]
ps.top4_family_total_counts <- prune_taxa(top4_family_total_counts, ps.family)
ps.top4_family_rel <- transform_sample_counts(ps.top4_family_total_counts, function(x) x / sum(x))
plot_bar(ps.top4_family_rel, x="SampleID", fill="Family") + 
  labs(title = "Relative Abundance of Top 4 Families by Total Family Counts")

# Top 4 ASVs relative abundance (by total ASV/OTU counts)
top4 <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:4]
ps.top4 <- transform_sample_counts(ps, function(OTU) OTU / sum(OTU))  # Convert to relative abundance
ps.top4 <- prune_taxa(top4, ps.top4)  # Keep only top 4 taxa
plot_bar(ps.top4, x = "SampleID", fill = "Family") + 
  labs(title = "Relative Abundance of Top 4 Families by Total Counts")

# Extract and save relative abundance for top 4 ASVs
ps.top4_asv_total_counts <- prune_taxa(top4, ps)
top4_asv_df <- as.data.frame(otu_table(ps.top4_asv_total_counts))
write.csv(top4_asv_df, "C:/Users/joaqu/Documents/ee282_project_analysis/CSV_data/top4_asv_relative_abundance_table.csv")
View(top4_asv_df) # View top 4 ASVs relative abundance data in a table format

# General family-level relative abundance
ps.family <- tax_glom(ps, taxrank="Family")
ps.family.rel <- transform_sample_counts(ps.family, function(x) x / sum(x))
plot_bar(ps.family.rel, fill="Family") + 
  theme(axis.text.x = element_text(angle=90, hjust=1)) + 
  labs(title = "Relative Abundance of Taxa by Family")

# Save OTU data as CSVs
otu_data <- as.data.frame(otu_table(ps))
write.csv(otu_data, "C:/Users/joaqu/Documents/ee282_project_analysis/CSV_data/otu_table_raw.csv")
View(otu_data) # View OTU table raw data in a table format

ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))
rel_otu_df <- as.data.frame(otu_table(ps.rel))
write.csv(rel_otu_df, "C:/Users/joaqu/Documents/ee282_project_analysis/CSV_data/otu_table_relative_abundance.csv")
View(rel_otu_df) # View relative abundance data in a table format

percentage_otu_df <- rel_otu_df
percentage_otu_df[] <- percentage_otu_df * 100
write.csv(percentage_otu_df, "C:/Users/joaqu/Documents/ee282_project_analysis/CSV_data/otu_table_percentage.csv")
View(percentage_otu_df) # View percentage abundance data in a table format

# Extract and save relative abundance for top 4 families
ps.top4_family_total_counts <- prune_taxa(top4_family_total_counts, ps.family)
top4_family_otu_df <- as.data.frame(otu_table(ps.top4_family_total_counts))
write.csv(top4_family_otu_df, "C:/Users/joaqu/Documents/ee282_project_analysis/CSV_data/top4_family_relative_abundance_table.csv")
View(top4_family_otu_df) # View top 4 families relative abundance data in a table format

# Extract and save relative abundance for top 4 ASVs
ps.top4 <- prune_taxa(top4, ps)
top4_asv_df <- as.data.frame(otu_table(ps.top4))
write.csv(top4_asv_df, "C:/Users/joaqu/Documents/ee282_project_analysis/CSV_data/top4_asv_relative_abundance_table.csv")
View(top4_asv_df) # View top 4 ASVs relative abundance data in a table format

# Extract and save relative abundance at Family level
family_otu_df <- as.data.frame(otu_table(ps.family.rel))
write.csv(family_otu_df, "C:/Users/joaqu/Documents/ee282_project_analysis/CSV_data/family_level_relative_abundance_table.csv")
View(family_otu_df) # View family-level relative abundance data in a table format

########### table

# Load required libraries
library(ggplot2)
library(gridExtra)  # For tableGrob()
library(grid)       # For grid.draw() and layout

# **1. Read Processing and Filtering Table**
read_processing_data <- data.frame(
  `Sample ID` = c('mR406-L1-P184-ACCAGAAATGTC', 
                  'mR411-L1-P123-GCGTGCCCGGCC', 
                  'mR411-L1-P171-GTCACGGACATT', 
                  'mR411-L1-P306-GAACAAAGAGCG'),
  `Input Reads` = c(25253, 13945, 22645, 44788),
  `Filtered Reads` = c(21898, 8625, 14513, 33572),
  `Denoised Reads` = c(21789, 8525, 14403, 33328),
  `Non-Chimeric Reads` = c(20721, 8522, 13907, 28156)
)

# Create table plot
table1 <- tableGrob(read_processing_data)

# Save the table as a PNG file
png(filename = "C:/Users/joaqu/Documents/ee282_project_analysis/read_processing_and_filtering.png", width = 800, height = 200)
grid.draw(table1)
dev.off()



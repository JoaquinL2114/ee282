#!/bin/bash

# Step 1: Download required files
echo "Downloading genome FASTA file and checksum..."
wget -q https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-chromosome-r6.48.fasta.gz
wget -q https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/md5sum.txt

# Step 2: Verify the integrity of the FASTA file
echo "Verifying file integrity..."
grep "dmel-all-chromosome-r6.48.fasta.gz" md5sum.txt > relevant_md5sum.txt
md5sum -c relevant_md5sum.txt
if [ $? -ne 0 ]; then
    echo "File integrity check failed. Exiting."
    exit 1
fi
echo "File integrity verified."

# Step 3: Unzip the FASTA file
echo "Unzipping genome FASTA file..."
gunzip -k dmel-all-chromosome-r6.48.fasta.gz

# Step 4: Summarize genome partitions and extract data for plots
echo "Processing genome partitions..."
bioawk -c fastx '
    {
        if (length($seq) <= 100000) {
            print length($seq), gc($seq) > "small_seqs_data.txt"
        } else {
            print length($seq), gc($seq) > "large_seqs_data.txt"
        }
    }
' dmel-all-chromosome-r6.48.fasta

# Step 5: Create plots using R
echo "Creating plots..."
Rscript - <<EOF
# Load data
small_data <- read.table("small_seqs_data.txt", col.names = c("Length", "GC"))
large_data <- read.table("large_seqs_data.txt", col.names = c("Length", "GC"))

# Function to save plots
save_plot <- function(data, column, title, xlab, file, log_scale=FALSE, breaks=30) {
  png(file)
  hist(data[[column]], breaks=breaks, main=title, xlab=xlab, col="skyblue", border="white")
  if (log_scale) {
    axis(1, at=log10(pretty(range(data[[column]]))), labels=pretty(range(data[[column]])))
    title(xlab=paste(xlab, "(Log Scale)"))
  }
  dev.off()
}

# Length distribution (histograms)
save_plot(small_data, "Length", "Sequence Length Distribution (≤ 100kb)", "Length", "small_length_hist.png", log_scale=TRUE)
save_plot(large_data, "Length", "Sequence Length Distribution (> 100kb)", "Length", "large_length_hist.png", log_scale=TRUE)

# GC% distribution (histograms)
save_plot(small_data, "GC", "GC% Distribution (≤ 100kb)", "GC%", "small_gc_hist.png", log_scale=FALSE)
save_plot(large_data, "GC", "GC% Distribution (> 100kb)", "GC%", "large_gc_hist.png", log_scale=FALSE)

# Cumulative size (CDF)
plot_cdf <- function(data, title, file) {
  data <- data[order(data$Length, decreasing = TRUE), "Length"]
  cumulative <- cumsum(data)
  png(file)
  plot(cumulative, type="l", col="blue", main=title, xlab="Sequence Index", ylab="Cumulative Size")
  dev.off()
}

plot_cdf(small_data, "Cumulative Size (≤ 100kb)", "small_cdf.png")
plot_cdf(large_data, "Cumulative Size (> 100kb)", "large_cdf.png")
EOF

# Step 6: Assemble genome with hifiasm
echo "Assembling genome with PacBio HiFi reads..."
SRUN_COMMAND="srun --cpus-per-task=16 --mem=64G --time=01:00:00"
READS_FILE="/data/class/ee282/public/ISO1_Hifi_AdaptorRem.40X.fasta.gz"

$SRUN_COMMAND hifiasm -o ISO1_assembly -t 16 $READS_FILE
if [ $? -ne 0 ]; then
    echo "Genome assembly failed. Exiting."
    exit 1
fi
echo "Genome assembly completed."

# Step 7: Convert GFA to FASTA
echo "Converting GFA to FASTA..."
awk '/^S/{print ">"$2"\n"$3}' ISO1_assembly.bp.p_ctg.gfa > ISO1_assembly.bp.p_ctg.fasta
if [ $? -ne 0 ]; then
    echo "GFA to FASTA conversion failed. Exiting."
    exit 1
fi
echo "FASTA conversion completed. Main assembly: ISO1_assembly.bp.p_ctg.fasta"

# Step 8: Calculate N50
echo "Calculating N50..."
bioawk -c fastx '{ print length($seq) }' ISO1_assembly.bp.p_ctg.fasta | sort -nr > contig_lengths.txt
N50=$(awk 'BEGIN { total=0 } { total+=$1; lengths[NR]=$1 } END { half=total/2; cumulative=0; for (i=1; i<=NR; i++) { cumulative+=lengths[i]; if (cumulative>=half) { print lengths[i]; exit } } }' contig_lengths.txt)
echo "N50: $N50"

# Step 9: Create Contiguity Plot
echo "Creating contiguity plot..."
bioawk -c fastx '{ print length($seq) }' dmel-all-chromosome-r6.48.fasta | sort -nr > reference_contigs.txt
if ! command -v plotCDF &> /dev/null; then
    echo "plotCDF command not found. Please ensure it is installed and accessible."
    exit 1
fi
plotCDF contig_lengths.txt reference_contigs.txt > contiguity_comparison.png || {
    echo "plotCDF failed. Exiting."
    exit 1
}

# Step 10: BUSCO Analysis
echo "Running BUSCO analysis on assemblies..."
source $(mamba info --base)/etc/profile.d/mamba.sh
mamba activate busco_env

if ! command -v busco &> /dev/null; then
    echo "BUSCO command not found. Please ensure it is installed in the busco_env."
    mamba deactivate
    exit 1
fi

busco -i ISO1_assembly.bp.p_ctg.fasta -o busco_your_assembly -l diptera_odb10 -m genome -c 16 || {
    echo "BUSCO analysis on your assembly failed. Exiting."
    mamba deactivate
    exit 1
}

busco -i dmel-all-chromosome-r6.48.fasta -o busco_reference -l diptera_odb10 -m genome -c 16 || {
    echo "BUSCO analysis on FlyBase assembly failed. Exiting."
    mamba deactivate
    exit 1
}
mamba deactivate

# Cleanup
echo "Cleaning up temporary files..."
rm dmel-all-chromosome-r6.48.fasta small_seqs_data.txt large_seqs_data.txt relevant_md5sum.txt contig_lengths.txt reference_contigs.txt

echo "Script completed. Assembly results, N50, contiguity plot, and BUSCO results are ready."

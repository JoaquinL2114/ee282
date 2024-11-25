#move to working directory
cd myrepos/ee282/code/scripts

#URLs for genome FASTA files
fasta_urls=(
    "https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-CDS-r6.48.fasta.gz"
    "https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-aligned-r6.48.fasta.gz"
    "https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-chromosome-r6.48.fasta.gz"
    "https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-clones-r6.48.fasta.gz"
    "https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-exon-r6.48.fasta.gz"
    "https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-intergenic-r6.48.fasta.gz"
    "https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-intron-r6.48.fasta.gz"
    "https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-miRNA-r6.48.fasta.gz"
    "https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-miscRNA-r6.48.fasta.gz"
    "https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-ncRNA-r6.48.fasta.gz"
    "https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-predicted-r6.48.fasta.gz"
    "https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-pseudogene-r6.48.fasta.gz"
    "https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-sequence_features-r6.48.fasta.gz"
)

#URL for GTF annotation file
gtf_url="https://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/gtf/dmel-all-r6.48.gtf.gz"
gtf_file="dmel-all-r6.48.gtf.gz"

#make directory to store downloaded files
mkdir -p dmel_genome
cd dmel_genome

#download FASTA files and verify integrity
for url in "${fasta_urls[@]}"; do
    wget "$url"
    echo "Verifying integrity of $(basename $url)"
    sha256sum "$(basename $url)"
done

#download GTF file and verify its integrity
wget "$gtf_url"
echo "Verifying integrity of $(basename $gtf_file)"
sha256sum "$gtf_file"

#initialize variables for total genome summary counts
total_genome_bases=0
total_genome_ns=0
total_genome_sequences=0

#summarize FASTA files
echo "Summarizing FASTA files..."
for fasta_file in *.fasta.gz; do
    echo "Processing $fasta_file..."
    
    #use faSize to get the total bases
    fa_size_output=$(faSize "$fasta_file" | head -n 1)  # Capture only the first line for total bases
    file_bases=$(echo "$fa_size_output" | awk '{print $2}')
    
    #count number of Ns and sequences
    file_ns=$(zcat "$fasta_file" | grep -o 'N' | wc -l)
    file_sequences=$(grep -c "^>" <(zcat "$fasta_file"))

    #output summary for each file
    echo "Summary for $fasta_file:"
    echo "Total bases: ${file_bases:-Unavailable}"
    echo "Total Ns: $file_ns"
    echo "Total sequences: $file_sequences"
    echo "------------------------------------"

    #accumulate totals for entire genome  
    total_genome_bases=$((total_genome_bases + file_bases))
    total_genome_ns=$((total_genome_ns + file_ns))
    total_genome_sequences=$((total_genome_sequences + file_sequences))
done

#summarize GTF file for features and genes
echo "Summarizing GTF annotation..."
bioawk -c gff '{print $3}' "$gtf_file" | sort | uniq -c | sort -nr > feature_summary.txt
bioawk -c gff '{if ($3 == "gene") print $1}' "$gtf_file" | sed -E 's/[^0-9XYLR]//g' | sort | uniq -c > gene_counts_by_chromosome.txt

#display feature summary
echo "Feature Summary (sorted by most common feature type):"
cat feature_summary.txt

#display gene counts per chromosome arm
echo "Gene Counts per Chromosome Arm:"
cat gene_counts_by_chromosome.txt

#compile final report
echo "============================= Final Report ============================="
echo "Total Genome Summary:"
echo "Total bases across all files: $total_genome_bases"
echo "Total Ns across all files: $total_genome_ns"
echo "Total sequences across all files: $total_genome_sequences"
echo "----------------------------------------------------------------------"
echo "Total number of features of each type (sorted from most to least common):"
cat feature_summary.txt
echo "----------------------------------------------------------------------"
echo "Total number of genes per chromosome arm (X, Y, 2L, 2R, 3L, 3R, 4):"
cat gene_counts_by_chromosome.txt
echo "========================================================================"


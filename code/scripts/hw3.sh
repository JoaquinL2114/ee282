#Directory path

GENOME_DIR="/mnt/c/Users/joaqu/Documents/d_melanogaster_genome"

#fasta files list

FASTA_FILES=(
    "dmel-all-aligned-r6.48.fasta"
    "dmel-all-CDS-r6.48.fasta"
    "dmel-all-chromosome-r6.48.fasta"
    "dmel-all-clones-r6.48.fasta"
    "dmel-all-exon-r6.48.fasta"
    "dmel-all-intergenic-r6.48.fasta"
    "dmel-all-intron-r6.48.fasta"
    "dmel-all-miRNA-r6.48.fasta"
    "dmel-all-miscRNA-r6.48.fasta"
    "dmel-all-ncRNA-r6.48.fasta"
    "dmel-all-predicted-r6.48.fasta"
    "dmel-all-pseudogene-r6.48.fasta"
    "dmel-all-sequence_features-r6.48.fasta"
)

#report file

REPORT_FILE="genome_assemblu_summary.txt"
echo "Genome Assembly Summary Report" > "$REPORT_FILE"
echo "=============================" >> "$REPORT_FILE"
echo "" >> "$REPORT_FILE"

#loop each FASTA file and generate summaries

for FILE in "${FASTA_FILES[@]}"; do
    FILE_PATH="$GENOME_DIR/$FILE"

    if [[ -f "$FILE_PATH" ]]; then
        echo "Processing $FILE..."
        echo "File: $FILE" >> "$REPORT_FILE"

        #summary for each file using faSize

        faSize -detailed "$FILE_PATH" > temp_summary.txt

        #extract key details

        total_bases=$(faSize "$FILE_PATH" | grep "total" | awk '{print $2}')
        total_ns=$(faSize "$FILE_PATH" | grep "totalN" | awk '{print $2}')
        total_sequences=$(faSize "$FILE_PATH" | grep "sequences" | awk '{print $2}')

        #report details

        echo "Total nucleotides: $total_bases" >> "$REPORT_FILE"
        echo "Total Ns: $total_ns" >> "$REPORT_FILE"
        echo "Total sequences: $total_sequences" >> "$REPORT_FILE"
        echo "" >> "$REPORT_FILE"

        #clean temp file

        rm temp_summary.txt
    else
        echo "$FILE_PATH not found!"
    fi
done

#display final report

echo "Report generated: $REPORT_FILE"
cat "$REPORT_FILE"


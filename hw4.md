# Homework 4: Pipelines, Plotting, Genome Assembly

## Author: Joaquin Lopez

## Summary Partitions of Genome Assembly

The genome assembly was partitioned into sequences ≤ 100kb and > 100kb using bioawk to process the sequences and calculate:

Sequences ≤ 100kb:

1. Total number of nucleotides: 6178042
2. Total number of Ns: 662593
3. TOtal number of sequences: 1863

Sequences > 100kb:

1. Total number of nucleotides: 137547960
2. Total number of Ns: 490385
3. Total number of sequences: 7

Histograms were created to visualize the distributions of sequence lengths and GC% for both partitions using R:

Plots for all sequences ≤ 100kb and all sequences > 100kb:

### Small Sequences (≤ 100kb)
![Small Sequence Length Distribution](code/scripts/small_length_hist.png)

### Large Sequences (> 100kb)
![Large Sequence Length Distribution](code/scripts/large_length_hist.png)

### GC% Distribution:
#### Small Sequences (≤ 100kb)
![Small GC% Distribution](code/scripts/small_gc_hist.png)

#### Large Sequences (> 100kb)
![Large GC% Distribution](code/scripts/large_gc_hist.png)

### Cumulative Size Plots:
#### Small Sequences (≤ 100kb)
![Cumulative Size (Small Sequences)](code/scripts/small_cdf.png)

#### Large Sequences (> 100kb)
![Cumulative Size (Large Sequences)](code/scripts/large_cdf.png)

## Genome Assembly

The genome was assembled using PacBio HiFi reads and the hifiasm tool. The primary contigs were extracted from the GFA file and converted into FASTA format

### N50 Calculation:

The N50 of the assembly was calculated to measure assembly contiguity:
- N50 Value: 217115751

Contiguity Plot:

A contiguity plot was made comparing the assembled contigs to the reference assembly using plotCDF:




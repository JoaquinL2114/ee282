# Homework 4: Pipelines, Plotting, Genome Assembly

## Author: Joaquin Lopez

## Summary Partitions of Genome Assembly

The genome assembly was partitioned into sequences ≤ 100kb and > 100kb using bioawk toanalyze the genome's composition and structure. Partitioning helps identify characteristics of smaller and larger sequences separately, such as nucleotide counts and sequence complexity:

Sequences ≤ 100kb:

1. Total number of nucleotides: 6178042
2. Total number of Ns: 662593
3. TOtal number of sequences: 1863

Sequences > 100kb:

1. Total number of nucleotides: 137547960
2. Total number of Ns: 490385
3. Total number of sequences: 7

To better understand the genome composition, histograms were generated for sequence length and GC% distributions. These visualizations allow for a quick assessment of sequence variability and GC content.

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

### Cumulative Size Plots

Cumulative size plots visualize the total sequence size sorted from largest to smallest sequences:

#### Small Sequences (≤ 100kb)
![Cumulative Size (Small Sequences)](code/scripts/small_cdf.png)

#### Large Sequences (> 100kb)
![Cumulative Size (Large Sequences)](code/scripts/large_cdf.png)

## Genome Assembly

The genome was assembled using PacBio HiFi reads and the hifiasm tool. TThis process generates contigs from high-quality long reads, which are then analyzed for assembly quality.

### N50 Calculation:

The N50 of the assembly was calculated to measure assembly contiguity:
- N50 Value: 217115751

Contiguity Plot:

A contiguity plot was made comparing the assembled contigs to the reference assembly using plotCDF

![Contiguity Comparison](code/scripts/contiguity_comparison.png)

### BUSCO Analysis:

The BUSCO tool was used to assess the completeness of the genome assemblies against the diptera_odbio lineage dataset. This analysis measures the proption of conserved orthologous genes in the assembly, providing a benchmark for assembly quailty. The BUSCO analysis is currently running for both the assembled genome and the reference genome. Due to the large assembles this process is a taking a significant amount of time.

Results for assembly:

1. Complete BUSCOSs:TBD
2. Fragmented BUSCOS:TBD
3. Missing BUSCOs:TBD

Results for reference assembly:

1. Complete BUSCOs:TBD
2. Fragmented BUSCOs:TBD
3. Missing BUSCOs:TBD





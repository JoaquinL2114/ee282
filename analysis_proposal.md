# Analysis Proposal: 
From Microbes to Phenotype: Characterizing Microbial Diversity and Community Differences Between Early and Late Flowering _Arabidopsis thaliana_

Author: Joaquin Lopez, Alejandra Hernandez Teran, Maria Rebolleda-Gomez

## Introduction

Microbial communities influence plant phenotypes by shaping ontogeny through processes such as nutrient acquisition and stress response, in turn affecting life history traits such as growth rate, reproductive timing, and overall fitness (Metcalf et. al 2019). A growing body of research has introduced the concept of microbiome-dependent ontogeny timing (MiDOT), where host species rely on microbial cues to regulate critical life history transitions (Metcalf et. al 2019). MiDOT also includes both accelerated and delayed developmental events  in response  to microbial communities, highlighting relationships between hosts and their microbial consortia (Metcalf et. al 2019). Previous research studies have shown that plant-associated microbiomes can influence traits such as disease resistance, growth, and abiotic stress tolerance across multiple generations (Puy et al., 2022; Lu et al. 2018), showing MiDOT can also play an important role in these processes. For instance, a recent experiment using multi-generational approaches revealed that microbiomes selected for specific plant traits, like flowering time, can reproducibly alter plant development and function, showing how microbial communities drive rapid phenotypic changes in plant populations over time (Lu et. al 2018). Highlighting the role of microbiota in regulating key life history transitions—such as flowering in plants—by either accelerating or delaying ontogenetic events, allowing hosts to adapt developmentally to microbial cues that ultimately maximize their fitness.

## Methods

Building on these findings, we designed an experimental study aimed at exploring how MiDOT shapes plant development and function, influencing various aspects of life history. This study began with the collection of rhizosphere samples from _Brassica nigra_ (_B. nigra_) found across the University of California, Irvine (UCI) campus. Soil slurries were prepared to inoculate Col-0 genotype _Arabidopsis thaliana_ (_A. thaliana_) Col-0 mutant seed types. After inoculation, _A. thaliana_ seedlings were grown within a germination chamber using Murashige and Skoog (MS) media lacking sucrose. Once true leaves emerged, the seedlings were transferred into falcon tubes containing sterile sand and slow release fertilizer. Phenotypic traits, including flowering time, were measured every other day. 

Based on results, microbial communities promoting early and late flowering were selected for an experimental evolution trial. In this trial, we used selected inoculates to grow 48 replicates of Col-0 _A. thaliana_ seedlings per treatment, including early flowering (EF) and late flowering (LF) groups, as well as controls. Once seedlings from transfer 1 germinated, they were transferred to soil, and phenotypic traits were measured routinely. Rhizosphere samples were then taken from EF and LF plant groups to inoculate a new batch of Col-0 _A. thaliana_ seedlings for transfer 2. This process continued until 4 transfers were completed. Additionally, rhizosphere samples from transfer 1 rhizosphere were used to extract microbial DNA using a soil DNA extraction kit, followed by the amplification of the 16s rRNA using PCR. Illumina next generation sequencing (NGS) was then used to perform high throughput sequencing and obtain 16s rRNA sequencing data. 

## Topic Proposal

Using this data I aim to explore the role of microbial community composition in regulating flowering time in _A. thaliana_. Specifically, I will investigate how shifts in microbial community composition and diversity influence the timing of flowering by analyzing 16s rRNA sequences extracted from the rhizosphere of microbial communities associated with early and late flowering plant types. Key questions guiding this project include: How does the microbial community composition differ between early and late stage flowering types? Are there specific taxa that exhibit higher abundance in one community compared to the other? Additionally, I seek to understand how the alpha and beta diversity of these microbial communities vary between early and late flowering types. I aim to assess how diverse the microbial populations are within each community, and how these differences contribute to the developmental outcomes observed in _A. thaliana_. 

## Analysis Proposal

To analyze and visualize the taxonomic composition and diversity of EF and LF microbial communities, I will be using QIIME 2, an open-source platform for microbiome data analysis. The workflow will begin by importing the FASTQ files, followed by demultiplexing and quality filtering with DADA2 to generate a feature table and representative sequences. Next, a pre-trained classifier, such as Greengenes, will be used to assign taxonomy, allowing for bar plots to be used to visualize the taxonomic composition across EF and LF groups. For alpha diversity analysis, Shannon's diversity index will assess within-sample richness, and a Kruskal-Wallis test will evaluate significant differences in richness between EF and LF microbial communities. In beta diversity analysis, Bray-Curtis dissimilarity will measure the differences between microbial communities, and a principal coordinates analysis (PCoA) plot will allow for the visualization of these differences. Additionally, PERMANOVA will test whether community compositions differ significantly between EF and LF groups. 

## Feasibility

This is a feasible project as I already have the FASTQ files of the 16s rRNA sequencing data and will be able to effectively conduct the taxonomic analysis as I continue to learn how to use mamba and R throughout the course. QIMME2 additionally has a tutorial named “Moving Pictures” to help guide me in this process. However, if I have difficulties using QIIME2 in mamba, I could also use DADA2 in R as they have a similar workflow.

## References

Metcalf CJE, Henry LP, Rebolleda-Gómez M, Koskella B.2019.Why Evolve Reliance on the Microbiome for Timing of Ontogeny?. mBio10:10.1128/mbio.01496-19.https://doi.org/10.1128/mbio.01496-19

Puy, J., Carmona, C. P., Hiiesalu, I., Öpik, M., de Bello, F., & Moora, M. (2022). Mycorrhizal symbiosis alleviates plant water deficit within and across generations via phenotypic 
plasticity. Journal of Ecology, 110, 262–276. https://doi.org/10.1111/1365-2745.13810 

Lu, T., Ke, M., Lavoie, M. et al. Rhizosphere microorganisms can influence the timing of plant flowering. Microbiome 6, 231 (2018). https://doi.org/10.1186/s40168-018-0615-0



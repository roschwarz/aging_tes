# Principles_of_tissue-specific_regulation_of_transposable_elements_during_aging

Code for the Manuscript "Principles of tissue-specific regulation of transposable elements during aging".

# Annotations

## Create TE annotation as an R object

```{bash}
Rscript src/annotation/build_TE_annotation.R 
```

The TE annotation in R does have more entries than the bed file from RepeatMasker. The reason for that is
that one TE can overlap with multiple genes. The build TE annotation script stores also the genomic 
context of TEs, e.g. exonic, intronic, and intergenic.

## Build TE Island annotation 

```{bash}

bash src/annotation/build_TE_island_annotation.sh \
    data/shared/mm10_transposome_sorted.bed \  # te instance annotation from repeatmasker
    data/shared/GRCm38.p6.chrom.sizes \ # Chromosome sizes
    data/processed/annotation   \ # output directory

```

# Functions

## Overlap of features

process_overlapping
Purpose:
This generic function finds overlaps between two sets of genomic ranges (e.g., gene TSS and transposable 
elements), flexibly summarizes overlap statistics, and lets you annotate overlapping intervals with relevant 
metadata.

What you get:

A list containing:

query_with_hit: All query features (as GRanges) that overlap with any subject feature.

subject_hits: The matching subject features (as GRanges), with selected metadata columns optionally 
transferred from queries.

stats: A summary list with

the total number of query features

the number of unique entities with at least one overlap (based on the identifier column if specified)

and the percentage these represent of all queries.

How to interpret:

Run the function to get which (and how many) query features overlap with your subject features.
The id_col argument lets you supply an identifier (e.g., gene ID) to count unique biological entities with 
overlaps, rather than just counting ranges. The extend_bp argument lets you define a symmetric window around 
each query feature—helpful for less precise features (e.g., ±200bp around a TSS). If you want to bring 
information from the query over to the subject (e.g., gene ID or score), provide the relevant columns via 
col_query_to_subj. The function prints the core overlap statistics, and returns both overlap pairs and a 
stats summary for downstream aggregation, table construction, QC, or plotting.

Typical use case:

Overlap annotated TSS with TE loci to determine what fraction of genes have TE-associated transcription 
start sites, and retrieve for each hit both TSS and TE information.

Example interpretation:

If you get pct_with_hit = 33.3%, then 1/3 of your unique query features are associated with at least one 
subject interval under your chosen window.

Switching the order of query/subject reverses the direction (“of TSS, what fraction overlap TE” vs. “of TEs, 
what fraction are near a TSS”).

Input requirements:

Both arguments must be GRanges objects; identifier columns (for id_col/col_query_to_subj) must exist in your 
query’s mcols.

Flexible integration:

Output structure and column-forwarding enable you to chain this result directly into summarization, 
annotation, and plotting workflows for genomics research.
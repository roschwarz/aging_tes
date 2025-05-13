# Principles_of_tissue-specific_regulation_of_transposable_elements_during_aging

Code for the Manuscript "Principles of tissue-specific regulation of transposable elements during aging".

# Create TE annotation as an R object

```{bash}
Rscript src/annotation/build_TE_annotation.R 
```

Stores also the genomic context of TEs, e.g. exonic, intronic, and intergenic.

# Build TE Island annotation 

```{bash}

bash src/annotation/build_TE_island_annotation.sh \
    data/shared/mm10_transposome_sorted.bed \  # te instance annotation from repeatmasker
    data/shared/GRCm38.p6.chrom.sizes \ # Chromosome sizes
    data/processed/annotation   \ # output directory

```

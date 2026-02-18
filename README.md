# Classification of Disease-Associated Genomic Variants in Rice Crops using Machine Learning.


## Team Members

- Aditya R Auradkar - BL.SC.U4AIE24005
- Bhuvaneswar Reddy - BL.SC.U4AIE24019
- P Varun Sathwik - BL.SC.U4AIE24035




## Project Overview

This project focuses on the systematic curation and validation of a genomic dataset for analyzing disease-associated nucleotide variation in rice (*Oryza sativa, Indica group*). The dataset consists of full-length healthy gene sequences and corresponding variant-bearing sequences constructed using biologically realistic mutation modeling.

The primary objective of this work is to build a balanced, validated, and structurally consistent genomic dataset suitable for downstream computational analysis.



# Dataset Curation

## Healthy Gene Sequences

- 3000 full-length rice gene sequences
- Average sequence length: ~7833 base pairs
- Stored in FASTA format
- Represent baseline (healthy) genomic references



## Diseased Gene Sequences

- 3000 variant-bearing sequences derived from healthy genes
- Mutation burden: 8–20 nucleotide substitutions per gene
- Average mutations per gene: ~13.9
- Sequence length preserved

Mutation characteristics:

- 70% transition mutations (A↔G, C↔T)
- 30% transversion mutations
- No insertions or deletions introduced
- Mutations distributed across valid nucleotide positions

The mutation model reflects biologically realistic nucleotide substitution behavior observed in plant genomes.



## Dataset Summary

| Category | Count |
|----------|-------|
| Healthy Sequences | 3000 |
| Diseased Sequences | 3000 |
| Total Sequences | 6000 |
| Average Length | ~7833 bp |



# Redundancy Analysis

A redundancy check was performed across all sequences to ensure uniqueness and prevent dataset bias.

Results:

Total Samples: 6000  
Unique Sequences: 5996  
Duplicate Sequences: 4  

Duplicate rate ≈ 0.067%

This minimal redundancy confirms high dataset integrity.



# Alignment Verification

Pairwise sequence alignment was conducted between representative healthy and diseased sequences to confirm mutation incorporation.

Sequence Length: 7833  
Matching Positions: 7819  
Differing Positions: 14  

This confirms:

- Correct mutation integration
- Preserved sequence length
- No structural distortion

The full alignment output is available in:

alignment_example.txt



# Feature Extraction

Each sequence was transformed into numerical features for computational analysis.

## 3-mer Frequency Representation

- k = 3
- Total possible 3-mers = 64
- Normalized frequency per sequence

## GC Content

- Proportion of G and C nucleotides per sequence

Total Features per Sequence:

- 64 k-mer features  
- 1 GC content feature  
- Total = 65 features  



#  Generated Files

- feature_extraction.py
- rice_healthy_3000.fasta
- rice_diseased_3000.fasta
- rice_full_dataset.csv
- train_features.csv
- test_features.csv
- alignment_example.txt

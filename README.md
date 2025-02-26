# Genome Assembly

## Overview

This repository contains a Python implementation of a genome assembly algorithm. Genome assembly involves reconstructing the original genome sequence from short DNA fragments (reads) obtained through sequencing technologies.

## Files in the Repository

- **Genome Assembly.py**: The main Python script that performs the genome assembly process.
- **genome_assembly.txt**: A sample input file containing DNA reads to test the assembly algorithm.

## How the Genome Assembly Algorithm Works

The algorithm follows these steps:

1. **Input Reads**: Load DNA sequences (reads) from the input file.
2. **Overlap Detection**: Identify overlaps between the suffix of one read and the prefix of another.
3. **Construct Overlap Matrix**: Build a matrix that quantifies the overlap between each pair of sequences.
4. **Sequence Ordering**: Determine the best order to merge sequences based on overlap scores.
5. **Sequence Merging**: Merge reads in the correct order to construct the contiguous genome sequence.

### Example

Given the following sequences:

```
1 : GGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTCGTCCAGACCCCTAGC
3 : GTCTTCAGTAGAAAATTGTTTTTTTCTTCCAAGAGGTCGGAGTCGTGAACACATCAGT
2 : CTTTACCCGGAAGAGCGGGACGCTGCCCTGCGCGATTCCAGGCTCCCCACGGG
5 : CGATTCCAGGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTC
4 : TGCGAGGGAAGTGAAGTATTTGACCCTTTACCCGGAAGAGCG
6 : TGACAGTAGATCTCGTCCAGACCCCTAGCTGGTACGTCTTCAGTAGAAAATTGTTTTTTTCTTCCAAGAGGTCGGAGT
```

#### Step 1: Construct the Overlap Matrix

Each cell in the matrix represents the length of the overlap between the suffix of one sequence and the prefix of another:

```
    1   2   3   4   5   6
1   -   1   0   0   1  29
2  13   -   1   0  21   0
3   0   0   -   1   0   1
4   1  17   1   -   2   0
5  39   1   0   0   -  14
6   0   0  43   1   0   -
```

#### Step 2: Determine the Best Order to Merge Sequences

Based on overlap scores, the best order to assemble sequences is:

```
4 → 2 → 5 → 1 → 6 → 3
```

#### Step 3: Merge Sequences

By merging sequences in this order while considering the overlaps, we get the final assembled genome:

```
TGCGAGGGAAGTGAAGTATTTGACCCTTTACCCGGAAGAGCGGGACGCTGCCCTGCGCGATTCCAGGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTCGTCCAGACCCCTAGCTGGTACGTCTTCAGTAGAAAATTGTTTTTTTCTTCCAAGAGGTCGGAGTCGTGAACACATCAGT
```

## Limitations of This Technique

- **Error Sensitivity**: The approach assumes that sequencing errors are minimal. Errors can lead to incorrect overlaps and misassemblies.
- **Short Read Lengths**: If the reads are too short, the overlap detection may be insufficient for proper assembly.
- **Repetitive Sequences**: The algorithm struggles with highly repetitive regions, which can cause misassemblies due to multiple valid overlaps.

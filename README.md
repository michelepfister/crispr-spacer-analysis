# CRISPR Spacer Analysis Pipeline (Phytobacter)

## Overview
This project implements a Python-based pipeline for analyzing CRISPR spacer sequences from bacterial isolates. It was developed as part of my Bachelor thesis, focusing on genomic characterization of *Phytobacter* isolates and their CRISPR spacer diversity.

The pipeline parses spacer sequences from FASTA files, groups similar spacers based on sequence similarity, and generates presence–absence matrices for comparative analysis across isolates.

---

## Key Features
- FASTA parsing of CRISPR spacer sequences
- Grouping of similar spacers using Hamming distance (≤ 2 mismatches)
- Identification of unique and shared spacers across isolates
- Generation of:
  - Spacer definition tables
  - Presence–absence matrices
- Modular pipeline design (CLI + processing + parsing)

---

## Workflow

1. **Input**
   - FASTA file containing CRISPR spacer sequences from multiple isolates

2. **Processing**
   - Parse sequences into structured format
   - Compare spacers using Hamming distance
   - Group similar spacers into clusters

3. **Output**
   - CSV file with spacer clusters
   - Presence/absence matrix across isolates

---

## Example Output

### Spacer Definition Table
| ID   | Sequences                  |
|------|---------------------------|
| S001 | ACGTACC; ACGTGCC          |
| S002 | CGGCTTA                   |

### Presence/Absence Matrix
| Sample   | S001 | S002 |
|----------|------|------|
| Isolate1 | x    | -    |
| Isolate2 | x    | x    |

---

## Project Structure
main.py      # CLI entry point
process.py   # pipeline orchestration
parser.py    # parsing and sequence analysis logic

---

## Methods

- **Sequence comparison**: Hamming distance
- **Clustering threshold**: ≤ 2 mismatches
- **Data format**: FASTA → structured Python objects → CSV

---

## Technologies
- Python
- Bioinformatics data processing
- Sequence analysis

---

## Context (Bachelor Thesis)

This pipeline was developed in the context of my Bachelor thesis on the genomic analysis of *Phytobacter*, with a focus on CRISPR spacer diversity and strain comparison.

---

## Michèle Pfister
MSc Applied Computational Life Sciences  
Bioinformatics | Genomics | Health Data Science

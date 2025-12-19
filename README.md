# DNA/RNA Sequence Analyzer

DNA/RNA Sequence Analyzer is a C++ console-based bioinformatics tool that analyzes nucleotide sequences from direct user input or FASTA files. It automatically detects whether the sequence is DNA or RNA, calculates GC content, and identifies Open Reading Frames (ORFs) across all six reading frames (three forward and three reverse complement).

---

## Features
- Supports both DNA and RNA sequences
- Automatic RNA detection and conversion to DNA
- GC content calculation with base count
- ORF detection in six reading frames
- FASTA file input support
- Console output or save results to a text file
- Object-Oriented Programming (OOP) based design

---

## Input Options
1. Direct sequence input via console  
2. FASTA file input  

---

## Output Details
- Detected sequence type (DNA or RNA)
- Base composition (A, T/U, G, C)
- GC percentage
- ORF information:
  - Start position
  - End position
  - Length
  - Nucleotide sequence

---

## How to Compile and Run

### Compile
```bash
g++ main.cpp -o analyzer

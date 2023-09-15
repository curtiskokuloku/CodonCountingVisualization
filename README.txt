Codon Counting and Visualization

Author: Curtis Kin Kokuloku
x500: kokul003

Description:

    - This program counts codon usage in DNA sequences, converts them to amino acid counts, and visualizes the results using a barplot.
    - The coding sequences show the correct reading frame codon usage. The whole genome contains random frame shifts so is not an accurate reflection of codon usage.

Usage:
The program can be run in two modes:

    - To generate a single codon count CSV file:
    
        python count_codons.py <input.fna> <output.csv>

    - To run the full workflow of generating codon and amino acid counts and visualization:

        python count_codons.py

Prerequisites:
The program requires Python 3 and the following packages:

    - Pandas
    - Matplotlib

Install required packages using the following command in the command-line/terminal:

    pip install pandas matplotlib

The key steps are:

    - Count codon usage in input FASTA files and general CSV files with the codon counts
    - Convert codon counts to amino acid counts and generate CSV files containing amino acid counts
    - Generate visualizations comparing coding sequences vs whole genomes (for both codon counts and amino acid counts)

I certify that the information contained in this README file is complete and accurate. I have both read and followed
the course policies in the ‘Academic Integrity - Course Policy’ section of the course syllabus.
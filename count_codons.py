import matplotlib.pyplot as plt                                     # For plotting the barplots
import pandas as pd                                                 # For table manipulation
import logging                                                      # For saving progress, modifications, and debugging information
import sys                                                          # For parsing command line arguments

# Configure logging
logging.basicConfig(filename='logs.txt', level=logging.INFO)
LOGGER = logging.getLogger(__name__)

# Codon to Amino Acid conversion table
CODON_TO_AMINO_TABLE = {
    #    T               C               A               G
    "TTT": "Phe",   "TTC": "Phe",   "TTA": "Leu",   "TTG": "Leu",
    "TCT": "Ser",   "TCC": "Ser",   "TCA": "Ser",   "TCG": "Cys",
    "TAT": "Tyr",   "TAC": "Tyr",   "TAA": "stop",  "TAG": "stop",
    "TGT": "Cys",   "TGC": "Cys",   "TGA": "stop",  "TGG": "Trp",
    "CTT": "Leu",   "CTC": "Leu",   "CTA": "Leu",   "CTG": "Leu",
    "CCT": "Phe",   "CCC": "Phe",   "CCA": "Phe",   "CCG": "Phe",
    "CAT": "His",   "CAC": "His",   "CAA": "Gln",   "CAG": "Gln",
    "CGT": "Arg",   "CGC": "Arg",   "CGA": "Arg",   "CGG": "Arg",
    "ATT": "Ile",   "ATC": "Ile",   "ATA": "Ile",   "ATG": "Met",
    "ACT": "Thr",   "ACC": "Thr",   "ACA": "Thr",   "ACG": "Thr",
    "AAT": "Asn",   "AAC": "Asn",   "AAA": "Lys",   "AAG": "Lys",
    "AGT": "Ser",   "AGC": "Ser",   "AGA": "Arg",   "AGG": "Agr",
    "GTT": "Val",   "GTC": "Val",   "GTA": "Val",   "GTG": "Val",
    "GCT": "Ala",   "GCC": "Ala",   "GCA": "Ala",   "GCG": "Ala",
    "GAT": "Asp",   "GAC": "Asp",   "GAA": "Glu",   "GAG": "Glu",
    "GGT": "Gly",   "GGC": "Gly",   "GGA": "Gly",   "GGG": "Gly"
}


def count_codons(fasta_file, csv_file):
    '''
    Purpose:
        Count how many times each 3-character codon appears in the given file
    Parameter(s):
        fasta_file - The given file being read from
        csv_file - The file to be written to
    Return Value:
        A dictionary representing the codons and their corresponding counts
    '''
    LOGGER.info(f'Counting codons in {fasta_file}')
    codon_counts = {}                                               # Initialize codon counter

    with open(fasta_file) as f:                                     # Open input file
        for line in f:                                              # Read sequences
            line = line.rstrip()
            if line.startswith('>'):                                # Lines that start with '>' are header lines that contain metadata about the sequence.
                continue
            if (len(line)-2) % 3 != 0:                              # Accounts for lines with length not evenly divisible by 3
                line = line[:len(line)-1]
            for i in range(0, len(line)-2, 3):                      # Iterate over sequence by codons
                codon = line[i:i+3]
                if codon in codon_counts:                           # Update codon count
                    codon_counts[codon] += 1
                else:
                    codon_counts[codon] = 1

    with open(csv_file, 'w') as f:                                  # Write counts to output CSV file
        f.write("Codon,Count\n")                                    # Write the header
        for codon, count in codon_counts.items():
            f.write(f"{codon},{count}\n")                           # Write codon and count
        
    LOGGER.info(f'Finished counting, saving results to {csv_file}')
    
    return codon_counts


def convert_to_amino(codon_counts, output_csv_file):
    '''
    Purpose:
        Convert codon counts to amino acid counts and save them to a CSV file.
    Parameters:
        codon_counts - dictionary containing codon counts
        output_csv_file - file to save the amino acid counts in CSV format
    Return Value:
        None
    '''
    LOGGER.info(f'Counting amino acid in Codon Counts Dictionary')
    amino_acid_counts = {}
    
    for codon, count in codon_counts.items():
        amino_acid = CODON_TO_AMINO_TABLE.get(codon, "Unknown")                                             # Retrieves codon value from conversion table; Returns "Unknown" if not found
        amino_acid_counts[amino_acid] = amino_acid_counts.get(amino_acid, 0) + count                        # Updates dictionary with the count of the amino acid.
        
    amino_acid_counts_df = pd.DataFrame(list(amino_acid_counts.items()), columns=["Amino Acid", "Count"])   # Create a DataFrame from the amino acid counts
    
    amino_acid_counts_df.to_csv(output_csv_file, index=False)                                               # Save the amino acid counts to a CSV file
    
    LOGGER.info(f'Saved amino acid counts to {output_csv_file}')


def plot_codons(gene_codons_csv, genome_codons_csv, output_png_file):
    '''
    Purpose:
        Generate a barplot comparing side-by-side the counts of each codon in the two different files
    Parameter(s):
        gene_codons_csv - input codon counts csv file for separate genes
        genome_codons_csv - input codon counts csv file for whole genome
        output_png_file - output png file for barplot
    Return Value:
        None
    '''
    # DataFrame of codon counts for separate genes and whole genome
    gene_codons_csv = pd.read_csv(gene_codons_csv)
    genome_codons_csv = pd.read_csv(genome_codons_csv) 
    
    LOGGER.info(f'Loading CSV files: {gene_codons_csv}, {genome_codons_csv}')
    
    # Sort gene counts and genome counts in descending order by 'Count' column
    gene_codons_csv = gene_codons_csv.sort_values('Count', ascending=False)
    genome_codons_csv = genome_codons_csv.sort_values('Count', ascending=False)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Adds a label for this set of bars, and Plots a bar chart using the 'Codon' column as x value and 'Frequency' as height.
    ax.bar(gene_codons_csv['Codon'], gene_codons_csv['Count'], label='Coding Sequences (correct frame shift)')
    ax.bar(genome_codons_csv['Codon'], genome_codons_csv['Count'], label='Whole Genome (random frame shift)')
    ax.set_xlabel('Codon')
    ax.set_ylabel('Frequency')
    ax.set_title('Codon Counts')
    ax.legend()
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=90)
    
    # Automatically adjusts subplot parameters for better spacing, and saves the figure to provided file path
    plt.tight_layout()
    plt.savefig(output_png_file)
    plt.close()
    
    # Save progress to logger file
    LOGGER.info(f'Saved plot to {output_png_file}')
    
def plot_amino_acid(gene_amino_acid_csv, genome_amino_acid_csv, output_png_file):
    '''
    Purpose:
        Generate a barplot comparing side-by-side the counts of each amino acid in the two different files
    Parameter(s):
        gene_amino_csv_file - CSV file containing amino acid counts for separate genes
        genome_amino_csv_file - CSV file containing amino acid counts for the whole genome
        output_png_file - output png file for barplot
    Return Value
        None
    '''
    gene_amino_acid_csv = pd.read_csv(gene_amino_acid_csv)
    genome_amino_acid_csv = pd.read_csv(genome_amino_acid_csv)
    
    LOGGER.info(f'Loading CSV files: {gene_amino_acid_csv}, {genome_amino_acid_csv}')
    
    gene_amino_acid_csv = gene_amino_acid_csv.sort_values('Count', ascending=False)
    genome_amino_acid_csv = genome_amino_acid_csv.sort_values('Count', ascending=False)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    ax.bar(gene_amino_acid_csv['Amino Acid'], gene_amino_acid_csv['Count'], label='Coding Sequences (correct frame shift)')
    ax.bar(genome_amino_acid_csv['Amino Acid'], genome_amino_acid_csv['Count'], label='Whole Genome (random frame shift)')
    ax.set_xlabel('Amino Acid')
    ax.set_ylabel('Frequency')
    ax.set_title('Amino Acid Counts')
    ax.legend()
    
    plt.xticks(rotation=90)
    
    plt.tight_layout()
    plt.savefig(output_png_file)
    plt.close()
    
    LOGGER.info(f'Saved amino acid plot to {output_png_file}')

def test():
    '''Test Driver for the program'''
    LOGGER.info('Running the program on a small fake genome file...')
    test_fasta = "test_genome.fna"
    test_genome_csv = "test_genome.csv"
    test_codon_counts = count_codons(test_fasta, test_genome_csv)
    

def main():
    '''Run entire program to generate codon and amino acid counts csv files.
       Also generate a barplot for both genome counts and amino acid counts.'''
    
    # Generate gene and genome codon CSVs
    LOGGER.info('Generating CSV files...')
    gene_fasta = "SARS-CoV-2_separate_genes.fna"
    gene_csv = "separate_genes.csv"
    genome_fasta = "SARS-CoV-2_whole_genome.fna"
    genome_csv = "whole_genome.csv"
    gene_codons = count_codons(gene_fasta, gene_csv)
    genome_codons = count_codons(genome_fasta, genome_csv)
    
    # Codon plot
    LOGGER.info('Generating codon usage barplot...')
    codon_png = "codon_counts.png"
    plot_codons(gene_csv, genome_csv, codon_png)
    
    # Convert codons counts to amino acid counts
    LOGGER.info('Converting codons to amino acid...')
    LOGGER.info('Generating amino acid CSV files...')
    gene_amino_csv = "gene_amino_counts.csv"
    genome_amino_csv = "genome_amino_counts.csv"
    convert_to_amino(gene_codons, gene_amino_csv)
    convert_to_amino(genome_codons, genome_amino_csv)
    
    # Amino acid plot
    LOGGER.info('Generating amino acid usage barplot...')
    amino_png = "amino_acids_counts.png"
    plot_amino_acid(gene_amino_csv, genome_amino_csv, amino_png)

       
if __name__ == '__main__':
    try:
        if len(sys.argv) == 3:
            # Codon CSV generation
            LOGGER.info('Generating Codon counts CSV file...')
            count_codons = count_codons(sys.argv[1], sys.argv[2])
        else:
            LOGGER.info('Running the entire program...')
            # test()
            main()
    except Exception as e:
        LOGGER.exception('Error running program')
        
        
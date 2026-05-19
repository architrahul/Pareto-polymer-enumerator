# File: process_hilbert_basis.py
#
# This script processes two input files: 'monomers.txt' and 'hilbert_basis.txt'.
# It reads the contents of 'monomers.txt', applies specific processing rules,
# and uses the matrix in 'hilbert_basis.txt' to generate output based on the following rules:
# 1. Each row in 'hilbert_basis.txt' represents a multiset of lines from 'monomers.txt'.
# 2. The script outputs the processed lines from 'monomers.txt' according to the counts
#    in each row of 'hilbert_basis.txt'.
# 3. If there are more columns in 'hilbert_basis.txt' than lines in 'monomers.txt',
#    the extra columns are ignored.
#
# Processing rules for 'monomers.txt':
# - If a line starts with text followed by a colon, that part (including the colon) is removed.
# - If a line ends with a comma followed by text, the comma and everything after it is removed.
#
# The output format is:
# nx [processed_line_content]
# Where 'n' is the count from 'hilbert_basis.txt' and 'processed_line_content' is from 'monomers.txt'.
# Each multiset is separated by a blank line.

import re

def process_monomer_line(line):
    """Process a single line from the monomers file."""
    # Remove leading text and colon
    line = re.sub(r'^[^:]*:\s*', '', line)
    # Remove trailing comma and text
    line = re.sub(r',.*$', '', line)
    return line.strip()

def read_monomers_file(filename):
    """Read and process the contents of the monomers file."""
    with open(filename, 'r') as file:
        return [process_monomer_line(line) for line in file]

def read_hilbert_basis(filename):
    """Read the Hilbert basis matrix from the file."""
    with open(filename, 'r') as file:
        return [list(map(int, line.split())) for line in file]

def process_files(monomers_filename, hilbert_basis_filename):
    """Process the monomers and Hilbert basis files and generate output."""
    monomer_lines = read_monomers_file(monomers_filename)
    hilbert_basis = read_hilbert_basis(hilbert_basis_filename)
    
    output = []
    for row in hilbert_basis:
        multiset = []
        for i, count in enumerate(row[:len(monomer_lines)]):
            if count > 0:
                multiset.append(f"{count}x {monomer_lines[i]}")
        if multiset:
            output.append("\n".join(multiset))
    
    return "\n\n".join(output)

def main():
    monomers_filename = 'monomers.txt'
    hilbert_basis_filename = 'hilbert_basis.txt'
    
    try:
        result = process_files(monomers_filename, hilbert_basis_filename)
        print(result)
    except FileNotFoundError as e:
        print(f"Error: {e}. Please make sure both input files exist.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
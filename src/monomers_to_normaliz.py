# File name: multiset_to_vector.py
#
# This script reads a text file (monomers.txt) where each line represents a multiset of domains (called a monomer),
# converts each monomer into a vector representation, and outputs the results to two files:
# 1. A space-separated text file (vectors.txt) with vector representations and configurable singleton vectors.
# 2. An 'eqs.in' file with the transposed matrix and additional information.
# After processing, it outputs the number of domains and monomers to the command line.
#
# Input: A text file where each line contains space-separated domain names.
#        - Lines may start with text followed by a colon, which is ignored.
#        - Lines may end with a comma followed by text, which is ignored.
#        - Domains followed by a star (e.g., 'dom1*') represent negative counts.
# Output: 
#   - A space-separated text file (vectors.txt) where:
#     - The first line contains the domain names (header)
#     - Subsequent lines contain the vector representation of the input monomers
#     - The final set of lines contains singleton vectors (configurable)
#   - An 'eqs.in' file containing:
#     - The transposed matrix (excluding the header)
#     - Additional information as specified
#   - Command line output with the number of domains and monomers
#
# The script handles cases where some monomers may not contain all domains.
# The order of domains in the output is based on their first appearance in the input file.
# Domains with stars are treated as negative counts in the vector representation.

from collections import OrderedDict

# Flags to enable/disable singleton types
INCLUDE_NEGATIVE_SINGLETONS = True
INCLUDE_POSITIVE_SINGLETONS = True

def read_input_file(file_path):
    """
    Read the input file and return a list of monomers.
    Ignore text before a colon at the start of a line and after a comma at the end.
    """
    monomers = []
    with open(file_path, 'r') as file:
        for line in file:
            # Remove leading/trailing whitespace
            line = line.strip()
            
            # Ignore empty lines
            if not line:
                continue
            
            # Remove text before colon if present
            if ':' in line:
                line = line.split(':', 1)[1].strip()
            
            # Remove text after comma if present
            if ',' in line:
                line = line.rsplit(',', 1)[0].strip()
            
            monomers.append(line)
    return monomers

def get_unique_domains(monomers):
    """
    Extract unique domains from all monomers and maintain their order of appearance.
    Remove stars from domain names for consistency.
    """
    unique_domains = OrderedDict()
    for monomer in monomers:
        for domain in monomer.split():
            unique_domains[domain.rstrip('*')] = None
    return list(unique_domains.keys())

def create_vector_representation(monomer, domain_order):
    """
    Convert a monomer to its vector representation based on the domain order.
    Domains with stars are counted as negative.
    """
    domain_counts = {}
    for domain in monomer.split():
        base_domain = domain.rstrip('*')
        count = -1 if domain.endswith('*') else 1
        domain_counts[base_domain] = domain_counts.get(base_domain, 0) + count
    
    return [domain_counts.get(domain, 0) for domain in domain_order]

def create_singleton_vectors(domain_order):
    """
    Create singleton vectors for each domain based on the configuration flags.
    """
    singleton_vectors = []
    for i in range(len(domain_order)):
        if INCLUDE_NEGATIVE_SINGLETONS:
            negative_vector = [0] * len(domain_order)
            negative_vector[i] = -1
            singleton_vectors.append(negative_vector)
        if INCLUDE_POSITIVE_SINGLETONS:
            positive_vector = [0] * len(domain_order)
            positive_vector[i] = 1
            singleton_vectors.append(positive_vector)
    return singleton_vectors

def write_vectors_txt(file_path, domain_order, vectors):
    """
    Write the domain order and vectors to a space-separated text file.
    """
    with open(file_path, 'w') as file:
        file.write(' '.join(domain_order) + '\n')
        for vector in vectors:
            file.write(' '.join(map(str, vector)) + '\n')

def write_eqs_in(file_path, vectors):
    """
    Write the eqs.in file with the transposed matrix and additional information.
    """
    transposed = list(map(list, zip(*vectors)))
    with open(file_path, 'w') as file:
        file.write(f"amb_space {len(vectors)}\n")
        file.write(f"equations {len(transposed)}\n")
        for row in transposed:
            file.write(' '.join(map(str, row)) + '\n')
        file.write("HilbertBasis\n")

def process_file(input_file_path, vectors_txt_path, eqs_in_path):
    """
    Process the input file and write the results to the vectors.txt file and eqs.in file.
    Returns the number of domains and monomers.
    """
    # Read input file
    monomers = read_input_file(input_file_path)
    
    # Get unique domains in order of appearance
    domain_order = get_unique_domains(monomers)
    
    # Create vector representations and singleton vectors
    vectors = [create_vector_representation(monomer, domain_order) for monomer in monomers]
    singleton_vectors = create_singleton_vectors(domain_order)
    all_vectors = vectors + singleton_vectors
    
    # Write vectors.txt
    write_vectors_txt(vectors_txt_path, domain_order, all_vectors)
    
    # Write eqs.in
    write_eqs_in(eqs_in_path, all_vectors)

    return len(domain_order), len(monomers)

if __name__ == "__main__":
    input_file = "tmp_monomers.txt"  # Input file path
    vectors_txt = "vectors.txt"  # Output vectors file path
    eqs_in_file = "eqs.in"  # Path for the eqs.in file
    
    num_domains, num_monomers = process_file(input_file, vectors_txt, eqs_in_file)
    print(f"Processing complete. Results written to {vectors_txt} and {eqs_in_file}")
    print(f"domains: {num_domains}")
    print(f"monomers: {num_monomers}")
    
    # Print information about included singleton types
    singleton_types = []
    if INCLUDE_NEGATIVE_SINGLETONS:
        singleton_types.append("negative")
    if INCLUDE_POSITIVE_SINGLETONS:
        singleton_types.append("positive")
    print(f"Included singleton types: {' and '.join(singleton_types)}")
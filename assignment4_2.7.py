import requests
import os
import re
from Bio.Seq import Seq

working_dir = "."
fasta_path = os.path.join(working_dir, "mc1r.fasta")
txt_path = os.path.join(working_dir, "mc1r_homology_list.txt")

# Define the mygene URL
server_mygene = 'https://mygene.info/v3'
# Define the API endpoint to get the ensembl id of the MC1R gene in humans
endpoint_mygene = '/query?q=MC1R&species=human&fields=name,ensembl.gene'

# Send a GET request to the mygene using the specified server and endpoint
data = requests.get(server_mygene + endpoint_mygene)
data_json = data.json()

# Extract the ensembl gene ID for the MC1R gene from the JSON response
gene_id = data_json['hits'][0]['ensembl']['gene']

# Write DNA sequence of MC1R gene into a fasta file
with open(fasta_path, 'w') as fasta:
    json_header = {'Content-Type': 'text/x-fasta'}
    ensembl_server = 'http://rest.ensembl.org'
    endpoint_2 = '/sequence/id/{0}?type=genomic'.format(gene_id)
    dna_fasta = requests.get(ensembl_server + endpoint_2, headers=json_header)
    fasta.write(dna_fasta.text)

# Extract the sequence from fasta file
with open(fasta_path, 'r') as f:
    fasta_content = f.read()
    sequence = "".join(fasta_content.split("\n")[1:])

# Define a function for finding the longest ORF and translating it into an amino acid sequence
def find_longest_orf_and_translate(sequence):
    longest_orf = ''
    for i in range(3):
        for orf in re.findall(r'ATG(?:[ATGC]{3})*?(?:TAA|TAG|TGA)', sequence[i:]):
            if len(orf) > len(longest_orf):
                longest_orf = orf
    longest_orf_seq = Seq(longest_orf)
    amino_acid_seq = longest_orf_seq.translate()
    return longest_orf_seq, amino_acid_seq

longest_orf, amino_acid_seq = find_longest_orf_and_translate(sequence)

# Append amino acid sequence to fasta file
with open(fasta_path, 'a') as fasta:
    fasta.write("\n*** longest orf dna sequence ***\n\n")
    fasta.write(str(longest_orf))
    fasta.write("\n\n*** amino acid sequence received from the longest orf dna conversion ***\n\n")
    fasta.write(str(amino_acid_seq))

# Get information about species that have genes homologous to MC1R
endpoint_homologous = '/homology/id/human/{0}?content-type=application/json;type=orthologues;format=condensed'.format(gene_id)
data_homologous = requests.get(ensembl_server + endpoint_homologous)
data_homologous_json = data_homologous.json()

# Initialize an empty set to store unique species names
unique_species = set()

# Extract unique species names
for orthologues in data_homologous_json['data'][0]['homologies']:
    species = orthologues['species']
    unique_species.add(species)

# Write unique species names into the txt file
with open(txt_path, 'w') as txt:
    txt.write("Other species that have genes homologous to MC1R are:\n")
    for each in unique_species:
        txt.write(each + '\n')
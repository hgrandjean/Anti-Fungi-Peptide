from kmer_parser import find_kmer, parse_fasta_file
from generate_peptide import score_kmers
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO

peptides = pd.read_excel('results/de_novo_peptide_library.xlsx') #for excel files
file_name = "resources/filtered_positive_db.fasta" #for fasta files

reduce = 6

'''
Scoring of each descriptor found in the given peptide using previously computed scores. 
Uses find_kmer function from kmer_parser
'''

# Loading score from computed .tsv file
score_file = 'results/descriptors_activity_scores.tsv'

print (f"Loading descriptors scores from file: {score_file}")
score_dict: dict = {}
with open(score_file, "r" ) as scores:
    for line in scores:
        key, value = line.removesuffix('\n').split('\t')
        value =  value.strip('][').split(', ')
        value = [float(x) for x in value]
        score_dict[key] = value
print ("Finished loading scores")
print('####################################################################################################################################################### \n \n \n')

scores = []

# For .fasta files:
multi_fasta = parse_fasta_file ("resources/filtered_positive_db.fasta")
for fasta in multi_fasta:
    seq_db = fasta.seq
    scores.append(score_kmers(seq_db, r_dict = reduce, score_dictionary = score_dict))

sns.distplot(scores, bins = 20).set(title = "Distribution of peptide scores from positive DB", xlabel = "Scores", ylabel = "Count")
plt.ylim(0, 1)
plt.xlim(-5, 25)

plt.show()

scores = []

# For Excel files
for seq in peptides["Peptide sequence"]:
    scores.append(score_kmers(seq, r_dict=reduce ,score_dictionary = score_dict))

# Peptides["Activity score"] = scores

sns.distplot(scores, bins = 20).set(title="Distribution of generated peptides scores", xlabel="Scores", ylabel="Count")
plt.ylim(0, 1)
plt.xlim(-5, 25)

plt.show()
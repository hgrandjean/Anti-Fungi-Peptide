from kmer_parser import find_kmer
import numpy as np
import random
import sys
from scipy.signal import find_peaks
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import ProtParamData
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
def parse_fasta_file (file_name):
    multi_fasta = [record for record in SeqIO.parse(file_name, "fasta")]
    print(f"Ended parsing of {file_name}")
    return multi_fasta

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

def score_kmers(pep_seq: str, r_dict: int, score_dictionary = None) -> float:
    """
    pep_seq: peptide sequence
    r_dict: variable precising the reduction of the amino acid dictionary (20 AA) to a specific one (see RED dictionaries in kmer_parser.py)
    score_dictionary: dictionary with scores and kmer sequences
    """

    if score_dictionary is None:
        score_dictionary = score_dict

    kmer_score = 0
    size = min(len(pep_seq), 5)

    if size <= 2:
        gap = 0
    else:
        gap = size - 2

    kmers = find_kmer (sequence = pep_seq, kmer_size = size, n_gap= gap, r_dict = reduce)

    for kmer in kmers:
        if kmer in score_dictionary.keys() :
            kmer_score += score_dictionary[kmer][2]

    return kmer_score/len(pep_seq)

scores= []
'''
#for fasta files
multi_fasta = parse_fasta_file (file_name)
for fasta in multi_fasta:
    seq = fasta.seq
'''

for seq in peptides["Peptide sequence"]: #for excel files
    scores.append(score_kmers(seq, r_dict=reduce ,score_dictionary = score_dict))

#peptides["Activity score"]=scores

sns.distplot(scores, bins=20).set(title="Distribution of generated peptides scores", xlabel="Scores", ylabel="Count")

'''
plt.hist(scores, bins=1)
plt.plot(scores)
'''
plt.show()
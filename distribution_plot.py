from kmer_parser import parse_fasta_file
from generate_peptide import score_kmers, load_descriptors_score
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Selected reduction dictionary
REDUCE = 6

PEPTIDES = pd.read_excel('results/de_novo_peptide_library.xlsx') # For Excel files
POSITIVE_DB = "resources/filtered_positive_db.fasta"

# Loading score from computed .tsv file
SCORE_FILE = 'results/descriptors_activity_scores.tsv'

'''
Scoring of each descriptor found in the given peptide using previously computed scores. 
Uses find_kmer function from kmer_parser
'''

print (f"Loading descriptors scores from file: {SCORE_FILE}")
score_dictionary = load_descriptors_score(SCORE_FILE)
print ("Finished loading scores")
print('####################################################################################################################################################### \n \n \n')

scores_DB = []

# For .fasta files:
multi_fasta = parse_fasta_file (POSITIVE_DB)
for fasta in multi_fasta:
    seq_db = fasta.seq
    scores_DB.append(score_kmers(seq_db, r_dict = REDUCE, score_dict = score_dictionary))

'''
2 separate plots for each database

sns.histplot(scores_DB, bins = 20, kde = True).set(title ="Distribution of peptide scores from positive DB", xlabel ="Scores", ylabel ="Count")
sns.histplot(PEPTIDES["Activity score"], bins = 20, kde = True).set(title="Distribution of generated peptides scores", xlabel="Scores", ylabel="Count")
'''

'''
The plots for the peptides of a positive database and generated set in a single frame
'''

fig, ax = plt.subplots()
for a in [scores_DB, PEPTIDES["Activity score"]]:
    sns.distplot(a, bins = 20, ax = ax, kde = True).set(title="Distribution of scores", xlabel="Scores", ylabel="Density")

plt.legend(labels = ['Positive database', 'Generated peptides'])
plt.show()

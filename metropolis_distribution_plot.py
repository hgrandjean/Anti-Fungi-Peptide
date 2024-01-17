from kmer_parser import parse_fasta_file
from metropolis_generate_peptide import score_kmers, load_descriptors_score
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Selected reduction dictionary
REDUCE = 6

"""
Load best peptides from the positive database used for the de novo generation of peptides, and the best
de novo generated peptides
"""

BEST_POSITIVE_DB = pd.read_excel("results/de_novo_old_peptide.xlsx")
BEST_PEPTIDES = pd.read_excel("results/de_novo_best_peptide.xlsx")

# Loading score from computed .tsv file
SCORE_FILE = "results/descriptors_activity_scores.tsv"

"""
The plots for the peptides of a positive database and generated set in a single frame
"""

fig, ax = plt.subplots()
for a in [BEST_POSITIVE_DB["Activity score"], BEST_PEPTIDES["Activity score"]]:
    sns.distplot(a, bins=20, ax=ax, kde=True).set(title="Distribution of scores", xlabel="Scores", ylabel="Density")

plt.legend(labels=["Positive database", "Generated peptides"])
plt.show()

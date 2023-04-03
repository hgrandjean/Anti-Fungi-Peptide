### import dependencies ####
from kmer_parser import reduce_seq , gap_kmer , find_kmer
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns 
import math 
import random

### scoring each descriptor found in the given peptide using score previously computed ###
### uses find_kmer function from kmer_parser ###
def score_kmers(pep_seq , reduce, score_dict ):
      '''
      pep_seq : str of peptide sequence
      reduce : bool precise reduction of the amino acid dictionnary (20 AA) to 
                a specific one (see REDdictionnaries) 
      score_dict : dict loadded scores and kmer sequences 
      '''
      kmer_score = 0
      size = min(len(pep_seq), 5)
      if size <= 2 : gap = 0 
      else : gap = size - 2 
      kmers = find_kmer (sequence = pep_seq , kmer_size = size , ngap = gap , reduce = reduce )
      for kmer in kmers:
          if kmer in score_dict.keys() : 
              kmer_score += score_dict[kmer][2]
              # print(kmer," ", score_dict[kmer][2])
      return kmer_score/len(pep_seq)
    
    
score_file = "unique_set.tsv"


### loading score from computed tsv file ###
print (f"Loading descriptors scores from file : {score_file}")
score_dict = {}
with open(score_file, "r" ) as scores :
    for line in scores:
        key, value = line.removesuffix('\n').split('\t')
        value =  value.strip('][').split(', ')
        value = [float(x) for x in value ]
        score_dict[key] = value 
print ("Finished loading scores \nStarting scoring peptides")


### loading of Real IC50 database ### 
data= pd.read_excel("antimicrobial_peptide_with_IC_50.xlsx")  # SMAP-18 RGLRRLGRKIAHGVKKYG peptide 

### computation of scores for each peptide ###
score =[]
IC50 = []
for index, row in data.iterrows():
  peptide = row['Sequence']  
  IC50.append( row['Rel IC50'])
  score.append(score_kmers(peptide,6,score_dict))
data['score']=score


### visualisation of multidimensionnal data ###
fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')
n=100
ax.scatter(data["score"],data["Rel IC50"],data["hydrophobic_moment "],color="red")
ax.set_xlabel("score")
ax.set_ylabel("Rel IC50")
ax.set_zlabel("hydrophobic_moment")
plt.show()


### generation and optimisation of a peptide sequence
peptide = "AAAAAAAAAAAAAAAAAA"
score_kmers(peptide,6,score_dict)
bootstrap_iterations = 100
score_evolution = []
bootstrap = range(bootstrap_iterations)


AA_list = {"A": 1, "C": 2, "D": 3, "E": 4, "F": 5, "G": 6, "H": 7, "I": 8, "K": 9, "L": 10, "M": 11, "N": 12, "P": 13, "Q": 14, "R": 15, "S": 16, "T": 17, "V": 18, "W": 19, "Y": 20}

for i in bootstrap:
    # randomisation of mutation in the peptide sequence should be applied to biological form (To develop)
    random_index = random.randint(0, len(peptide) - 1)
    
    #replacing the amino acid selected to a knew one 
    random_amino_acid = peptide[random_index]
    new_amino_acid = random.choice(list(AA_list.keys()))
    new_peptide = peptide[:random_index] + new_amino_acid + peptide[random_index+1:]

    # Calculating scores of previous and new peptides sequences 
    peptide_score = score_kmers(peptide,6,score_dict)
    score_evolution.append(peptide_score)
    new_peptide_score = score_kmers(new_peptide,6,score_dict)

    # The peptide is selected if new score is higher 
    score_difference = new_peptide_score - peptide_score

    if score_difference > 0:
        peptide = new_peptide

    #addtion of randomness : if the mutation is not too much unfavored by the env then it can appen to be selected (To develop)
    
    #else:
    #   probability_of_acceptance = 1**(score_difference/100000)
    #    if random.random() < probability_of_acceptance:
    #        peptide = new_peptide

   
plt.scatter(bootstrap, score_evolution)   
plt.show() 
'''
'''

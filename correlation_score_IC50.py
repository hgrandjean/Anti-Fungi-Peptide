from generate_peptide import score_kmers
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from sklearn import linear_model
import statsmodels.api as sm

# Selected reduction dictionary
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

    kmers = find_kmer (pep_seq,size, gap , reduce) #find_kmer (sequence, kmer_size, ngap, reduce):

    for kmer in kmers:
        if kmer in score_dictionary.keys() :
            kmer_score += score_dictionary[kmer][2]

    return kmer_score/len(pep_seq)


AMPs_DB = pd.read_excel('resources/AMPs_DB_IC50.xlsx')


scores= []
for seq in AMPs_DB["Sequence"]:
  scores.append(score_kmers(seq, r_dict=reduce ,score_dictionary = score_dict))

AMPs_DB["score"]=scores
AMPs_DB['log_IC50'] = np.log10(AMPs_DB['Rel IC50'])
AMPs_DB = AMPs_DB[AMPs_DB['select']==1]


train_x = AMPs_DB[['score','hydrophobic fraction']] #'hydrophobic fraction'
train_y = AMPs_DB['log_IC50']

regr = linear_model.LinearRegression()
regr.fit(train_x,train_y)

# with sklearn
regr = linear_model.LinearRegression()
regr.fit(train_x, train_y)

print('Intercept: \n', regr.intercept_)
print('Coefficients: \n', regr.coef_)

# with statsmodels
train_x = sm.add_constant(train_x) # adding a constant
 
model = sm.OLS(train_y, train_x).fit()
predictions = model.predict(train_x) 
 
print_model = model.summary()
print(print_model)

fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')
n=100
ax.scatter(AMPs_DB["score"],AMPs_DB["hydrophobic moment"],AMPs_DB["log_IC50"],color="red")
ax.set_xlabel("score")
ax.set_ylabel("hydrophobic moment")
ax.set_zlabel("log_IC50")
plt.show() 

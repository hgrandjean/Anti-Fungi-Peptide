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
from sklearn import linear_model
import statsmodels.api as sm
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.svm import SVR
from sklearn.metrics import classification_report, confusion_matrix

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
for seq in AMPs_DB["sequence"]:
  scores.append(score_kmers(seq, r_dict=reduce ,score_dictionary = score_dict))

AMPs_DB["score"]=scores
AMPs_DB['log_IC50'] = np.log10(AMPs_DB['rel_IC50'])

pred_score_0=[]
pred_score_1=[]
for i in range(0,100):
  #split and prepare dataset
  train = AMPs_DB.sample(frac = 0.75)
  x_train_data= train[['score','hydrophobic_moment','a3v_Sequence_Average']]
  y_train_data=train['select']
  x_train_col_list =  x_train_data.values.tolist()
  y_train_col_list =  y_train_data.values.tolist()
  
  
  test = AMPs_DB.drop(train.index)
  x_test_data= test[['score','hydrophobic_moment','a3v_Sequence_Average']]
  y_test_data=test['select']
  x_test_col_list =  x_test_data.values.tolist()
  y_test_col_list =  y_test_data.values.tolist()
  
  
  
  classifier = SVC(kernel='linear', random_state = 0) #kernel = 'linear' ; 'poly'; 'rbf' or gamma = 'auto'
  classifier.fit(x_train_col_list, y_train_col_list)
  #Prediction sur le Test set
  y_pred = classifier.predict(x_test_col_list)
  
  
  cm = confusion_matrix(y_test_col_list,y_pred)
  #sns.heatmap(cm, annot=True, fmt='d').set_title('Confusion matrix of linear SVM')
  #print(classification_report(y_test_col_list, y_pred, output_dict=True)['0']['precision'])
  pred_score_0.append(classification_report(y_test_col_list, y_pred, output_dict=True)['0']['precision']*100)
  pred_score_1.append(classification_report(y_test_col_list, y_pred, output_dict=True)['1']['precision']*100)  
  '''
  
  regr = SVR()
  regr.fit(x_train_col_list, y_train_col_list)
  y_pred= regr.predict(x_test_col_list)
  
  plt.scatter(y_test_col_list, y_pred)
  '''
  
  #plt.show()
df={
  "pred_score_0" : pred_score_0,
  "pred_score_1" : pred_score_1
}

df=pd.DataFrame(df)
print(df)
sns.displot(df , x="pred_score_0", y="pred_score_1", binwidth=(5, 5), cbar=True)
plt.show()

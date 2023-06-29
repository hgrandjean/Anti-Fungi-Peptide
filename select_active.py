import joblib
import os
import shutil 
import sys



import pandas as pd 
from datetime import date 



today = date.today()
today = today.strftime("%Y%M%D")

def generate_fasta_file(sequences, names, scores,  output_file):
    with open(output_file, 'w') as fasta_file:
        for i in range(len(sequences)):
            fasta_file.write(f'>P{today}_{names[i]}_S_{scores[i]}\n')
            fasta_file.write(f'{sequences[i]}\n')

# define path to files 
model_path = 'resources/SVC_model.sav'
path_aggrescan  = sys.argv[1]
path_tango = sys.argv[2]
path_peptide_library = 'results/de_novo_peptide_library.xlsx'
output_file = 'results/selected_peptides.fasta'

# load the model from disk
model_SVC = joblib.load(model_path)

aggrescan_df = pd.read_csv(path_aggrescan, sep= '\t')
tango_df = pd.read_csv(path_tango , sep= '\t')
peptide_library = pd.read_excel(path_peptide_library)


# reformat database 
tango_df = tango_df.drop('Unnamed: 0', axis=1)
aggrescan_df = aggrescan_df.iloc[:, :-2]
aggrescan_df= aggrescan_df.set_index('Sequence Name')
aggrescan_df=aggrescan_df.T
colnames= [ 'a3vSA',	'nHS',	'NnHS',	'AAT',	'THSA',	'TA',	'AATr',	'THSAr',	'Na4vSS']
aggrescan_df.columns = colnames


# get all used variables 
peptide_library['agg']= tango_df['AGG']
peptide_library['a3vSA'] = aggrescan_df['a3vSA'].values.tolist()

peptide_library_to_list = peptide_library[["score", "hydrophobicity_average", "a3vSA","agg"]].values.tolist()

# run prediction 
y_pred = model_SVC.predict(peptide_library_to_list)
peptide_library['selection']=y_pred

# save prediction 
peptide_library.to_excel('results/de_novo_peptide_library.xlsx')
selected=peptide_library[peptide_library['selection']== 1 ]
names = range(len(selected['peptide_sequence'].values.tolist()))
generate_fasta_file(selected['peptide_sequence'].values.tolist(),names,selected['score'].values.tolist() ,output_file)

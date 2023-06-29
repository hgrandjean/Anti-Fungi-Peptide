from generate_peptide import score_kmers, load_descriptors_score
import numpy as np
import math
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import ProtParamData
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.svm import SVC
from sklearn.metrics import classification_report, confusion_matrix
import joblib

# Selected reduction dictionary
REDUCE = 6
SCORE_FILE = "results/descriptors_activity_scores.tsv"

"""
Scoring of each descriptor found in the given peptide using previously computed scores. 
Uses find_kmer function from kmer_parser
"""

# Loading score from computed .tsv file
print(f"Loading descriptors scores from file: {SCORE_FILE}")
score_dictionary = load_descriptors_score(SCORE_FILE)
print("Finished loading scores")

print(
    "####################################################################################################################################################### \n \n \n"
)


# Global descriptor scoring function
def pep_hydrophobicity_analysis(pep_seq: str) -> list:
    """
    Physical analysis of peptide sequence based on residues compositions
    arg:
    pep_seq: peptide sequence
    return:
    Hydrophobicity GRAVY (GRand AVerage of hYdrophaty) for complete peptide
    Alternative return: hydrophobicity profile of the peptide with the scale of 2
    """
    pa = ProteinAnalysis(pep_seq)

    """
    Kyte-Doolitle hydrophobicity profile grand average 
    """
    hydrophobicity = pa.gravy()
    hydrophobicity_scale = pa.protein_scale(ProtParamData.kd, 2, edge=1.0)
    return hydrophobicity


def calculate_moment(array, angle=100):
    """Calculates the hydrophobic dipole moment from an array of hydrophobicity
    values. Formula defined by Eisenberg, 1982 (Nature). Returns the average
    moment (normalized by sequence length)

    uH = sqrt(sum(χcos(i*d))**2 + sum(χsin(iδ))**2),
    where i is the amino acid index and δ is an angular value in
    º (100º for α-helix, 180º for β-sheet).

    Extracted from: https://github.com/JoaoRodrigues/hydrophobic_moment/blob/main/hydrophobic_moment.py

    arg
    array:  is a scale of hydrophobicity for peptide sequence
    return
    hydrophobic:  average moment based on the Eisenberg formula
    """

    sum_cos, sum_sin = 0.0, 0.0
    for i, hv in enumerate(array):
        rad_inc = ((i * angle) * math.pi) / 180.0
        sum_cos += hv * math.cos(rad_inc)
        sum_sin += hv * math.sin(rad_inc)
    # print(sum_cos, sum_sin, rad_inc)

    return math.sqrt(sum_cos**2 + sum_sin**2) / len(array)


# Import database of IC50 published in Fjell, C. D. et al. Identification of Novel Antibacterial Peptides by Chemoinformatics and Machine Learning. J. Med. Chem. 52, 2006–2015  (2009).
AMPs_DB = pd.read_excel("resources/AMPs_DB_IC50.xlsx")

# compute scores and hydrophobicity for all peptides and add them to the dataframe
scores = []
hydrophobicity_profile = []
for seq in AMPs_DB["sequence"]:
    scores.append(score_kmers(seq, r_dict=REDUCE, score_dict=score_dictionary))
    hydrophobicity_profile.append(pep_hydrophobicity_analysis(seq))
    # hydrophobicity_profile.append(calculate_moment(pep_physical_analysis(seq)))

AMPs_DB["score"] = scores
AMPs_DB["log_IC50"] = np.log10(AMPs_DB["rel_IC50"])
AMPs_DB["hydrophobicity_average"] = hydrophobicity_profile


# Perform multiple iteration of the SVC to see false and true positive discovery

pred_score_0 = []
pred_score_1 = []
for i in range(0, 1000):
    # Split and prepare dataset for train and test in a specified ratio
    # Train dataset

    train = AMPs_DB.sample(frac=0.75)
    x_train_data = train[["score", "hydrophobicity_average", "a3vSA",'agg']]
    y_train_data = train["select"]
    x_train_col_list = x_train_data.values.tolist()
    y_train_col_list = y_train_data.values.tolist()

    # Test dataset
    test = AMPs_DB.drop(train.index)
    x_test_data = test[["score", "hydrophobicity_average", "a3vSA","agg"]]
    y_test_data = test["select"]
    x_test_col_list = x_test_data.values.tolist()
    y_test_col_list = y_test_data.values.tolist()

    # Model creation and training
    classifier = SVC(kernel="linear", random_state=0)  # kernel = 'linear' ; 'poly'; 'rbf' or gamma = 'auto'
    classifier.fit(x_train_col_list, y_train_col_list)

    # Prediction on the test dataset
    y_pred = classifier.predict(x_test_col_list)

    # Model evaluation based on test dataset
    cm = confusion_matrix(y_test_col_list, y_pred)

    # sns.heatmap(cm, annot=True, fmt='d').set_title('Confusion matrix of linear SVM')
    # print(classification_report(y_test_col_list, y_pred, output_dict=True)['0']['precision'])
    if i == 0:
        actual_pred_score_0 = classification_report(y_test_col_list, y_pred, output_dict=True)["0"]["precision"] * 100
        actual_pred_score_1 = classification_report(y_test_col_list, y_pred, output_dict=True)["1"]["precision"] * 100
    else:
        new_pred_score_0 = classification_report(y_test_col_list, y_pred, output_dict=True)["0"]["precision"] * 100
        new_pred_score_1 = classification_report(y_test_col_list, y_pred, output_dict=True)["1"]["precision"] * 100

        if new_pred_score_0 >= actual_pred_score_0 and new_pred_score_1 >= actual_pred_score_0:
            actual_pred_score_0 = new_pred_score_0
            actual_pred_score_1 = new_pred_score_1
            decision = classifier.decision_function(x_test_col_list)
            parameters = classifier.get_params()
            training_set_x = x_train_data
            training_set_y = y_train_data
            testing_set_x = x_test_data
            testing_set_y = y_test_data

    pred_score_0.append(classification_report(y_test_col_list, y_pred, output_dict=True)["0"]["precision"] * 100)
    pred_score_1.append(classification_report(y_test_col_list, y_pred, output_dict=True)["1"]["precision"] * 100)


# save the model to disk
filename = 'resources/SVC_model.sav'
joblib.dump(classifier, filename)


with open("results/best_svc.txt", "w") as file:
    file.write(
        f"Best confidence scores:\n\tConfidence of inactive attribution: {actual_pred_score_0}%\n\n\tConfidence of active attribution: {actual_pred_score_1}%"
    )
    file.write(f"\n\nDecision function: \n{decision}\n\nParameters: \n{parameters}")
    file.write(
        f"\n\nTraining set:\n\tX: \n{training_set_x}\n\n\tY: \n{training_set_y}\n\nTesting set:\n\tX: \n{testing_set_x}\n\n\tY: \n{testing_set_y}"
    )


"""
Linear regression model based on SVM did not show linear relationship between analysed values

regr = SVR()
regr.fit(x_train_col_list, y_train_col_list)
y_pred= regr.predict(x_test_col_list)
  
plt.scatter(y_test_col_list, y_pred)

"""

# Plot result of all performed iterations and confidence value for true and false positive discovery
df = {"pred_score_0": pred_score_0, "pred_score_1": pred_score_1}

df = pd.DataFrame(df)
print(df)
sns.histplot(df, x="pred_score_0", y="pred_score_1", binwidth=(5, 5), cbar=True).set(
    title="Distribution of false and true discovery rate\nin active and inactive peptides",
    xlabel="Confidence of inactive attribution",
    ylabel="Confidence of active attribution",
)
plt.xlim(None, 100)
plt.ylim(None, 100)
plt.show()

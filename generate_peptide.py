from kmer_parser import find_kmer
import numpy as np
import random
import sys
from scipy.signal import find_peaks
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import ProtParamData
import pandas as pd

'''
Loading interface
'''

print( '                                                                                                                                                       ' )
print( '                                                                                                                                                       ' )
print( '                                                                                                                    ==========----                     ' )
print( '                                                                                                                    ++++++++++******++=                ' )
print( '                                                                                                                    -+++++++++********##-              ' )
print( '  ---------  ---      --- ----------- ---  ---------- ---      --- ---     ---  ---------   ----------  ---          ++++++++++*********+              ' )
print( ' #@@@@@@@@@#-#@@#=   +@@*=@@@@@@@@@@% @@@ *@@@@@@@@@@=*@@+    +@@#=@@@*    @@@ #@@@@@@@@@#=+%@@@@@@@@%+=@@#          =+++++++++*******+++=             ' )
print( ' @@@=====%@@-#@@@#=  +@@*-===+@@%==== @@@ %@@+=======-*@@+    +@@#=@@@%*   @@@ @@@+====#%%=#@@*====*@@*=@@#           +=======++***+++++++             ' )
print( ' @@@     %@@-#@@@@#- +@@*    -@@#     @@@ %@@-        *@@+    +@@#=@@@@%*  @@@ @@@-        #@@+    +@@*=@@#                  -=+++++++++++=            ' )
print( ' @@@@@@@@@@@-#@@*#@#=+@@*    -@@#     @@@ %@@@@@@@=   *@@+    +@@#=@@%+%@*-@@@ @@@-  #@@@#-#@@@@@@@@@@*=@@#               =+***+++++++++++             ' )
print( ' @@@+++++%@@-#@@= *@@@@@*    -@@#     @@@ %@@+++++-   *@@+    +@@#=@@% -@@@@@@ @@@-  =+%@@=#@@#++++#@@*=@@#            =+******++++++++=-              ' )
print( ' @@@     %@@-#@@=  *@@@@*    -@@#     @@@ %@@-        *@@+    +@@#=@@%   @@@@@ @@@-    #@@=#@@+    +@@*=@@#         -+**********+++++-                 ' )
print( ' @@@     %@@-#@@=   *@@@*    -@@#     @@@ %@@-        *@@@@@@@@@@#=@@%    @@@@ @@@@@@@@@@@=#@@+    +@@*=@@@@@@@@@@@-##***********+-                    ' )
print( ' +++     +++-=++-    =++=    -+++     +++ +++          =++++++++= -+++     +++  +++++++++- =++-    -++= -++++++++++ *#********+-                       ' )
print( '                                                                                                                    =##******++***++++++===-           ' )
print( '                                                                                                                     *##**+++++++**********##+         ' )
print( '                                                                                                                     -#*++++++++++***********+-        ' )
print( '                                        -*******--*******-*******-+*******-*+-******+ =******+-******+                ==++++++++++********+++++        ' )
print( '                                        -@%+++#@*#@*+++++-@%+++#@*=++@@*++=@%+@#+++%@=%@+++++++@#+++#%=                   --===+++*****++++++++-       ' )
print( '                                        -@%+++#@*#@*+++  -@%+++#@*   %@-  =@%+@*   %@=%@*++=  +@#++++=                            =**++++++++++=       ' )
print( '                                        -@@*****-#@#**+  -@%*****-   %@-  =@%+@*   %@=%@***=   *****%@=                        =+***+++++++++++        ' )
print( '                                        -@%      #@*=====-@%         %@-  =@%+@#===%@=%@+======#*===%@=                     =+*******+++++++=-         ' )
print( '                                        -#*      -#######-#*         *#   -#*=######*-+######*-*#####*                    +***********++++-            ' )
print( '                                                                                                                         *#***********=-               ' )
print( '                                                                                                                         +##********+=----             ' )
print( '                                                                                                                         -##******++++*********++=     ' )
print( '                                                                                                                          +###*++++++++**********#*-   ' )
print( '                                                                                                                           #*++++++++++***********++   ' )
print( '                                                                                                                           -====++++++++********++++=  ' )
print( '                                                                                                                                     --==******++++++  ' )
print( '                                                                                                                                         *#****++++++  ' )
print( '                                                                                                                                          #*****+++-   ' )
print( '                                                                                                                                          =#****+-     ' )
print( '                                                                                                                                           **+-        ' )
print( '                                                                                                                                                       ' )
print( '                                                                                                                                                       ' )
print('\n \n \n #######################################################################################################################################################')

# Selected reduction dictionary
reduce = 6

def progress_bar(count, total, size = 100, sides = "[]", full = '#', empty = '.', prefix = ""):
    x = int(size * count/total)
    sys.stdout.write("\r" + prefix + sides[0] + full * x + empty * (size - x) + sides[1] + ' ' + str(count).rjust(len(str(total)),' ') + "/" + str(total))
    if count == total:
        sys.stdout.write("\n")

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

    kmers = find_kmer (sequence = pep_seq, kmer_size = size, n_gap= gap, r_dict = reduce)

    for kmer in kmers:
        if kmer in score_dictionary.keys() :
            kmer_score += score_dictionary[kmer][2]

    return kmer_score/len(pep_seq)

### Computation of peptide's physical properties
def pep_physical_analysis(pep_seq: str) -> list [str, float, float, float]:
    """
    pep_seq: peptide sequence
    """
    pa = ProteinAnalysis(pep_seq)
    helix_propensity = pa.secondary_structure_fraction()[0]
    charge = pa.charge_at_pH(7)
    
    """
    Kyte-Doolitle hydrophobicity profile
    """

    hydrophobicity_profile = pa.protein_scale(ProtParamData.kd, 2, edge = 1.0)

    # Identify the minimal distance between 2 peaks
    peaks, _ = find_peaks (hydrophobicity_profile, distance = 2)

    # Mean distance between hydrophilic maximums
    periodicity = np.mean (np.diff(peaks))
    return [pep_seq, helix_propensity * 100, charge, periodicity]


'''
### loading of Real IC50 database ### 
data= pd.read_excel("antimicrobial_peptide_with_IC_50.xlsx")  # SMAP-18 RGLRRLGRKIAHGVKKYG peptide 

### computation of scores for each peptide ###
score =[]
IC50 = []
helix_prob=[]
computed_charge =[]
space =[]
for index, row in data.iterrows():
  peptide = row['Sequence']  
  IC50.append( row['Rel IC50'])
  score.append(score_kmers(peptide,6,score_dict))
  physical_analysis = pep_physical_analysis(peptide)
  helix_prob.append(physical_analysis[1])
  computed_charge.append(physical_analysis[2])
  space.append(physical_analysis[3])

data["helix_prob"]=helix_prob
data["computed_charge"]=computed_charge
data["space"]=space
data["score"]=score
print(data)
data.to_excel("computed_db.xlsx")
'''

'''
Import substitution probabilities from PAM2 Dataframe relative to mutation frequencies with conservation excluded
'''

pam2_prob_matrix = pd.read_excel('resources/PAM_2_substitution_probabilities_formated.xlsx')
print (f"Loading PAM substitution probabilities")
AA_order = 'ARNDCQEGHILKMFPSTWYV'
print (f"Finished loading PAM substitution probabilities \n \n \n")

# Convert the Dataframe into a dictionary
pam2_prob_dict: dict = {}
for row in pam2_prob_matrix.to_numpy():
    pam2_prob_dict[row[0]]=row[1:]


''' Clustering -> not successful
### visualisation of multidimensional data ###
fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')
n=100
ax.scatter(data["score"],data["Rel IC50"],data["helix_prob"],color="red")
ax.set_xlabel("score")
ax.set_ylabel("Rel IC50")
ax.set_zlabel("helix_prob")
plt.show()
'''

'''
Generation and optimisation of a peptide sequence
'''

print('####################################################################################################################################################### \n \n \n')
n_pep = 30
pep_seq = "RGLRRLGRKIAHGVKKYG"
score_kmers (pep_seq, reduce, score_dict)
bootstrap_iterations = 1000
score_evolution: list[float] = []
helix_propensity_evolution: list[float] = []
charge_evolution: list[float] = []
periodicity_evolution: list[float] = []
peptides_generated: list[str] = []
bootstrap = range(bootstrap_iterations)
print("Starting generation of {npep}peptides from {bootstrap_iterations} bootstrap iterations \n ")
 
#Test
for p in range(0,n_pep):
    pep_seq = "RGLRRLGRKIAHGVKKYG"
    for i in bootstrap:
        # Randomisation of mutation location in the peptide sequence should be applied to biological form (To develop)
        random_index = random.randint(0, len(pep_seq) - 1)
        
        # Replacing the amino acid selected to a knew one
        random_amino_acid = pep_seq[random_index]
        prob = pam2_prob_dict[random_amino_acid]
        new_amino_acid = random.choices(AA_order, prob, k = 1)[0]
        new_peptide = pep_seq[:random_index] + new_amino_acid + pep_seq[random_index + 1:]
    
        # Calculating scores of previous and new peptide sequences
        peptide_score = score_kmers(pep_seq, reduce, score_dict)
        physical_analysis: list [str, float, float, float] = pep_physical_analysis(pep_seq)
        new_peptide_score = score_kmers(new_peptide, reduce, score_dict)
    
        # The peptide is selected if the new score is higher
        score_difference = new_peptide_score - peptide_score
    
        if score_difference > 0:
            pep_seq = new_peptide
    
        # Addition of randomness: if the mutation is not unfavored by the environment then it can happen to be selected (to develop)
        
        # else:
        #   probability_of_acceptance = 1**(score_difference/100000)
        #    if random.random() < probability_of_acceptance:
        #        peptide = new_peptide
        progress_bar (count = i + 1, total = bootstrap_iterations, size = 100, sides = "||", full = '#', empty = ' ', prefix = "Performing bootstraps... ")
    

    print("\nFinal peptide: \n")
    print(pep_seq)
    print(f"final score: {score_kmers(pep_seq, reduce, score_dict)}")
    print(f"final helix propensity: {round(physical_analysis[1], 2)} %")
    print(f"final global charge Q: {physical_analysis[2]}")
    print(f"final hydrophobicity frequency: {physical_analysis[3]}")
    print(f"{p+1} sequence generated over {n_pep} \n \n " )
    score_evolution.append(new_peptide_score)
    
    helix_propensity_evolution.append(physical_analysis[1])
    charge_evolution.append(physical_analysis[2])
    periodicity_evolution.append(physical_analysis[3])
    peptides_generated.append(pep_seq)

print(score_evolution)
final_df=pd.DataFrame(list(zip(peptides_generated, score_evolution, charge_evolution, periodicity_evolution)), columns = ["Peptide sequence", "Activity score" , "Net charge" , "Hydrophobicity-based periodicity"])
print(final_df)
final_df.to_excel('results/de_novo_peptide_library.xlsx')

#plt.plot(bootstrap, score_evolution, label = "Score")
#plt.plot(bootstrap, helix_propensity_evolution, label = "Helix propensity")
#plt.plot(bootstrap, charge_evolution, label = "Net charge")
#plt.plot(bootstrap, periodicity_evolution, label = "Hydrophobicity-based periodicity")
#plt.legend()
#plt.show() 
'''
'''

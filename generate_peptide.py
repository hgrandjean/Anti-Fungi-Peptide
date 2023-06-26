### import dependencies ####
import random
import sys

import numpy as np
import pandas as pd
from Bio.SeqUtils import (
    ProtParamData,  # Local https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/ProtParamData.py
)
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from scipy.signal import find_peaks

from kmer_parser import find_kmer

from rich import box, print as rich_print
from rich.panel import Panel
from rich.progress import BarColumn, Progress, SpinnerColumn, MofNCompleteColumn, TaskID, TextColumn, TimeElapsedColumn
from rich.table import Table
from rich.align import Align

print(
    "                                                                                                                                                       "
)
print(
    "                                                                                                                                                       "
)
print(
    "                                                                                                                    ==========----                     "
)
print(
    "                                                                                                                    ++++++++++******++=                "
)
print(
    "                                                                                                                    -+++++++++********##-              "
)
print(
    "  ---------  ---      --- ----------- ---  ---------- ---      --- ---     ---  ---------   ----------  ---          ++++++++++*********+              "
)
print(
    " #@@@@@@@@@#-#@@#=   +@@*=@@@@@@@@@@% @@@ *@@@@@@@@@@=*@@+    +@@#=@@@*    @@@ #@@@@@@@@@#=+%@@@@@@@@%+=@@#          =+++++++++*******+++=             "
)
print(
    " @@@=====%@@-#@@@#=  +@@*-===+@@%==== @@@ %@@+=======-*@@+    +@@#=@@@%*   @@@ @@@+====#%%=#@@*====*@@*=@@#           +=======++***+++++++             "
)
print(
    " @@@     %@@-#@@@@#- +@@*    -@@#     @@@ %@@-        *@@+    +@@#=@@@@%*  @@@ @@@-        #@@+    +@@*=@@#                  -=+++++++++++=            "
)
print(
    " @@@@@@@@@@@-#@@*#@#=+@@*    -@@#     @@@ %@@@@@@@=   *@@+    +@@#=@@%+%@*-@@@ @@@-  #@@@#-#@@@@@@@@@@*=@@#               =+***+++++++++++             "
)
print(
    " @@@+++++%@@-#@@= *@@@@@*    -@@#     @@@ %@@+++++-   *@@+    +@@#=@@% -@@@@@@ @@@-  =+%@@=#@@#++++#@@*=@@#            =+******++++++++=-              "
)
print(
    " @@@     %@@-#@@=  *@@@@*    -@@#     @@@ %@@-        *@@+    +@@#=@@%   @@@@@ @@@-    #@@=#@@+    +@@*=@@#         -+**********+++++-                 "
)
print(
    " @@@     %@@-#@@=   *@@@*    -@@#     @@@ %@@-        *@@@@@@@@@@#=@@%    @@@@ @@@@@@@@@@@=#@@+    +@@*=@@@@@@@@@@@-##***********+-                    "
)
print(
    " +++     +++-=++-    =++=    -+++     +++ +++          =++++++++= -+++     +++  +++++++++- =++-    -++= -++++++++++ *#********+-                       "
)
print(
    "                                                                                                                    =##******++***++++++===-           "
)
print(
    "                                                                                                                     *##**+++++++**********##+         "
)
print(
    "                                                                                                                     -#*++++++++++***********+-        "
)
print(
    "                                        -*******--*******-*******-+*******-*+-******+ =******+-******+                ==++++++++++********+++++        "
)
print(
    "                                        -@%+++#@*#@*+++++-@%+++#@*=++@@*++=@%+@#+++%@=%@+++++++@#+++#%=                   --===+++*****++++++++-       "
)
print(
    "                                        -@%+++#@*#@*+++  -@%+++#@*   %@-  =@%+@*   %@=%@*++=  +@#++++=                            =**++++++++++=       "
)
print(
    "                                        -@@*****-#@#**+  -@%*****-   %@-  =@%+@*   %@=%@***=   *****%@=                        =+***+++++++++++        "
)
print(
    "                                        -@%      #@*=====-@%         %@-  =@%+@#===%@=%@+======#*===%@=                     =+*******+++++++=-         "
)
print(
    "                                        -#*      -#######-#*         *#   -#*=######*-+######*-*#####*                    +***********++++-            "
)
print(
    "                                                                                                                         *#***********=-               "
)
print(
    "                                                                                                                         +##********+=----             "
)
print(
    "                                                                                                                         -##******++++*********++=     "
)
print(
    "                                                                                                                          +###*++++++++**********#*-   "
)
print(
    "                                                                                                                           #*++++++++++***********++   "
)
print(
    "                                                                                                                           -====++++++++********++++=  "
)
print(
    "                                                                                                                                     --==******++++++  "
)
print(
    "                                                                                                                                         *#****++++++  "
)
print(
    "                                                                                                                                          #*****+++-   "
)
print(
    "                                                                                                                                          =#****+-     "
)
print(
    "                                                                                                                                           **+-        "
)
print(
    "                                                                                                                                                       "
)
print(
    "                                                                                                                                                       "
)
print(
    "\n \n \n #######################################################################################################################################################"
)

SCORE_FILE = "results/descriptors_activity_scores.tsv"
PAM2_EXCEL_FILE = "resources/PAM_2_substitution_probabilities_formated.xlsx"
AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
DEFAULT_PEPTIDE = "RGLRRLGRKIAHGVKKYG"
NB_PEPTIDE = 50
NB_ITERATIONS = 1000


def default_progress():
    return Progress(
        SpinnerColumn(),
        TimeElapsedColumn(),
        BarColumn(),
        MofNCompleteColumn(),
        TextColumn("[bold blue]({task.description})"),
        transient=False,
    )


### scoring each descriptor found in the given peptide using score previously computed ###
### uses find_kmer function from kmer_parser ###
def score_kmers(pep_seq, reduce, score_dict):
    """
    pep_seq : str of peptide sequence
    reduce : bool precise reduction of the amino acid dictionnary (20 AA) to
              a specific one (see REDdictionnaries)
    score_dict : dict loadded scores and kmer sequences
    """
    kmer_score = 0
    size = min(len(pep_seq), 5)
    if size <= 2:
        gap = 0
    else:
        gap = size - 2
    kmers = find_kmer(sequence=pep_seq, kmer_size=size, ngap=gap, reduce=reduce)
    for kmer in kmers:
        if kmer in score_dict.keys():
            kmer_score += score_dict[kmer][2]
            # print(kmer," ", score_dict[kmer][2])
    return kmer_score / len(pep_seq)


def load_descriptors_score(score_file: str = SCORE_FILE) -> dict:
    score_dict = {}
    with open(score_file, "r") as scores:
        for line in scores:
            key, value = line.removesuffix("\n").split("\t")
            value = value.strip("][").split(", ")
            value = [float(x) for x in value]
            score_dict[key] = value
    return score_dict


### computation of peptide physical properties ###
def pep_physical_analysis(peptide):
    pa = ProteinAnalysis(str(peptide))
    helix_prob = pa.secondary_structure_fraction()[0]
    charge = pa.charge_at_pH(7)

    # get location of hydrophoilic residues
    hydropho = pa.protein_scale(
        ProtParamData.kd, 2, edge=1.0
    )  # hydrophobicity using a step of one kd = Kyle Doolittle score
    size = len(hydropho)
    peaks, _ = find_peaks(hydropho, distance=2)  # identify maximums miniaml distance between two peaks = 2
    space = np.mean(np.diff(peaks))  # mean space between hydrophilic maximums
    return [peptide, helix_prob * 100, charge, space]


def generate_prob_dict_from_excel(file_name: str = PAM2_EXCEL_FILE) -> dict:
    pam2_prob_matrix = pd.read_excel(file_name)

    # Convert the DataFrame to a dictionary
    pam2_prob_dict = {}
    for row in pam2_prob_matrix.to_numpy():
        pam2_prob_dict[row[0]] = row[1:]
    return pam2_prob_dict


def generate_peptides(
    aa_probs: dict, nb_peptide: int = NB_PEPTIDE, default_peptide: str = DEFAULT_PEPTIDE, bootstrap: int = NB_ITERATIONS
):
    Align(rich_print(Panel(f"Generation of {nb_peptide} peptides", style="light_slate_blue", expand=False)))
    # Test
    score_evolution = []
    helix_proba_evol = []
    charge_evol = []
    space_evolution = []
    peptides_generated = []

    for p in range(0, nb_peptide):
        peptide = default_peptide
        with default_progress() as progress:
            task_id: TaskID = progress.add_task(f"Generating peptide number {p+1} out of {nb_peptide}", total=bootstrap)
            for i in range(bootstrap):
                # randomisation of mutation location in the peptide sequence should be applied to biological form (To develop)
                random_index = random.randint(0, len(peptide) - 1)

                # replacing the amino acid selected to a knew one
                random_amino_acid = peptide[random_index]
                prob = aa_probs[random_amino_acid]
                new_amino_acid = random.choices(AA_ORDER, prob, k=1)[0]
                new_peptide = peptide[:random_index] + new_amino_acid + peptide[random_index + 1 :]

                # Calculating scores of previous and new peptides sequences
                peptide_score = score_kmers(peptide, 6, score_dict)
                # score_evolution.append(peptide_score)

                physical_analysis = pep_physical_analysis(peptide)
                # helix_proba_evol.append(physical_analysis[1])
                # charge_evol.append(physical_analysis[2])
                # space_evolution.append(physical_analysis[3])

                new_peptide_score = score_kmers(new_peptide, 6, score_dict)

                # The peptide is selected if new score is higher
                score_difference = new_peptide_score - peptide_score
                #print(f"current score: {peptide_score} ; new score: {new_peptide_score}")

                if new_peptide_score - peptide_score > 0:
                    peptide = new_peptide

                # addtion of randomness : if the mutation is not too much unfavored by the env then it can appen to be selected (To develop)

                # else:
                #   probability_of_acceptance = 1**(score_difference/100000)
                #    if random.random() < probability_of_acceptance:
                #        peptide = new_peptide

                # progress bar HERE
                progress.update(task_id, advance=1)

        '''
        rich_print("\nFinal peptide : \n")
        rich_print(peptide)
        rich_print(f"final score : {score_kmers(peptide,6,score_dict)}")
        rich_print(f"final helix probability : {round(physical_analysis[1], 2)} %")
        rich_print(f"final global charge Q : {physical_analysis[2]}")
        rich_print(f"final hydrophobicity frequency : {physical_analysis[3]}")
        '''

        table = Table(
            title=f"[i]Result of peptide[/] [b royal_blue1]{peptide}",
            show_header=False,
            box=box.ROUNDED,
            style="cyan",
            title_style=""
        )
        table.add_row("Final score", f"{score_kmers(peptide,6,score_dict)}")
        table.add_row("Final helix probability", f"{round(physical_analysis[1], 2)}%")
        table.add_row("Final global charge Q", f"{physical_analysis[2]}")
        table.add_row("Final hydrophobicity frequency", f"{physical_analysis[3]}")

        rich_print("\n", table, "\n")


        score_evolution.append(new_peptide_score)
        helix_proba_evol.append(physical_analysis[1])
        charge_evol.append(physical_analysis[2])
        space_evolution.append(physical_analysis[3])
        peptides_generated.append(peptide)
    #print(score_evolution)
    final_df = pd.DataFrame(
        list(zip(peptides_generated, score_evolution, charge_evol, space_evolution)),
        columns=["peptide_sequence", "score", "charge", "hydro_frequency"],
    )
    rich_print(Align(Panel("Generated peptides and their respectives properties", style="light_slate_blue", expand=False), align="center"))
    print(final_df)
    final_df.to_excel("results/de_novo_peptide_library.xlsx")


if __name__ == "__main__":
    # Getting scores from CSV file
    score_dict = load_descriptors_score()

    # Import substitution probabilities from PAM2 data frame relative to mutation frequencies with conservation excluded
    pam2_probs = generate_prob_dict_from_excel()

    ### generation and optimisation of a peptide sequence
    generate_peptides(pam2_probs)


"""
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
"""

""" Clustering -> not succesful
### visualisation of multidimensionnal data ###
fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')
n=100
ax.scatter(data["score"],data["Rel IC50"],data["helix_prob"],color="red")
ax.set_xlabel("score")
ax.set_ylabel("Rel IC50")
ax.set_zlabel("helix_prob")
plt.show()
"""

# plt.plot(bootstrap, score_evolution ,label = "score")
# plt.plot(bootstrap, helix_proba_evol ,label = "helix probability")
# plt.plot(bootstrap, charge_evol ,label = "charge")
# plt.plot(bootstrap, space_evolution ,label = "hydrophobicity space")
# plt.legend()
# plt.show()
"""
"""

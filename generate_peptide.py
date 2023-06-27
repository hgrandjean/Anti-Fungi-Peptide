import random

import numpy as np
import pandas as pd
from Bio.SeqUtils import (
    ProtParamData,  # Local https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/ProtParamData.py
)
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

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

# Selected reduction dictionary
REDUCE = 6

# Loading interface with the progress bar

def default_progress():
    return Progress(
        SpinnerColumn(),
        TimeElapsedColumn(),
        BarColumn(),
        MofNCompleteColumn(),
        TextColumn("[bold blue]({task.description})"),
        transient = False,
    )

def load_descriptors_score(score_file: str = SCORE_FILE) -> dict:
    with open(score_file, "r") as scores:
        score_dictionary: dict = {}
        for line in scores:
            key, value = line.removesuffix("\n").split("\t")
            value = value.strip("][").split(", ")
            value = [float(x) for x in value]
            score_dictionary[key] = value
    return score_dictionary

"""
Scoring of each descriptor found in the given peptide using previously computed scores. 
Uses find_kmer function from kmer_parser
"""

def score_kmers(pep_seq, r_dict: int, score_dict) -> float:

    """
    pep_seq: peptide sequence
    r_dict: variable precising the reduction of the amino acid dictionary (20 AA) to a specific one (see RED dictionaries in kmer_parser.py)
    score_dictionary: dictionary with scores and kmer sequences
    """

    if score_dict is None:
        score_dict = score_dictionary

    kmer_score = 0
    size = min(len(pep_seq), 5)

    if size <= 2:
        gap = 0
    else:
        gap = size - 2

    kmers = find_kmer(sequence = pep_seq, kmer_size = size, n_gap = gap, r_dict = REDUCE)
    for kmer in kmers:
        if kmer in score_dict.keys():
            kmer_score += score_dict[kmer][2]

    return kmer_score / len(pep_seq)


# Computation of peptide physical properties
def pep_physical_analysis(pep_seq) -> list[float]:
    """
    pep_seq: peptide sequence
    """
    pa = ProteinAnalysis(str(pep_seq))
    helix_propensity = pa.secondary_structure_fraction()[0]
    charge = pa.charge_at_pH(7)

    """
    Kyte-Doolitlle hydrophobicity profile
    """

    hydrophobicity_profile = pa.protein_scale(ProtParamData.kd, 2, edge=1.0)

    # Identify the minimal distance between 2 peaks
    peaks, _ = find_peaks(hydrophobicity_profile, distance=2)

    # Mean distance between hydrophilic maximums
    periodicity = np.mean(np.diff(peaks))
    return [pep_seq, helix_propensity * 100, charge, periodicity]

def generate_prob_dict_from_excel(file_name: str = PAM2_EXCEL_FILE) -> dict:
    pam2_prob_matrix = pd.read_excel(file_name)

    # Convert the DataFrame to a dictionary
    pam2_prob_dict: dict = {}
    for row in pam2_prob_matrix.to_numpy():
        pam2_prob_dict[row[0]] = row[1:]
    return pam2_prob_dict

def generate_peptides(
    aa_probs: dict, nb_peptide: int = NB_PEPTIDE, default_peptide: str = DEFAULT_PEPTIDE, bootstrap: int = NB_ITERATIONS
):
    Align(rich_print(Panel(f"Generation of {nb_peptide} peptides", style = "light_slate_blue", expand = False)))
    # Test
    score_evolution = []
    helix_propensity_evolution = []
    charge_evolution = []
    periodicity_evolution = []
    peptides_generated = []

    for p in range(0, nb_peptide):
        pep_seq = default_peptide
        with default_progress() as progress:
            task_id: TaskID = progress.add_task(f"Generating peptide number {p+1} out of {nb_peptide}", total=bootstrap)
            for i in range(bootstrap):
                # Randomisation of mutation location in the peptide sequence should be applied to biological form (To develop)
                random_index = random.randint(0, len(pep_seq) - 1)

                # Replacing the amino acid selected to a knew one
                random_amino_acid = pep_seq[random_index]
                prob = aa_probs[random_amino_acid]
                new_amino_acid = random.choices(AA_ORDER)[0]
                new_peptide = pep_seq[:random_index] + new_amino_acid + pep_seq[random_index + 1 :]

                # Calculating scores of previous and new peptides sequences
                peptide_score = score_kmers(pep_seq, REDUCE, score_dictionary)
                physical_analysis: list[float] = pep_physical_analysis(pep_seq)
                new_peptide_score = score_kmers(new_peptide, REDUCE, score_dictionary)

                # The peptide is selected if new score is higher
                score_difference = new_peptide_score - peptide_score

                if score_difference > 0:
                    pep_seq = new_peptide

                # Progress bar
                progress.update(task_id, advance = 1)

        table = Table(
            title = f"[i]Result of peptide[/] [b royal_blue1]{pep_seq}",
            show_header = False,
            box = box.ROUNDED,
            style = "cyan",
            title_style = "",
        )
        table.add_row("Final score", f"{score_kmers(pep_seq,REDUCE,score_dictionary)}")
        table.add_row("Final helix probability", f"{round(physical_analysis[1], 2)}%")
        table.add_row("Final global charge Q", f"{physical_analysis[2]}")
        table.add_row("Final hydrophobicity frequency", f"{physical_analysis[3]}")

        rich_print("\n", table, "\n")

        score_evolution.append(new_peptide_score)
        helix_propensity_evolution.append(physical_analysis[1])
        charge_evolution.append(physical_analysis[2])
        periodicity_evolution.append(physical_analysis[3])
        peptides_generated.append(pep_seq)

    final_df = pd.DataFrame(
        list(zip(peptides_generated, score_evolution, charge_evolution, periodicity_evolution)),
        columns=["peptide_sequence", "Activity score", "Net charge", "Hydrophobicity-based periodicity"],
    )
    rich_print(
        Align(
            Panel("Generated peptides and their respectives properties", style = "light_slate_blue", expand = False),
            align="center",
        )
    )
    print(final_df)
    final_df.to_excel("results/de_novo_peptide_library.xlsx")

    '''
    Evolution plots
    
    plt.plot(NB_ITERATIONS, score_evolution, label = "score")
    plt.plot(NB_ITERATIONS, helix_propensity_evolution, label = "helix propensity")
    plt.plot(NB_ITERATIONS, charge_evolution, label = "charge")
    plt.plot(NB_ITERATIONS, periodicity_evolution, label = "hydrophobicity space")
    plt.legend()
    plt.show()
    '''

if __name__ == "__main__":
    # Getting scores from CSV file
    score_dictionary = load_descriptors_score()

    # Import substitution probabilities from PAM2 data frame relative to mutation frequencies with conservation excluded
    pam2_probs = generate_prob_dict_from_excel()

    ### generation and optimisation of a peptide sequence
    generate_peptides(pam2_probs)




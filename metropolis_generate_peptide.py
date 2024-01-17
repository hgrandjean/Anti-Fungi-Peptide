import random
import math
import os
import shutil
import sys

import numpy as np
import pandas as pd
import seaborn as sns
from Bio.SeqUtils import (
    ProtParamData,  # Local https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/ProtParamData.py
)
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

from kmer_parser import find_kmer, parse_fasta_file

from rich import box, print as rich_print
from rich.panel import Panel
from rich.progress import BarColumn, Progress, SpinnerColumn, MofNCompleteColumn, TaskID, TextColumn, TimeElapsedColumn
from rich.table import Table
from rich.align import Align
from rich.prompt import Confirm

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

# inputs and outputs
SCORE_FILE = "results/descriptors_activity_scores.tsv"
PAM2_EXCEL_FILE = "resources/PAM_2_substitution_probabilities_formated.xlsx"
POSITIVE_DB_FILE = "resources/filtered_pos_db.fasta"
output_file = "results/generated_peptides.fasta"
tango_dir = "".join(os.getcwd() + "/results/tango_results/")
tango_output = "".join(tango_dir + "generated_peptides_tango.sh")
AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
DEFAULT_PEPTIDES = ["RGLRRLGRKIAHGVKKYG"]
NB_PEPTIDE = 50
NB_ITERATIONS = 50
NB_BOOTSTRAPS = 500
KB = 0.001987204  # Boltzmann constant
T = 1500  # Temperature for Metropolis algorithm


if len(sys.argv) > 1:
    NB_PEPTIDE = int(sys.argv[1])
    NB_ITERATIONS = int(sys.argv[2])

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
        transient=False,
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


def setup_directory(dir_name):
    if os.path.exists(dir_name):
        answer = input(f"Found {dir_name}\nAre you sure that you want to delete it? [y, n]\n")
        if answer == "y":
            shutil.rmtree(dir_name)
            print(f"{dir_name} deleted.")
        else:
            print("Operation canceled")
            os._exit(1)

    os.makedirs(dir_name)
    print(f"Created {dir_name}")


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

    kmers = find_kmer(sequence=pep_seq, kmer_size=size, n_gap=gap, r_dict=REDUCE)
    for kmer in kmers:
        if kmer in score_dict.keys():
            kmer_score += score_dict[kmer][2]

    return kmer_score


# Computation of peptide physical properties
def calculate_moment(array, angle=100):
    """Calculates the hydrophobic dipole moment from an array of hydrophobicity
    values. Formula defined by Eisenberg, 1982 (Nature). Returns the average
    moment (normalized by sequence length)

    uH = sqrt(sum(Hi cos(i*d))**2 + sum(Hi sin(i*d))**2),
    where i is the amino acid index and d (delta) is an angular value in
    degrees (100 for alpha-helix, 180 for beta-sheet).

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


### Computation of peptide's physical properties
def pep_physical_analysis(pep_seq: str) -> list[float]:
    """
    pep_seq: peptide sequence
    """
    pa = ProteinAnalysis(pep_seq)
    helix_propensity = pa.secondary_structure_fraction()[0]
    charge = pa.charge_at_pH(7)

    """
    Kyte-Doolitle hydrophobicity profile
    """

    hydrophobicity_profile = pa.protein_scale(ProtParamData.kd, 2, edge=1.0)

    # Identify the minimal distance between 2 peaks
    peaks, _ = find_peaks(hydrophobicity_profile, distance=2)

    # Mean distance between hydrophilic maximums
    periodicity = np.mean(np.diff(peaks))
    return [pep_seq, helix_propensity * 100, charge, periodicity, calculate_moment(hydrophobicity_profile), pa.gravy()]


"""
def generate_prob_dict_from_excel(file_name: str = PAM2_EXCEL_FILE) -> dict:
    pam2_prob_matrix = pd.read_excel(file_name)

    # Convert the DataFrame to a dictionary
    pam2_prob_dict: dict = {}
    for row in pam2_prob_matrix.to_numpy():
        pam2_prob_dict[row[0]] = row[1:]
    return pam2_prob_dict
"""


def generate_fasta_file(sequences, names, output_file):
    with open(output_file, "w") as fasta_file:
        for i in range(len(sequences)):
            fasta_file.write(f">{names[i]}\n")
            fasta_file.write(f"{sequences[i]}\n")


def generate_tango_script(peptides: str, names: str, output_file):
    setup_directory(tango_dir)
    with open(output_file, "w") as file:
        for peptide, name in zip(peptides, names):
            line = f'./Tango P{name}_ nt="N" ct="N" ph="7.4" te="303" io="0.05" tf="0" stab="-4" seq="{peptide}" >> results/tango_results/peptide_agregg.txt\n'
            file.write(line)


def keep_best_peptides(peptides: dict[str, float], nb_best_pep: int = NB_PEPTIDE) -> dict[str, float]:
    return dict(sorted(peptides.items(), key=lambda x: x[1], reverse=True)[:nb_best_pep])


def get_first_peptides(file_name=POSITIVE_DB_FILE) -> dict[str, float]:
    best_peptides: dict[str, float] = {}
    peptides = parse_fasta_file(file_name, True)
    for item in peptides:
        best_peptides[str(item.seq)] = score_kmers(item.seq, REDUCE, score_dictionary)
    return keep_best_peptides(best_peptides)


def transform_dataframe(df, nb_peptides: int = NB_PEPTIDE):
    transformed_df = pd.DataFrame(columns=["Iterations", "Values", "Global descriptors"])
    for col in df.columns:
        variable_name = col
        values = df[col].values.tolist()
        print(df[col].values)
        x_values = list(range(NB_BOOTSTRAPS)) * nb_peptides
        print(len(x_values), len(values))
        temp_df = pd.DataFrame({"Iterations": x_values, "Values": values, "Global descriptors": variable_name})
        transformed_df = pd.concat([transformed_df, temp_df], ignore_index=True)
    return transformed_df


def generate_peptides(nb_peptide: int = NB_PEPTIDE, bootstrap: int = NB_BOOTSTRAPS, iterations: int = NB_ITERATIONS):
    Align(
        rich_print(
            Panel(
                f"Generation of {nb_peptide} peptides over {iterations} iterations",
                style="light_slate_blue",
                expand=False,
            )
        )
    )
    # Test
    score_evolution: list[float] = []
    helix_propensity_evolution: list[float] = []
    charge_evolution: list[float] = []
    periodicity_evolution: list[float] = []
    hydrophobic_moment_evolution: list[float] = []
    gravy_evolution: list[float] = []
    peptides_generated: list[str] = []
    best_peptides: dict[str, float] = {}
    all_old_pep = []
    metropolis: bool = False
    selected = 0
    rejected = 0
    all_best_pep = []

    # Variables for convergence plot
    plot_score_evolution: list[float] = []
    plot_helix_propensity_evolution: list[float] = []
    plot_charge_evolution: list[float] = []
    plot_periodicity_evolution: list[float] = []
    plot_hydrophobic_moment_evolution: list[float] = []
    plot_gravy_evolution: list[float] = []
    plot_peptides_generated: list[str] = []

    # get first pep
    best_peptides = get_first_peptides()
    all_old_pep.append(best_peptides)

    for i in range(iterations):
        rich_print(
            Align(
                Panel(f"Beginning of iteration {i+1}", style="light_slate_blue", expand=False),
                align="center",
            )
        )

        tmp_best: dict[str, float] = {}
        for index, pep_seq in enumerate(best_peptides.keys()):
            with default_progress() as progress:
                task_id: TaskID = progress.add_task(
                    f"Generating peptide number {index+1} out of {len(best_peptides)} | Iteration {i+1}",
                    total=bootstrap,
                )
                for _ in range(bootstrap):
                    # Randomisation of mutation location in the peptide sequence should be applied to biological form (To develop)
                    random_index = random.randint(0, len(pep_seq) - 1)

                    # Replacing the amino acid selected to a knew one
                    new_amino_acid = random.choices(AA_ORDER)[0]
                    new_peptide = pep_seq[:random_index] + new_amino_acid + pep_seq[random_index + 1 :]

                    # Calculating scores of previous and new peptides sequences
                    peptide_score = score_kmers(pep_seq, REDUCE, score_dictionary)
                    best_physical_analysis: list[float] = pep_physical_analysis(pep_seq)
                    new_peptide_score = score_kmers(new_peptide, REDUCE, score_dictionary)

                    # Plot evolution of scores
                    if i == iterations - 1:
                        plot_score_evolution.append(peptide_score)
                        plot_helix_propensity_evolution.append(best_physical_analysis[1])
                        plot_charge_evolution.append(best_physical_analysis[2])
                        plot_periodicity_evolution.append(best_physical_analysis[3])
                        plot_hydrophobic_moment_evolution.append(best_physical_analysis[4])
                        plot_gravy_evolution.append(best_physical_analysis[5])
                        plot_peptides_generated.append(pep_seq)

                    # The peptide is selected if new score is higher
                    score_difference = new_peptide_score - peptide_score

                    if score_difference >= 0:
                        pep_seq = new_peptide
                    else:
                        metropolis = True
                        a = random.random()
                        if math.exp(score_difference / (KB * T)) > a:
                            pep_seq = new_peptide
                            selected += 1
                        else:
                            rejected += 1
                    # Progress bar
                    progress.update(task_id, advance=1)

            table = Table(
                title=f"[i]Result of peptide[/] [b royal_blue1]{pep_seq}",
                show_header=False,
                box=box.ROUNDED,
                style="cyan",
                title_style="",
            )
            table.add_row("Final score", f"{score_kmers(pep_seq,REDUCE,score_dictionary)}")
            table.add_row("Final helix probability", f"{round(best_physical_analysis[1], 2)}%")
            table.add_row("Final global charge Q", f"{best_physical_analysis[2]}")
            table.add_row("Final hydrophobicity frequency", f"{best_physical_analysis[3]}")
            table.add_row("Final hydrophobicity moment", f"{best_physical_analysis[4]}")
            table.add_row("Final Average hydrophobicity", f"{best_physical_analysis[5]}")

            rich_print("\n", table, "\n")

            tmp_best[pep_seq] = new_peptide_score

        best_peptides = keep_best_peptides(best_peptides | tmp_best)
    all_best_pep.append(best_peptides)

    for pep_seq, best_score in best_peptides.items():
        best_physical_analysis: list[float] = pep_physical_analysis(pep_seq)

        score_evolution.append(best_score)
        helix_propensity_evolution.append(best_physical_analysis[1])
        charge_evolution.append(best_physical_analysis[2])
        periodicity_evolution.append(best_physical_analysis[3])
        hydrophobic_moment_evolution.append(best_physical_analysis[4])
        gravy_evolution.append(best_physical_analysis[5])
        peptides_generated.append(pep_seq)

    rich_print(all_best_pep, all_old_pep)
    rich_print(f"Rejected: {rejected}, Selected: {selected}")

    final_df = pd.DataFrame(
        list(
            zip(
                peptides_generated,
                score_evolution,
                charge_evolution,
                periodicity_evolution,
                helix_propensity_evolution,
                gravy_evolution,
                hydrophobic_moment_evolution,
            )
        ),
        columns=[
            "peptide_sequence",
            "score",
            "net_charge",
            "hydrophobicity_periodicity",
            "helix_propensity",
            "hydrophobicity_average",
            "hydrophobic_moment",
        ],
    )

    df_best_histplot = pd.DataFrame(
        list(zip(all_best_pep[0].keys(), all_best_pep[0].values())), columns=["peptide", "Activity score"]
    )
    df_old_histplot = pd.DataFrame(
        list(zip(all_old_pep[0].keys(), all_old_pep[0].values())), columns=["peptide", "Activity score"]
    )

    rich_print(
        Align(
            Panel("Generated peptides and their respectives properties", style="light_slate_blue", expand=False),
            align="center",
        )
    )
    print(final_df)
    df_best_histplot.to_excel("results/de_novo_best_peptide.xlsx")
    df_old_histplot.to_excel("results/de_novo_old_peptide.xlsx")
    final_df.to_excel("results/de_novo_peptide_metropolis.xlsx")
    if metropolis:
        print(f"Metropolis selection: {int(selected*100/(selected+rejected))}%")
    print("Run in vivo aggregation study at http://bioinf.uab.es/aggrescan/ using generated fasta file in results/")
    print("Run in vitro aggregation study using Tango algorithm using generated bash file {}".format(tango_output))

    # Generate the FASTA file
    # generate_fasta_file(best_peptides.values(), range(0, nb_peptide), output_file)
    # generate_tango_script(best_peptides.keys(), range(0, nb_peptide), tango_output)

    """
    Convergence plots
    """
    plot_df = pd.DataFrame(
        list(
            zip(
                plot_score_evolution,
                plot_charge_evolution,
                plot_periodicity_evolution,
                plot_helix_propensity_evolution,
                plot_gravy_evolution,
                plot_hydrophobic_moment_evolution,
            )
        ),
        columns=[
            "Score",
            "Charge",
            "Periodicity",
            "Helix propensity",
            "Hydrophobicity average",
            "Hydrophobic moment",
        ],
    )

    answer = Confirm.ask("Do you want to generate a convergence plot?")
    if answer:
        df = transform_dataframe(plot_df)
        print(df)

        plot = sns.lineplot(x="Iterations", y="Values", hue="Global descriptors", data=df)
        sns.move_legend(
            plot,
            "lower center",
            bbox_to_anchor=(0.5, 1),
            ncol=3,
            title=None,
            frameon=False,
        )
        plt.savefig("results/metropolis_convergence_plot.png")


if __name__ == "__main__":
    # Getting scores from CSV file
    score_dictionary = load_descriptors_score()

    # Import substitution probabilities from PAM2 data frame relative to mutation frequencies with conservation excluded
    # pam2_probs = generate_prob_dict_from_excel()

    ### generate a directory for Tango outputs
    # setup_directory(tango_dir)

    ### generation and optimisation of a peptide sequence
    generate_peptides()

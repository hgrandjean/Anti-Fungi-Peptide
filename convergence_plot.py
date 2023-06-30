import random
import os
import shutil 

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from kmer_parser import setup_directory

from rich import box, print as rich_print
from rich.panel import Panel
from rich.progress import BarColumn, Progress, SpinnerColumn, MofNCompleteColumn, TaskID, TextColumn, TimeElapsedColumn
from rich.table import Table
from rich.align import Align

from generate_peptide import default_progress
from generate_peptide import load_descriptors_score, score_kmers, pep_physical_analysis
from generate_peptide import generate_prob_dict_from_excel, generate_tango_script, generate_fasta_file


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
    " #@@@@@@@@@#-#@@#=   +@@=@@@@@@@@@@% @@@ *@@@@@@@@@@=*@@+    +@@#=@@@    @@@ #@@@@@@@@@#=+%@@@@@@@@%+=@@#          =+++++++++*******+++=             "
)
print(
    " @@@=====%@@-#@@@#=  +@@-===+@@%==== @@@ %@@+=======-*@@+    +@@#=@@@%   @@@ @@@+====#%%=#@@*====*@@*=@@#           +=======++***+++++++             "
)
print(
    " @@@     %@@-#@@@@#- +@@*    -@@#     @@@ %@@-        @@+    +@@#=@@@@%  @@@ @@@-        #@@+    +@@*=@@#                  -=+++++++++++=            "
)
print(
    " @@@@@@@@@@@-#@@#@#=+@@    -@@#     @@@ %@@@@@@@=   *@@+    +@@#=@@%+%@*-@@@ @@@-  #@@@#-#@@@@@@@@@@*=@@#               =+***+++++++++++             "
)
print(
    " @@@+++++%@@-#@@= @@@@@    -@@#     @@@ %@@+++++-   *@@+    +@@#=@@% -@@@@@@ @@@-  =+%@@=#@@#++++#@@*=@@#            =+******++++++++=-              "
)
print(
    " @@@     %@@-#@@=  @@@@    -@@#     @@@ %@@-        *@@+    +@@#=@@%   @@@@@ @@@-    #@@=#@@+    +@@*=@@#         -+**********+++++-                 "
)
print(
    " @@@     %@@-#@@=   @@@    -@@#     @@@ %@@-        *@@@@@@@@@@#=@@%    @@@@ @@@@@@@@@@@=#@@+    +@@*=@@@@@@@@@@@-##***********+-                    "
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
    "                                        -@%+++#@#@*+++  -@%+++#@   %@-  =@%+@*   %@=%@*++=  +@#++++=                            =**++++++++++=       "
)
print(
    "                                        -@@****-#@#**+  -@%*****-   %@-  =@%+@   %@=%@***=   *****%@=                        =+***+++++++++++        "
)
print(
    "                                        -@%      #@*=====-@%         %@-  =@%+@#===%@=%@+======#*===%@=                     =+*******+++++++=-         "
)
print(
    "                                        -#*      -#######-#*         #   -#*=######*-+######*-*#####                    +***********++++-            "
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
output_file = 'results/generated_peptides.fasta'
tango_dir="".join(os.getcwd() +'/tango_results/')
tango_output="".join(tango_dir +'generated_peptides_tango.sh')
AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
DEFAULT_PEPTIDE = "RGLRRLGRKIAHGVKKYG"
NB_PEPTIDE = 20
NB_ITERATIONS = 1000

# Selected reduction dictionary
REDUCE = 6

def generate_peptides(
    aa_probs: dict, NB_PEPTIDE: int = NB_PEPTIDE, DEFAULT_PEPTIDE: str = DEFAULT_PEPTIDE, bootstrap: int = NB_ITERATIONS
):
    Align(rich_print(Panel(f"Generation of {NB_PEPTIDE} peptides", style ="light_slate_blue", expand = False)))
    # Test
    boot: list[int] = []
    score_evolution: list[float] = []
    helix_propensity_evolution: list[float] = []
    charge_evolution: list[float] = []
    periodicity_evolution: list[float] = []
    hydrophobic_moment_evolution: list[float] = []
    gravy_evolution: list[float] = []
    peptides_generated: list[str] = []


    for p in range(0, NB_PEPTIDE):
        pep_seq = DEFAULT_PEPTIDE
        with default_progress() as progress:
            task_id: TaskID = progress.add_task(f"Generating peptide number {p+1} out of {NB_PEPTIDE}", total=bootstrap)
            for i in range(bootstrap):
            
                # Randomisation of mutation location in the peptide sequence should be applied to biological form (To develop)
                random_index = random.randint(0, len(pep_seq) - 1)

                # Replacing the amino acid selected to a new one
                random_amino_acid = pep_seq[random_index]
                prob = aa_probs[random_amino_acid]
                new_amino_acid = random.choices(AA_ORDER)[0]
                new_peptide = pep_seq[:random_index] + new_amino_acid + pep_seq[random_index + 1 :]

                # Calculating scores of previous and new peptides sequences
                peptide_score = score_kmers(pep_seq, REDUCE, score_dictionary)
                physical_analysis: list[float] = pep_physical_analysis(pep_seq)
                new_peptide_score = score_kmers(new_peptide, REDUCE, score_dictionary)
               
                # Plot evolution of scores 
                boot.append(i)
                score_evolution.append(peptide_score)
                helix_propensity_evolution.append(physical_analysis[1])
                charge_evolution.append(physical_analysis[2])
                periodicity_evolution.append(physical_analysis[3])
                hydrophobic_moment_evolution.append(physical_analysis[4])
                gravy_evolution.append(physical_analysis[5])
                peptides_generated.append(pep_seq)

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
        table.add_row("Final hydrophobicity moment", f"{physical_analysis[4]}")
        table.add_row("Final Average hydrophobicity", f"{physical_analysis[5]}")

        rich_print("\n", table, "\n")

        

    final_df=pd.DataFrame(list(zip(peptides_generated, score_evolution, charge_evolution, periodicity_evolution , helix_propensity_evolution ,gravy_evolution , hydrophobic_moment_evolution)), columns = ["peptide_sequence", "activity_score" , "net_charge" , "hydrophobicity_periodicity" ,"helix_propensity","hydrophobicity_average", "hydrophobic_moment"])
    
    rich_print(
        Align(
            Panel("Generated peptides and their respective properties", style = "light_slate_blue", expand = False),
            align="center",
        )
    )
    print(final_df)
    final_df.to_excel("results/de_novo_peptide_library.xlsx")
    print("Run in vivo aggregation study at http://bioinf.uab.es/aggrescan/ using generated fasta file in results/")
    print("Run in vitro aggregation study using Tango algorithm using generated bash file {}".format(tango_output))
    
    # Generate the FASTA file
    #generate_fasta_file(peptides_generated, range(0, nb_peptide), output_file)
    #generate_tango_script(peptides_generated, range(0, nb_peptide), tango_output)
    
    '''
    Convergence plots
    '''
    
    #define dataframe 
    df= {
        "Score" : score_evolution , 
        "Charge": charge_evolution , 
        "Periodicity" : periodicity_evolution ,
        "Hydrophobic moment" : hydrophobic_moment_evolution,
        "Hydrophobicity average" : gravy_evolution,
        "Helix propensity" : helix_propensity_evolution }
    df = pd.DataFrame(df)
    
    def transform_dataframe(df):
        transformed_df = pd.DataFrame(columns=['Iterations', 'Values', 'Global descriptors'])
        for col in df.columns:
            variable_name = col
            values = df[col].values.tolist()
            x_values = list(range(NB_ITERATIONS))*NB_PEPTIDE
            temp_df = pd.DataFrame({'Iterations': x_values, 'Values': values, 'Global descriptors': variable_name})
            transformed_df = pd.concat([transformed_df, temp_df], ignore_index=True)
        return transformed_df
    
    
    df = transform_dataframe(df)
        
    plot = sns.lineplot(x="Iterations", y="Values",
                hue="Global descriptors",
                data=df)
    sns.move_legend(
        plot, "lower center",
        bbox_to_anchor=(.5, 1), ncol=3, title=None, frameon=False,
    )
    plt.savefig("results/Convergence_plot.png")
    







if __name__ == "__main__":
    # Getting scores from CSV file
    score_dictionary = load_descriptors_score()

    # Import substitution probabilities from PAM2 data frame relative to mutation frequencies with conservation excluded
    pam2_probs = generate_prob_dict_from_excel()
    

    ### generation and optimisation of a peptide sequence
    generate_peptides(pam2_probs)
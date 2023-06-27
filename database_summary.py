from collections import Counter, defaultdict

import matplotlib.pyplot as plt

from kmer_parser import parse_fasta_file

# Define the color properties and databases
PROPERTIES_AA = [["D", "E"], ["K", "R", "H"], ["N", "Q", "S","T", "Y"], ["F", "C", "W", "A", "V", "L", "I", "G", "M"], ["P"]]
COLOR_AA = ["grey", "blue", "green", "orange", "red"]

"""
The categories of the AA are set in accordance with the used reduction dictionary, with the exception of P, which is considered as a unique AA

Grey: hydrophilic, negatively charged
Blue: hydrophilic, positively charged
Green: hydrophilic, uncharged
Orange: hydrophobic
Red: proline 
"""

POS_DB_NAME = "resources/filtered_positive_db.fasta"
NEG_DB_NAME = "resources/filtered_negative_db.fasta"

# Define the graph
fig = plt.figure()
ax = fig.add_subplot()
fig.subplots_adjust(top=0.5)

# Get fasta file, sort sequences and plot
def sort_aa_positions(db_file: str):
    db_fastas = parse_fasta_file(db_file)
    print("Parsing fasta files")
    position_counts = defaultdict(list)
    for fasta in db_fastas:
        seq = fasta.seq
        for position, aa in enumerate(list(seq)):
            position_counts[position].append(aa)
    count = {k: Counter(v) for k, v in position_counts.items()}

    font_size = 0
    color_default = "black"
    for fasta in db_fastas:
        seq = fasta.seq
        for aa in seq:
            for position_aa in range(len(count)):
                if aa in count[position_aa]:
                    counted_value = count[position_aa][aa]
                    font_size = (counted_value / 40) * 35
                    for color in range(len(PROPERTIES_AA)):
                        if aa in PROPERTIES_AA[color]:
                            color_default = COLOR_AA[color]
                    ax.text(position_aa, y_pos, aa, color=color_default, fontsize=font_size)
                    del count[position_aa][aa]


if __name__ == "__main__":
    # Set titles for the figure and the subplot respectively
    ax.axis([0, 18, -1, 1])
    ax.text(-1, 1.1, "Summary of databases with amino acid occurences and positions", fontsize=10, fontweight="bold")
    ax.set_xlabel("Amino Acid positions")
    ax.set_ylabel("Sizes in terms of occurrence")

    # Generate the plot
    y_pos = 0.05
    sort_aa_positions(POS_DB_NAME)
    ax.text(5, 0.5, "Positive database", color="black", fontsize=10)
    y_pos = -0.15
    sort_aa_positions(NEG_DB_NAME)
    ax.text(5, -0.5, "Negative database", color="black", fontsize=10)
    plt.show()

import matplotlib
import matplotlib.pyplot as plt
from Bio import SeqIO
from collections import defaultdict, Counter

#Define the color properties and databases
properties_aa = [["D", "E", "K", "N"], ["H", "F", "Y", "W"], ["A", "V", "L", "I", "G"]]
color_aa = ["red", "blue", "green"]
pos_db_name = "resources/filtered_positive_db.fasta"
neg_db_name = "resources/filtered_negative_db.fasta"

#Define the graph
fig = plt.figure()
ax = fig.add_subplot()
fig.subplots_adjust(top=0.5)

#Parse the fasta file
def parse_fasta_file (file_name) -> list[str]:
    multi_fasta = [record for record in SeqIO.parse(file_name, "fasta")]
    print(f"Ended parsing of {file_name}")
    return multi_fasta

#Set titles for the figure and the subplot, respectively
ax.axis([0, 18, -1, 1])
ax.text(-1, 1.1, 'Summary of databases with amino acid frequencies and positions', fontsize=10, fontweight='bold')
ax.set_xlabel('Amino Acid positions')
ax.set_ylabel('Sizes in terms of occurrence')

#Get fasta file, sort sequences and plot
def sort_aa_positions(db_file):
    db_fastas = parse_fasta_file(db_file)
    print ("Parsing fasta files")
    position_counts = defaultdict(list)
    for fasta in db_fastas:
        seq = fasta.seq
        for position, aa in enumerate(list(seq)):
            position_counts[position].append(aa)
    count = {k:Counter(v) for k, v in position_counts.items()}

    font_size = 0
    color_default = "black"
    for fasta in db_fastas:
        seq = fasta.seq
        for aa in seq:
            for position_aa in range(len(count)):
                if aa in count[position_aa]:
                    counted_value = count[position_aa][aa]
                    font_size = (counted_value/40)*35
                    for color in range(len(properties_aa)):
                        if aa in properties_aa[color]: color_default = color_aa[color]
                    ax.text(position_aa, y_pos, aa, color=color_default, fontsize=font_size)  
                    del count[position_aa][aa]
y_pos = 0.05
sort_aa_positions(pos_db_name)
ax.text(5, 0.5, "Positive database", color="black", fontsize=10)
y_pos = -0.15
sort_aa_positions(neg_db_name)
ax.text(5, -0.5, "Negative database", color="black", fontsize=10)
plt.show()
from Bio import SeqIO
from kmer_parser import parse_fasta_file

# Database to be filtered
full_neg_file_name = "resources/full_negative_db.fasta"
full_pos_file_name = "resources/full_positive_db.fasta"

# Used databases
filtered_pos_file_name = "resources/filtered_pos_db.fasta"
filtered_neg_file_name = "resources/filtered_neg_db.fasta"

def clean_database(db_file_name, clean_db_file_name):
    print(f"Cleaning {db_file_name} to keep peptides between 3 and 50")
    multi_fasta = parse_fasta_file(db_file_name)
    multi_fasta_size = []
    for fasta in multi_fasta:
        seq = fasta.seq
        fasta.description = ""
        if fasta.id.find("|") != -1:
            fasta.id = "".join(fasta.id.split("|")[1])
        if len(seq) in range(3, 51):
            multi_fasta_size.append(fasta)

    SeqIO.write(multi_fasta_size, clean_db_file_name, "fasta")
    print(f"Output clean database in {clean_db_file_name}")

if __name__ == "__main__":
    print("Start filtering the databases")

# Select peptides with a length between 3 and 50 aa from a positive and a negative DB
    clean_database(full_neg_file_name, filtered_neg_file_name)
    clean_database(full_pos_file_name, filtered_pos_file_name)

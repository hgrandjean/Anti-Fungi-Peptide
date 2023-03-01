import os
from Bio import SeqIO
from collections import defaultdict, Counter
from itertools import combinations_with_replacement
import pandas as pd
import numpy as np
import hashlib
from matplotlib import pyplot as plt
import multiprocessing as mp
from collections import Counter
import threading
import time
from functools import partial


# reduce AA sequence complexity using different set of in-vitro/silico properties
# Reduction Encoding : 
# RED1 : Hydrophobicity A= hydrophobic ; B = hydrophilic ;
# RED2 : Physico-chemical   A= hydrophobic ; B = hydrophilic ; C = Aromatic ; D = Polar ; E = Acidic ; F = Basic ; G = Ionizable ;
# RED3 : Solvent accessibility ; A = Low ; B = Medium ; C = High
# RED4 : Hydrophobicity and charge; A = hydrophobic ; B = Hydrophilic : C = Charged
# RED5 : Hydrophobicity and structure;  A = Hydrophilic ; B = Hydrophobic : C = Structural
# RED6 : Hydrophobicity size and charge; A = Large and hydrophobic; B = small hydrophobic ; P = positive hydrophilic ; U = uncharged hydrophilic ; N = negative hydrophilic

reduction_dictionaries = {
    'A' :['A','A','B','B','B', 'B'], #Alanine
    'C' :['B','G','A','A','A', 'B'], #Cysteine
    'D' :['B','E','C','C','A', 'N'], #Aspartic acid
    'E' : ['B','E','C','C','A', 'N'], #Glutamic acid
    'F' : ['B','C','A','A','A', 'B'], #Phenylalanine
    'G' : ['A','A','B','B','C', 'B'], #Glycine
    'H' : ['B','B','B','A','A', 'P'], #Histidine
    'I' : ['A','A','A','B','B', 'A'], #Isoleucine
    'K' : ['B','F','C','C','A', 'P'], #Lysine
    'L' : ['A','A','A','B','B', 'A'], #Leucine
    'M' : ['A','A','A','B','B', 'A'], #Methionine
    'N' : ['B','D','C','A','A', 'U'], #Asparagine
    'P' : ['B','B','C','A','C', 'B'], #Proline
    'Q' : ['B','D','C','A','A', 'U'], #Glutamine
    'R' : ['B','F','C','C','A', 'P'], #Arginine
    'S' : ['B','D','B','A','A', 'U'], #Serine
    'T' : ['B','D','B','A','A', 'U'], #Threonine
    'V' : ['A','A','A','B','B', 'A'], #Valine
    'W' : ['B','-','A','A','A', 'A'], #Tryptophan
    'Y' : ['B','G','A','A','A', 'U'], #Tyrosine
    
    'r' : ['B','F','C','C','A', 'P'], #Arginine
    'J' : ['B','F','C','C','A', 'P'], #un-usual amino-acid
}
reduce = 6 
neg_temp_path = "".join(os.getcwd() + "/kmr_neg_temp/")
pos_temp_path = "".join(os.getcwd() + "/kmr_pos_temp/")

def reduce_seq(sequence, RED_dict ,r_dict = reduction_dictionaries):
    """ transform sequence using AA characteristics in proteins:
    __ Args __ 
    sequence (Seq): AA sequence in single letter codification 
    r_dict (dict) : transformation dictionary in single letter codification
    
    __ Returns __ 
    reduced AA sequence using transformation dictionary
    """
    reduced_seq = ""
    for aa in sequence:
        if aa not in r_dict.keys():
            pass
        else : 
            reduced_seq += r_dict[aa][RED_dict - 1]
    return reduced_seq 


def hash_kmer(kmer):
    """
    Hash a k-mer using the SHA-256 algorithm
    Args:
        kmer (str): The k-mer to hash
    Returns:
        str: The hashed k-mer
    """
    hashed_kmer = hashlib.sha256(kmer.encode()).hexdigest()
    return hashed_kmer


def gap_kmer(kmers):
    k_gap = []
    for kmer in kmers : 
        for z in range(0,len(kmer)) :
            if kmer[z] != "_" :
                k_gap.append("".join(kmer[:z] + "_" + kmer[z+1 :]))
    return set(k_gap) 

def find_kmer(sequence, kmer_size, ngap, reduce ):
    kmers = []
    if reduce != None :
        sequence = reduce_seq(sequence, RED_dict = reduce)
    for i in range(len(sequence)):
        if i+ kmer_size <= len(sequence):
                kmers .append (sequence[i:i+kmer_size])
    
    current_kmers=kmers           
    for k in range(ngap):
        current_kmers = gap_kmer(current_kmers)
        kmers += current_kmers

    #return [hash_kmer(kmer) for kmer in kmers] 
    return kmers


def get_kmers(seq_record , reduce, path):
    record , seq = seq_record.id , seq_record.seq
    with open("".join(path+f"{record}.kmr") , "w" ) as save:
        size = min(len(seq), 5)
        if size <= 2 : gap = 0 
        else : gap = size - 2 
        kmers = find_kmer(sequence= seq , kmer_size= size , ngap=gap , reduce= reduce )
        for kmer in kmers : 
            save.write("".join(str(kmer + '\n')))
            
def setup_directory(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
        print(f"Created {dir_name}")

def parse_fasta_file(file_name):
    multi_fasta = [record for record in SeqIO.parse(file_name, "fasta")]
    print(f"Ended parsing of {file_name}")
    return multi_fasta

def run(fastas, name):
    print(f"[{name}] Performing Gapped k-mer count on {len(fastas)} sequences ; reduction = {reduce} )")
    pool = mp.Pool(processes=4)

    # map the analyze_sequence function to the sequences
    main = partial(get_kmers, reduce = reduce  , path = neg_temp_path)
    results = pool.map(main , fastas)

    # close the pool and wait for the worker processes to finish
    pool.close()
    pool.join() 

    print(f"[{name}] Finished running")  

if __name__ == '__main__' : 

    # Create directories for stocking descriptors 
    setup_directory(neg_temp_path)
    setup_directory(pos_temp_path)

    # Get list of fastas
    neg_fastas = parse_fasta_file("negative_db_size.fasta")
    pos_fastas = parse_fasta_file("positive_db_size.fasta")
    
    # Create descriptors for each peptide
    run(neg_fastas, "Negative peptides")
    run(pos_fastas, "Positive peptides")
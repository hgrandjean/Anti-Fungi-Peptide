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



# reduce AA sequence complexity using different set of in-vitro/sillico properties 
# Reduction Encoding : 
# RED1 : Hydrophobicity A= hydrophocic ; B = hydrophilic ; 
# RED2 : Physiochemical   A= hydrophocic ; B = hydrophilic ; C = Aromatic ; D = Polar ; E = Acidic ; F = Basic ; G = Ionizable ; 
# RED3 : Solvent accessibility ; A = Low ; B = Medium ; C = High
# RED4 : Hydrophobicity and charge; A = hydrophobic ; B = Hydrophilic : C = Charged
# RED5 : Hydrophobicity and structure;  A = Hydrophilic ; B = Hydrophobic : C = Structural
# RED6 : Hydrophobicity size and charge; A = Large and hyphobic; B = small hydrophobic ; P = positive hydrophilic ; U = unchuraged hydrophilic ; N = negative hydrophilic 


reduction_dictionnaries = {  
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
    'Y' : ['B','G','A','A','A', 'Y'], #Tyrosine
}

def reduce_seq(sequence, RED_dict ,r_dict = reduction_dict):
    """ transform sequence using AA characteristics in proteins:
    __ Args __ 
    sequence (Seq): AA sequence in single letter codification 
    r_dict (dict) : transformation dictionnary in single letter codifiaction 
    
    __ Returns __ 
    reduced AA sequence using transformation dictionnary 
    """
    reduced_seq = ""
    for aa in sequence:
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
                k_gap.append(kmer[:z] + "_" + kmer[z+1 :])
    return k_gap

def find_kmer(sequence, kmer_size, ngap, reduce = reduce):
    kmers = []
    if reduce =! None :
        sequence = reduce_seq(sequence, RED_dict = reduce)
    for i in range(len(sequence)):
        if i+ kmer_size <= len(sequence):
                kmers .append (sequence[i:i+kmer_size])
    for i in range(ngap):
        kmers =gap_kmer(kmers)
    return [hash_kmer(kmer) for kmer in kmers] 
    #return kmers



def get_kmers(seq_record , kmer_size , reduce , path):
    record , seq = seq_record.id , seq_record.seq
    record = "".join(record.split('|')[1]+"_"+record.split('|')[2])
    with open("".join(path+"{}.kmr".format(record)) , "w" ) as save:
        for size in kmer_size :
          if size <= 2 : gap= 0 
          else : gap = size - 2 
          kmers = find_kmer(sequence= seq , kmer_size= size , ngap=gap , reduce= reduce )
          for kmer in kmers : 
              save.write("".join(str(kmer + '\n')))
    save.close()



if __name__ == '__main__' : 
    #parameter
    k = range(3,11)    
    reduce = 6 

    folder_path = "".join(os.getcwd()+"/kmr_temp/")

    if not os.path.exists(folder_path):

        os.makedirs(folder_path)

    multi_fasta = [record for record in SeqIO.parse("<db/path.fasta>", "fasta")]
    print("Perfomring Gapped k-mer count on {} sequences | parameters (k-mer size ={}; reduction = {} )".format(len(multi_fasta), k , str(reduce)))
    
    pool = mp.Pool(processes=4)

    # map the analyze_sequence function to the sequences
    main = partial(get_kmers, kmer_size = k , reduce = reduce  , path = folder_path)
    results = pool.map(main , multi_fasta)

    # close the pool and wait for the worker processes to finish
    pool.close()
    pool.join()




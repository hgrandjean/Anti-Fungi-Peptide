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
    'A' :['A','A','B','B','B'] ,
    'C' :['B','G','A','A','A'] ,
    'D' :['B','E','C','C','A'] ,
    'E' : ['B','E','C','C','A'] ,
    'F' : ['B','C','A','A','A'] ,
    'G' : ['A','A','B','B','C'] ,
    'H' : ['B','B','B','A','A'] ,
    'I' : ['A','A','A','B','B'] ,
    'K' : ['B','F','C','C','A'] ,
    'L' : ['A','A','A','B','B'] ,
    'M' : ['A','A','A','B','B'] ,
    'N' : ['B','D','C','A','A'] ,
    'P' : ['B','B','C','A','C'] ,
    'Q' : ['B','D','C','A','A'] ,
    'R' : ['B','F','C','C','A'] ,
    'S' : ['B','D','B','A','A'] ,
    'T' : ['B','D','B','A','A'] ,
    'V' : ['A','A','A','B','B'] ,
    'W' : ['B','-','A','A','A'] ,
    'Y' : ['B','G','A','A','A'] ,
   #unusual amino acids
    "X" : 'B' , 
    "U" : 'B' ,
    "B" : 'B' , 
    "O" : 'B' ,
    "X" : 'B' ,
    "Z" : 'B' ,
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

def find_kmer(sequence, kmer_size, ngap, reduce = False):
    kmers = []
    if reduce == True :
        sequence = reduce_seq(sequence)
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
          kmers = find_kmer(sequence= seq , kmer_size= size , ngap=gap , reduce= False )
          for kmer in kmers : 
              save.write("".join(str(kmer + '\n')))
    save.close()



if __name__ == '__main__' : 
    #parameter
    k = range(3,11)    
    reduce = False 

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




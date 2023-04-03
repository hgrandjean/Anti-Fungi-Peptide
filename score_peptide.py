from kmer_parser import reduce_seq , gap_kmer , find_kmer

def score_kmers(pep_seq , reduce, score_dict ):
      kmer_score = 0
      size = min(len(pep_seq), 5)
      if size <= 2 : gap = 0 
      else : gap = size - 2 
      kmers = find_kmer (sequence = pep_seq , kmer_size = size , ngap = gap , reduce = reduce )
      for kmer in kmers:
          if kmer in score_dict.keys() : 
              kmer_score += score_dict[kmer][2]
              # print(kmer," ", score_dict[kmer][2])
      return kmer_score
    
    
score_file = "unique_set.tsv"

### loadding score from computed tsv file ###
print (f"Loading descriptors scores from file : {score_file}")
score_dict = {}
with open(score_file, "r" ) as scores :
    for line in scores:
        key, value = line.removesuffix('\n').split('\t')
        value =  value.strip('][').split(', ')
        value = [float(x) for x in value ]
        score_dict[key] = value 
print ("Finished loading scores \n Starting scoring peptides")

score = score_kmers("RGLRRLGRKIAHGVKKYG ",6,score_dict) # SMAP-18
print(score)

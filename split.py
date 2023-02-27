

with open("kmr_temp/A2T3P0_NSP2_ROTSH.kmr" , "r") as pos : 
    positif = set (pos.readlines())
with open("kmr_temp/Q9FI78_HST_ARATH.kmr" , "r") as neg : 
    negatif = set (neg.readlines())

unique_positive_kmers=[]
for kmer in positif:
    if kmer not in negatif:
        unique_positive_kmers.append(kmer)

with open("unique_set.kmr" , "w") as save : 
    for kmer in unique_positive_kmers : 
              save.write("".join(str(kmer + '\n')))

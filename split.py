with open("kmr_temp/A2T3P0_NSP2_ROTSH.kmr" , "r") as pos : 
    positif = list (pos.readlines()) #replace set by something else 
with open("kmr_temp/Q9FI78_HST_ARATH.kmr" , "r") as neg : 
    negatif = list (neg.readlines()) #replace set by something else 

#count of descriptors in positive list 
kmers_counter={}
for kmer in positif:
    if kmer in kmers_counter.keys() : 
      kmers_counter[kmer][0] +=1 
    else : 
      kmers_counter[kmer] = [1,0]
for kmer in negatif:
    if kmer in kmers_counter.keys() : 
      kmers_counter[kmer][1] +=1 
    else : 
      kmers_counter[kmer] = [0,1]
          
#count of descriptors in negatif list 


with open("unique_set.tsv" , "w") as save : 
    for kmer in kmers_counter.keys() : 
      save.write(kmer+'\t'+str(kmers_counter[kmer])+'\n')
      

import math 
with open("kmr_temp/sp0.kmr" , "r") as pos : 
    positive = list (pos.readlines()) #replace set by something else
with open("kmr_temp/sp1.kmr" , "r") as neg : 
    negative = list (neg.readlines()) #replace set by something else

#count of descriptors in positive then in negative list
kmers_counter={}
for kmer in positive:
    if kmer in kmers_counter.keys() : 
      kmers_counter[kmer][0] +=1 
    else : 
      kmers_counter[kmer] = [1,0,0]
for kmer in negative:
    if kmer in kmers_counter.keys() : 
      kmers_counter[kmer][1] +=1 
    else : 
      kmers_counter[kmer] = [0,1, 0]

#score attribution to each descriptor
for kmer in kmers_counter.keys() : 
    kmers_counter[kmer][2] = math.log((kmers_counter[kmer][0]+1)/(kmers_counter[kmer][1]+1))

#save data to tsv file
with open("unique_set.tsv" , "w") as save : 
    for kmer in kmers_counter.keys() : 
      save.write(str(kmer).strip()+'\t'+str(kmers_counter[kmer])+'\n') 

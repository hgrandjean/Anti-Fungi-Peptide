import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import (
    ProtParamData,  # Local https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/ProtParamData.py
)
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from scipy.signal import find_peaks

multi_fasta = [record for record in SeqIO.parse("resources/filtered_positive_db.fasta", "fasta")]
space = np.array([])
loc = np.array([])

for seq in multi_fasta:
    pa = ProteinAnalysis(str(seq.seq))

    print(pa)
    print(pa.secondary_structure_fraction()[0])
    print(pa.charge_at_pH(7))
    print(pa.protein_scale(ProtParamData.kd, 2, edge=1.0))

    hydropho = pa.protein_scale(ProtParamData.kd, 2, edge=1.0)
    size = len(hydropho)

    # autocorrelation of signal
    """
    x = np.array(hydropho) 
    # Mean
    mean = np.mean(hydropho)
    # Variance
    var = np.var(hydropho)
    # Normalized data
    ndata = hydropho - mean
    acorr = np.correlate(ndata, ndata, 'full')[len(ndata)-1:] 
    acorr = acorr / var / len(ndata)
    """

    # get location of hydrophoilic residues
    peaks, _ = find_peaks(hydropho, distance=2)
    loc = np.append(loc, peaks)
    space = np.append(space, np.mean(np.diff(peaks)))  # space between hydrophilic residues

    # plt.plot(range(len(hydropho)), hydropho ) #
    # plt.plot(range(size), acorr )


loc = loc.flatten()
space = space.flatten()
print(loc)
print(space)
plt.hist(space)

# wheel =helixvis.draw_wheel(str(seq.seq))
# plt.title('Auto-correlation of hydrophobicity in the peptides sequences')
plt.show()
"""


for seq in seqs :  
  print('\n > ',seq)
  #code here 
  seq_ = seq+ "_"* (18 - len(seq))
  
  
  
  hydrophobicity = {"A": 0.62, "R": -2.53, "N": -0.78, "D": -0.9, "C": 0.29,
        "Q": -0.85, "E": -0.74, "G": 0.48, "H": -0.4, "I": 1.38,
        "L": 1.06, "K": -1.5, "M": 0.64, "F": 1.19, "P": 0.12,
        "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08 , "_": 0} #Eisenberg hydrophobicity score consensus computed from alpha helix and beta sheet 
    
  
  P1= np.array( [0, 4, 8 ,12, 16])
  P2= np.array( [1, 5, 9,13, 17])
  P3= np.array( [2, 6, 10,14])
  P4= np.array( [3, 7, 11,15])
  
  
  F1=  np.array([seq_[p] for p in P1])
  F2=  np.array([seq_[p] for p in P2])
  F3=  np.array([seq_[p] for p in P3])
  F4=  np.array([seq_[p] for p in P4])
  
  H1 = sum(hydrophobicity[aa] for aa in F1)
  H2 = sum(hydrophobicity[aa] for aa in F2)
  H3 = sum(hydrophobicity[aa] for aa in F3)
  H4 = sum(hydrophobicity[aa] for aa in F4)
  
 
#wheel =helixvis.draw_wheel(seq)
#plt.show()
 """

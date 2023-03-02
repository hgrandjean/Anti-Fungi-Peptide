from Bio.SeqUtils.ProtParam import ProteinAnalysis 
from Bio.SeqUtils import ProtParamData  # Local https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/ProtParamData.py
from Bio.SeqUtils import IsoelectricPoint  # Local
from Bio.Seq import Seq
from Bio.Data import IUPACData
from Bio.SeqUtils import molecular_weight
import helixvis
import matplotlib.pyplot as plt 
import numpy as np
import pylab as P
from scipy import signal


seq= "GILSSLWKKLKKIIAK"#"WLRAFRRLVRRLARLLRR"
pa = ProteinAnalysis(seq)
print(pa)
print(pa.secondary_structure_fraction()[0])
print(pa.charge_at_pH(7))
print(pa.protein_scale(ProtParamData.kd, 2, edge=1.0))
#wheel =helixvis.draw_wenxiang(seq)
hydropho = pa.protein_scale(ProtParamData.kd, 2, edge=1.0)
size = len(hydropho)
#autocorrelation of signal 
#x = np.array(hydropho) 
# Mean
#mean = np.mean(hydropho)
# Variance
#var = np.var(hydropho)
# Normalized data
#ndata = hydropho - mean
#acorr = np.correlate(ndata, ndata, 'full')[len(ndata)-1:] 
#acorr = acorr / var / len(ndata)
#plt.plot(range(len(ndata)), acorr )

corr = signal.correlate(hydropho, np.ones(size), mode='same') / size
plt.plot(corr)
plt.title('Cross-correlated hydrophobicity')
plt.show()
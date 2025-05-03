import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

###4/21 Ar 0.4 atm####
##inputfile:  C:\Data/MLF2024B/list/SDD_000006_01-00000000001.bin
##=====run0006 run6, 19MeV/c=====
##noise counts: 72.0
##11.12keV, nEv: 150, signal: 78.0, S/N: 2.0833333333333335
##20.48keV, nEv: 131, signal: 59.0, S/N: 1.8194444444444444
##inputfile:  C:\Data/MLF2024B/list/SDD_000007_01-00000000001.bin
##=====run0007 run7, 20MeV/c=====
##noise counts: 82.5
##11.12keV, nEv: 284, signal: 201.5, S/N: 3.4424242424242424
##20.48keV, nEv: 197, signal: 114.5, S/N: 2.3878787878787877
##inputfile:  C:\Data/MLF2024B/list/SDD_000008_01-00000000001.bin
##=====run0008 run8, 21MeV/c=====
##noise counts: 76.0
##11.12keV, nEv: 330, signal: 254.0, S/N: 4.342105263157895
##20.48keV, nEv: 217, signal: 141.0, S/N: 2.8552631578947367
##inputfile:  C:\Data/MLF2024B/list/SDD_000009_01-00000000001.bin
##=====run0009 run9, 22MeV/c=====
##noise counts: 94.0
##11.12keV, nEv: 331, signal: 237.0, S/N: 3.521276595744681
##20.48keV, nEv: 285, signal: 191.0, S/N: 3.0319148936170213
##inputfile:  C:\Data/MLF2024B/list/SDD_000010_01-00000000001.bin
##=====run0010 run10, 23MeV/c=====
##noise counts: 103.0
##11.12keV, nEv: 313, signal: 210.0, S/N: 3.0388349514563107
##20.48keV, nEv: 271, signal: 168.0, S/N: 2.6310679611650487

#mom=      [  19,   20,   21,   22,   23]
#sig_11keV=[78.0,201.5,254.0,237.0,210.0]
#sig_20keV=[59.0,114.5,141.0,191.0,168.0]   
#noi      =[72.0, 82.5, 76.0, 94.0,103.0]
#SN_11keV =[s/n for (s,n) in zip(sig_11keV,noi)]
#SN_20keV =[s/n for (s,n) in zip(sig_20keV,noi)]


###4/22 Ar 0.1 atm####
#Ar01atm_19MeV_20MeV_21MeV_22MeV_
#cal lines: [11.12, 20.48]
#inputfile:  C:\Data/MLF2024B/list/SDD_000026_01-00000000001.bin
#=====run0026 run26, 19MeV/c=====
#noise counts: 67.5
#11.12keV, nEv: 113, signal: 45.5, S/N: 1.674074074074074
#20.48keV, nEv: 102, signal: 34.5, S/N: 1.511111111111111
#inputfile:  C:\Data/MLF2024B/list/SDD_000027_01-00000000001.bin
#=====run0027 run27, 20MeV/c=====
#noise counts: 102.0
#11.12keV, nEv: 145, signal: 43.0, S/N: 1.4215686274509804
#20.48keV, nEv: 163, signal: 61.0, S/N: 1.5980392156862746
#inputfile:  C:\Data/MLF2024B/list/SDD_000028_01-00000000001.bin
#=====run0028 run28, 21MeV/c=====
#noise counts: 91.0
#11.12keV, nEv: 162, signal: 71.0, S/N: 1.7802197802197801
#20.48keV, nEv: 152, signal: 61.0, S/N: 1.6703296703296704
#inputfile:  C:\Data/MLF2024B/list/SDD_000029_01-00000000001.bin
#=====run0029 run29, 22MeV/c=====
#noise counts: 67.5
#11.12keV, nEv: 100, signal: 32.5, S/N: 1.4814814814814814
#20.48keV, nEv: 100, signal: 32.5, S/N: 1.4814814814814814
mom=      [  19,   20,   21,   22]
sig_11keV=[45.5, 43.0, 71.0, 32.5]
sig_20keV=[34.5, 61.0, 61.0, 32.5]   
noi      =[67.5,102.0, 91.0, 67.5]
SN_11keV =[s/n for (s,n) in zip(sig_11keV,noi)]
SN_20keV =[s/n for (s,n) in zip(sig_20keV,noi)]


fig,ax=plt.subplots(2,2,sharex='all')
ax[0][0].errorbar(mom,sig_11keV,yerr=np.sqrt(sig_11keV),color='k')
ax[0][1].errorbar(mom,sig_20keV,yerr=np.sqrt(sig_20keV),color='r')
ax[1][0].scatter(mom,SN_11keV,color='k')
ax[1][1].scatter(mom,SN_20keV,color='r')

ax[0][0].set_title('sig 11keV')
ax[0][1].set_title('sig 20keV')
ax[1][0].set_title('S/N 11keV')
ax[1][1].set_title('S/N 20keV')
ax[1][0].set_xlabel('Beam Mom [MeV/c]')
ax[1][1].set_xlabel('Beam Mom [MeV/c]')
plt.show()

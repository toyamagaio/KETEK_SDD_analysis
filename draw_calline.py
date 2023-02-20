import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scipy
from scipy.stats import norm
from scipy.optimize import curve_fit
import math
import pandas as pd

def draw_calline(inputfile = "cal_acgain20.csv", datadir="C:\Data/TMU2023A/calib/", outdir="C:\Data/TMU2023A/pdf/", outname="calib.pdf"):
  print('inputfile:',datadir+inputfile)
  pd_csv=pd.read_csv(datadir+inputfile, header=0,comment='#')
  print(pd_csv)
  print(pd_csv.columns)
  coeffs = np.polyfit(pd_csv["peakch"],pd_csv["energy"],1)
  a=coeffs[0]
  b=coeffs[1]
  x_space = np.linspace(min(pd_csv["peakch"])-100,max(pd_csv["peakch"])+100,100)
  #print(x_space)
  print('a: {0:.4e}, b: {1:.4e}'.format(a,b))

  #fig, ax=plt.subplots(2,1)
  #ax[0].scatter(pd_csv["peakch"],pd_csv["energy"])
  #ax[0].set_xlabel('pha ch')
  #ax[0].set_ylabel('Energy [keV]')

  #ax[0].plot(x_space, a*x_space+b, color='red',label='{0:.3e}x+{1:.3e}'.format(a,b))
  #ax[0].legend(ncol=1, fontsize=9, loc='upper left')

  fig=plt.figure(figsize=(8,6))
  ax1=fig.add_subplot(3,1,(1,2))
  ax1.scatter(pd_csv["peakch"],pd_csv["energy"])
  ax1.plot(x_space, a*x_space+b, color='red',label='{0:.3e}x+{1:.3e}'.format(a,b))
  ax1.legend(ncol=1, fontsize=9, loc='upper left')
  ax1.set_ylabel('Energy [keV]')
  ax2=fig.add_subplot(3,1,3, sharex=ax1)

  ax2.scatter(pd_csv["peakch"],pd_csv["energy"]-a*pd_csv["peakch"]+b)

  
  ax2.plot(x_space, np.zeros(len(x_space)), color='red',label='{0:.3e}x+{1:.3e}'.format(a,b))
  ax2.set_xlabel('ADC [ch]')
  ax2.set_ylabel('Residual [keV]')


  fig.tight_layout()
  fig.savefig(outdir+outname)

if __name__ == "__main__":
  #draw_calline(inputfile = "cal_acgain20.csv", datadir="C:\Data/TMU2023A/calib/", outdir="C:\Data/TMU2023A/pdf/", outname="calib.pdf")
  draw_calline(inputfile = "cal_acgain20_20230204.csv", datadir="C:\Data/TMU2023A/calib/", outdir="C:\Data/TMU2023A/pdf/", outname="calib_20230204.pdf")

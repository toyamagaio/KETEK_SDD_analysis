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

head_line=65

#def addcsv(inputfiles,outncv,datadir="C:\Data/TMU2022C/histograms/"):
#    for inputfile in inputfiles:
#      pd_csv=pd.read_csv(datadir+inputfile, header=head_line)
#      print(pd_csv)
#      print(pd_csv.columns)
#      print(pd_csv.columns[0])

def readcsv(inputfile = "list_000008_01-00000000001.bin", ENum = -1, datadir="C:\Data/TMU2022C/list/", outdir="C:\Data/TMU2022C/pdf/", outname="SDD_ene.pdf", cal_lines=[0], ranges=[[0,8000]]):
  print('inputfile: ', datadir+inputfile)
  pd_csv=pd.read_csv(datadir+inputfile, header=head_line)
  print(pd_csv)
  print(pd_csv.columns)
  print(pd_csv.columns[0])
  ch_or_ene=pd_csv.columns[0]

  if len(cal_lines) != len(ranges):
    print('inconsistent inputs, cal lines and ranges have different length.')
    return

  if ch_or_ene == 'ch':
      print('ch csv')
      xlabel="PHA [ch]"
      gene_cut = (pd_csv[ch_or_ene] > 500)
  elif ch_or_ene == 'eV':
      print('eV csv')
      xlabel="Energy [eV]"
      gene_cut = (pd_csv[ch_or_ene] >2000) & (pd_csv[ch_or_ene] < 100000)
  elif ch_or_ene == 'keV':
      print('keV csv')
      xlabel="Energy [keV]"
      gene_cut = (pd_csv[ch_or_ene] >2000) & (pd_csv[ch_or_ene] < 100000)
  else:
      print('unknown style')
  #print(pd_csv.values[0][0])
  print('max: {} '.format(max(pd_csv[ch_or_ene])))

  x_pos=[]
  peaks=[]
  for i, r in enumerate(ranges):
    cut = (pd_csv[ch_or_ene] > r[0]) & (pd_csv[ch_or_ene] < r[1])
    integral_in_r=(sum(pd_csv["CH1"][cut]))
    max_count   =(max(pd_csv["CH1"][cut]))
    idx_maxcount=np.argmax(pd_csv["CH1"][cut])
    aa =pd_csv[ch_or_ene][cut].values
    x_maxcount  =aa[idx_maxcount]
    print('range:',r[0],'--',r[1], xlabel, 'integral:',integral_in_r)
    #print(idx_maxcount, max_count, x_maxcount)
    print(cal_lines[i], x_maxcount)
    x_pos.append(x_maxcount)
    peaks.append(max_count)
  pp=PdfPages(outdir+outname)
  fig,ax=plt.subplots(2,1)
  ax[0].plot(pd_csv[ch_or_ene][gene_cut],pd_csv["CH1"][gene_cut])
  #ax[0].set_xlim(1800,41800)
  ax[0].scatter(x_pos,peaks,c='m')
  ax[1].scatter(x_pos,cal_lines)
  fig.tight_layout()

  fig2,ax2=plt.subplots(len(ranges),1)
  for i, r in enumerate(ranges):
    cut = (pd_csv[ch_or_ene] > r[0]) & (pd_csv[ch_or_ene] < r[1])
    ax2[i].plot(pd_csv[ch_or_ene][cut],pd_csv["CH1"][cut])
  fig2.tight_layout()

  pp.savefig(fig)
  pp.savefig(fig2)
  #plt.show()
  pp.close()

if __name__ == "__main__":
    datadir="C:\Data/TMU2023A/histograms/"

    #ranges=[[2000,3000], [3000,4000]]
    #readcsv(inputfile = "hist_230129_153628.csv", datadir=datadir, outdir="C:\Data/TMU2023A/pdf/", outname="SDD_ene_Eu152_test.pdf",  cal_source="Pd", ranges=ranges)

    #cal_lines=[39.522, 40.118, 45.294, 45.414] #keV
    #ranges=[[6800,6950],[6500,7500], [7600,8200], [8050,8200]]
    #readcsv(inputfile = "hist_230131_144931.csv", datadir=datadir, outdir="C:\Data/TMU2023A/pdf/", outname="SDD_cal_Eu152_test.pdf",  cal_lines=cal_lines, ranges=ranges)
    #readcsv(inputfile = "hist_230131_145350.csv", datadir=datadir, outdir="C:\Data/TMU2023A/pdf/", outname="SDD_ene_Eu152_test.pdf",  cal_source="Pd", ranges=ranges)
    #readcsv(inputfile = "hist_230204_235912.csv", datadir=datadir, outdir="C:\Data/TMU2023A/pdf/", outname="SDD_ene_Fe55_20230204.pdf",  cal_lines=cal_lines, ranges=ranges)


    #cal_lines=[14.41300] #keV
    #ranges=[[2000,3000]]
    #readcsv(inputfile = "hist_230129_132813.csv", datadir=datadir, outdir="C:\Data/TMU2023A/pdf/", outname="SDD_cal_Co57_test.pdf",  cal_lines=cal_lines, ranges=ranges)
    #readcsv(inputfile = "hist_230204_221718_ch.csv", datadir=datadir, outdir="C:\Data/TMU2023A/pdf/", outname="SDD_cal_Co57_0204.pdf",  cal_lines=cal_lines, ranges=ranges)

    #cal_lines=[5.899, 6.490] #keV
    #ranges=[[800,1050],[1000,1200]]
    #readcsv(inputfile = "hist_230204_235912.csv", datadir=datadir, outdir="C:\Data/TMU2023A/pdf/", outname="SDD_ene_Fe55_test.pdf",  cal_lines=cal_lines, ranges=ranges)


    ##for ene. hist
    #cal_lines=[22.162, 24.942]#Ag Kalpha1, Kbeta1
    #ranges=[[21.7,22.7],[24.4,25.4]]
    #readcsv(inputfile = "hist_230204_235912.csv", datadir=datadir, outdir="C:\Data/TMU2023A/pdf/", outname="SDD_ene_Ag_MeV.pdf",  cal_lines=cal_lines, ranges=ranges)

    #cal_lines=[6400, 211200., 211200., 23820.]#Pd Kalpha Kbeta
    #ranges=[[6000, 7500],[20600.,21600.],[20000.,25000.],[23300.,24300.]]
    ##readcsv(inputfile = "hist_230204_221718.csv", datadir=datadir, outdir="C:\Data/TMU2023A/pdf/", outname="SDD_ene_Co57.pdf",  cal_lines=cal_lines, ranges=ranges)
    #readcsv(inputfile = "hist_230205_220503.csv", datadir=datadir, outdir="C:\Data/TMU2023A/pdf/", outname="SDD_ene_MLFbeam27MeVc.pdf",  cal_lines=cal_lines, ranges=ranges)



    #datadir="C:\Data/RKN2024A/histograms/"
    #cal_lines=[39.522, 40.118, 45.294, 45.414] #keV
    ##ranges=[[39000, 40000],[39600,40600],[44800,45800],[35000,50000]]#[44900,45900]
    ##readcsv(inputfile = "20240203_Eu152_0_man.csv", datadir=datadir, outdir="C:\Data/RKN2024A/pdf/", outname="SDD_Eu152_0_man.pdf",  cal_lines=cal_lines, ranges=ranges)
    #ranges=[[3375, 3455],[3475,3550],[3850,4000],[4000,4100]]#[44900,45900]
    #readcsv(inputfile = "hist_240203_191350.csv", datadir=datadir, outdir="C:\Data/RKN2024A/pdf/", outname="SDD_Eu152_0_ch_lowgain.pdf",  cal_lines=cal_lines, ranges=ranges)
    

    datadir="C:\Data/MLF2024B/histograms/"
    cal_lines=[11.12, 20.48] #keV
    ranges=[[(11.12-0.5)/0.00725,(11.12+0.5)/0.00725 ],[(20.48-0.5)/0.00725 , (20.48+0.5)/0.00725 ]]#[44900,45900]
    readcsv(inputfile = "hist_momscan_18MeVc_250421_113136.csv", datadir=datadir, outdir="C:\Data/MLF2024B/pdf/", outname="Mom18MeV_Ar01atm.pdf",  cal_lines=cal_lines, ranges=ranges)

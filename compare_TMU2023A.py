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

PHAmsk=0x3fff
isGate=0x8000
CHmsk=0x000f

#def fit_gauss(data=[100, 101], f_range=[0, 100], d_range=[0,100]):
#  param = norm.fit(data)
#  return param
def get_energy(adc):
  #2023/01/31
  #a=5.7241e-03
  #b=2.5681e-2
  #a=5.7333e-03
  #b=5.6402e-2
  #2023/02/04
  a=5.748e-03
  b=1.077e-03
  return a*adc+b

def get_cal_lines(source="Fe55"):
  if source=="Fe55":
    ene_list=[5.899, 6.490]#Mn Kalpha, Kbeta
  elif source=="Co57":
    ene_list=[14.41300,14.41300+6.404]#epsilon, +Fe Kalpha
  elif source=="Eu152":
    ene_list=[39.522, 40.118, 45.294, 45.414, 46.578, 46.705]#Sm Kalpha2, Kalpha1, Kbeta3, Kbeta1, Kbeta2, Kbeta4
  elif source=="Pd":
    ene_list=[21.12, 23.82]#Pd Kalpha Kbeta
  elif source=="Ag":
    ene_list=[22.162, 21.990, 24.942]#Ag Kalpha1, Kalpha2, Kbeta1
    #ene_list=[22.162, 21.990, 24.942, 2.984, 2.978, 3.150, 3.347, 3.51959]#Ag Kalpha1, Kalpha2, Kbeta1, Lalpha1, Lalpha2, Lbeta1, Lgamma
  elif source=="AgPd":
    ene_list=[22.162, 21.990, 24.942, 21.12, 23.82]#Ag Kalpha1, Kalpha2, Kbeta1, Pd Kalpha Kbeta
  else :
    ene_list=[]
  return ene_list

def readfile(runs = [1], ENum = -1, datadir="C:\Data/TMU2022C/list/", outdir="C:\Data/TMU2022C/pdf/", outname="SDD_ene.pdf", cal_source=""):
  #bdata = open('C:\Data\list_000008_01-00000000001.bin','rb')
  counter=0
  trigIDs=list()
  phas=list()
  phas_wG=list()
  phas_woG=list()
  eners=list()
  eners_wG=list()
  eners_woG=list()
  gates=list()
  pufs=list()
  times=list()
  times_wG=list()
  times_woG=list()

  eners_MnKalpha=list()
  eners_MnKbeta =list()
  eners_Co57=list()
  eners_Eu152=list()

  colors=['b', 'k', 'm', 'g', 'y']
  
  cal_lines=get_cal_lines(cal_source)
  print("cal lines:",cal_lines)
  #Co57_ranges=[13.413+0.3,15.413-0.3]
  Co57_ranges=[14.2,14.65]
  fig,ax=plt.subplots(3,1)
  fig1,ax1=plt.subplots(3,1)
  fig2,ax2=plt.subplots(3,1)
  fig5, ax5=plt.subplots(2,1)
  for i, run in enumerate(runs):
    inputfile="list_{0:06d}_01-00000000001.bin".format(run)
    print('inputfile: ', datadir+inputfile)
    bdata = open(datadir+inputfile,'rb')
    while True:
      data = bdata.read(2)
      #print(data)
      #print(hex(data))
      idata=int.from_bytes(data, byteorder='big', signed=False)
      #print(hex(idata))

      if hex(idata) == '0x7fff':
        #print("header!!")
        data = bdata.read(2)#header
        #print("header2!!",data)
        data = bdata.read(4)#trig ID
        trigid=int.from_bytes(data, byteorder='big', signed=False)
        #print("trigID:", trigid)

        data = bdata.read(2)#U/G, U/V, Realtime1
        idata=int.from_bytes(data, byteorder='big', signed=False)
        gate=((idata>>15)&1) #gate flag
        veto=((idata>>14)&1) #veto flag
        time1=(idata&0x1fff) #not sure (will not be used anyway...)
        #print('time1:',time1)
        if gate > 0:
          #print('Gated!:',gate)
          gates.append(True)
        else:
          gates.append(False)
        #if veto > 0:
        #  print('Veto!:',veto)

        data = bdata.read(4)#Realtime2
        idata=int.from_bytes(data, byteorder='big', signed=False)
        time2=idata*0.01 #0.01: 10ns->us
        #print('time2:',time2)
        times.append(time2)
        #time2=((idata&0xfff)<<32) #not sure (will not be used anyway...)

        data = bdata.read(2)#Realtime3, unit, ch
        idata=int.from_bytes(data, byteorder='big', signed=False)
        ch=(idata&CHmsk)
        unit=((idata>>4)&CHmsk)

        data = bdata.read(2)#Pileup, pulse height
        idata=int.from_bytes(data, byteorder='big', signed=False)
        #print("pha:", hex(idata))
        pha=(idata&PHAmsk)
        ener=get_energy(pha) #calib. ch-> keV
        #ener=(pha-2500.)*(14.400-6.43)/(2500.-1110.)+14.400 #rough calib. ch-> keV  2022/03
        #ener=pha*14.400/2520. #rough calib. ch-> keV
        eners.append(ener)
        puf=((idata>>15)&1) #pileup flag
        #print("pha:", pha)
        #print("hex(pha):", hex(pha))
        #if(puf>0):
        #  print("pileup!:",puf)
        phas.append(pha)
        if(gate>0):
          phas_wG.append(pha)
          eners_wG.append(ener)
          times_wG.append(time2)
          if(ener>5.75) & (ener<6.25):
            eners_MnKalpha.append(ener)
          elif(ener>6.25) & (ener<6.75):
            eners_MnKbeta.append(ener)
          elif(ener>Co57_ranges[0]) & (ener<Co57_ranges[1]):
            eners_Co57.append(ener)
          elif(ener>38) & (ener<42):
            eners_Eu152.append(ener)
        else:
          phas_woG.append(pha)
          eners_woG.append(ener)
          times_woG.append(time2)

      if len(data) == 0:
        break
      if counter >= ENum & ENum > 0:
        break
      counter += 1
    bdata.close

                   # 3-2, 4-3,  5-4
    muonicO_xrays=[24.86, 8.7, 4.03]
                   #  5-4,   6-5,   7-6,  8-7,  9-8, 10-9,11-10
    muonicFe_xrays=[42.78, 23.23, 14.01, 9.09, 6.23, 4.46, 3.3]
                   #  5-4,   6-5,   7-6,   8-7,  9-8, 10-9,11-10
    muonicCu_xrays=[53.25, 28.91, 17.43, 11.31, 7.75, 5.55, 4.10]
                   #  4-3,   5-4,   6-5,  7-6,  8-7,  9-8
    muonicAr_xrays=[44.26, 20.48, 11.12, 6.71, 4.35, 2.98]
                     #   5-3,  6-4,   6-3,
    muonicO_xrays_sub=[12.73, 6.21, 14.91]
                     #    7-5,   8-6,   8-5,   9-7,   9-6,   9-5,
    muonicFe_xrays_sub=[37.23, 23.09, 46.32, 15.32, 29.32, 52.55 ]

                    #  2-1    3-2
    muonicBe_xrays=[33.378, 6.18]
                    #  4-3   5-4   6-5 
    muonicAl_xrays=[23.04, 10.66, 5.79 ]
                    # 7-6    8-7    9-8   10-9  11-10 12-11  13-12 14-13
    muonicAg_xrays=[45.86, 29.75, 20.39, 14.58, 10.79, 8.21, 6.39, 5.07]

    ebins1=np.arange( 2,30,0.1)
    ebins2=np.arange(30,48,0.1)
    ax[0].hist(phas    , bins=1000, range=[500,10000],histtype='step',color=color[i], label='run{0:04d}'.format(run)     )
    #ax[1].hist(phas_wG , bins=1700, range=[1000,2700],histtype='step',color='r', label='w/ Gate' )
    ax[1].hist(phas_wG , bins=1000, range=[500,10000],histtype='step',color=color[i])
    ax[2].hist(phas_woG, bins=1000, range=[500,10000],histtype='step',color=color[i])
    ax[2].set_xlabel('Pulse height [ch]')
    #ax[0].set_xlim(500,10000)
    #ax[1].set_xlim(500,10000)
    #ax[2].set_xlim(500,10000)
    #ax[1].axvline(1110,color='b',linewidth=1)
    #ax[1].axvline(2500,color='b',linewidth=1)
    ax[1].set_xticks(np.arange(1000,2700,100), minor=True)
    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    #plt.show()

    #ax1[0].hist(phas    , bins=200, range=[2250,2750],histtype='step',color='k', label='all'     )
    ax1[0].hist(phas    , bins=1000, range=[0,10000],histtype='step',color=color[i], label='run{0:04d}'.format(run)     )
    ax1[1].hist(phas_wG , bins=1000, range=[0,10000],histtype='step',color=color[i] )
    ax1[2].hist(phas_woG, bins=1000, range=[0,10000],histtype='step',color=color[i] )
    ax1[2].set_xlabel('Pulse height [ch]')
    ax1[0].legend()
    ax1[1].legend()
    ax1[2].legend()
    #plt.show()

    bins_full_range=np.arange(0,48,0.1)
    #Fe_xrays=[6.46
    ax2[0].hist(eners    , bins=bins_full_range,histtype='step',color=color[i], label='all'     )
    #ax2[0].hist(eners    , bins=550, range=[4,15],histtype='step',color='k', label='all'     )
    ax2[1].hist(eners_wG , bins=bins_full_range,histtype='step',color=color[i], label='w/ Gate' )
    ax2[2].hist(eners_wG, bins=100, range=[20.,30.],histtype='step',color=color[i], label='w/ Gate (zoom)')
    #ax2[2].hist(eners_wG, bins=100, range=[4.9,6.9],histtype='step',color='r', label='w/ Gate (zoom)')
    #ax2[2].hist(eners_wG, bins=100, range=[2,4],histtype='step',color='r', label='w/ Gate (zoom)')
    #ax2[2].hist(eners_woG, bins=580, range=[2,60],histtype='step',color='b', label='w/o Gate')

    ax2[0].set_yscale('log')
    ax2[1].set_yscale('log')
    ax2[1].grid(which='minor', axis='x', linestyle='dashed')
    ax2[1].grid(which='major', axis='x')
    ax2[1].set_xticks(np.arange(0,48,1), minor=True)
    ax2[0].grid(which='minor', axis='x', linestyle='dashed')
    ax2[0].grid(which='major', axis='x')
    ax2[0].set_xticks(np.arange(0,48,1), minor=True)
    #ax2[0].axvline(6.40,color='b',linewidth=1)
    #ax2[0].axvline(7.06,color='b',linewidth=1)
    #ax2[0].axvline(14.4,color='b',linewidth=1)

    #for muO in muonicO_xrays:
    #  ax2[1].axvline(muO,color='b',linewidth=1)
    #for muOsub in muonicO_xrays_sub:
    #  ax2[1].axvline(muOsub,color='b',linewidth=1,linestyle='dashed')
    #for muFe in muonicFe_xrays:
    #  ax2[1].axvline(muFe,color='m',linewidth=1)
    #for muFesub in muonicFe_xrays_sub:
    #  ax2[1].axvline(muFesub,color='m',linewidth=1, linestyle='dashed')
    #for muAr in muonicAr_xrays:
    #  ax2[1].axvline(muAr,color='m',linewidth=1)


    #for muCu in muonicCu_xrays:
    #  ax2[1].axvline(muCu,color='k',linewidth=1)
    ax2[2].set_xlabel('energy [keV] (do not trust me!)')
    ax2[0].legend()
    ax2[1].legend()
    ax2[2].legend()

    ax5[0].hist(eners_wG    , bins=ebins1,histtype='step',color=color[i], label='run{0:04d}'.format(run)     )
    ax5[1].hist(eners_wG    , bins=ebins2,histtype='step',color=color[i], label='run{0:04d}'.format(run)     )
    ax5[0].set_xlim(ebins1[0],ebins1[-1])
    ax5[1].set_xlim(ebins2[0],ebins2[-1])
    ax5[0].set_xticks(np.arange( 2,30,1), minor=True)
    ax5[1].set_xticks(np.arange(30,48,1), minor=True)
    ax5[1].set_xlabel('Energy [keV]')
    ax5[0].set_yscale('log')
    ax5[1].set_yscale('log')
    ax5[0].grid(which='minor', axis='x', linestyle='dashed')
    ax5[0].grid(which='major', axis='x')
    ax5[1].grid(which='minor', axis='x', linestyle='dashed')
    ax5[1].grid(which='major', axis='x')

  for cal_line in cal_lines:
    ax2[1].axvline(cal_line,color='b',linewidth=1, linestyle='dashed')
    ax2[2].axvline(cal_line,color='b',linewidth=1, linestyle='dashed')
    ax5[0].axvline(cal_line,color='b',linewidth=1, linestyle='dashed')
    ax5[1].axvline(cal_line,color='b',linewidth=1, linestyle='dashed')
  #for muX in muonicAg_xrays:
  #  if (muX < 30) & (muX > 5):
  #    ax5[0].axvline(muX,color='g',linewidth=1)
  #  elif (muX > 30) & (muX < 55):
  #    ax5[1].axvline(muX,color='g',linewidth=1)
  #for muX in muonicAl_xrays:
  #  if (muX < 30) & (muX > 5):
  #    ax5[0].axvline(muX,color='m',linewidth=1)
  #  elif (muX > 30) & (muX < 55):
  #    ax5[1].axvline(muX,color='m',linewidth=1)
  #for muX in muonicBe_xrays:
  #  if (muX < 30) & (muX > 5):
  #    ax5[0].axvline(muX,color='y',linewidth=1)
  #  elif (muX > 30) & (muX < 55):
  #    ax5[1].axvline(muX,color='y',linewidth=1)
  pp=PdfPages(outdir+outname)
  fig.tight_layout()
  pp.savefig(fig)
  fig1.tight_layout()
  #fig.savefig(outdir+outname)
  pp.savefig(fig1)
  fig2.tight_layout()
  #fig.savefig(outdir+outname)
  pp.savefig(fig2)
  fig5.tight_layout()
  pp.savefig(fig5)

  pp.close()
  print(outdir+outname,'is saved.')

if __name__ == "__main__":


    runlists=[]
    str_runlist=''
    runs= [13,15,16,17]
    for run in runs :
      str_runlist+="{0:04d}_".format(run)
    readfile(runs = runlists, ENum = -1, datadir="C:\Data/TMU2023A/list/", outdir="C:\Data/TMU2023A/pdf/", outname="SDD_ene_MLF_run{0:s}compare.pdf".format(str_runlist),  cal_source="AgPd")

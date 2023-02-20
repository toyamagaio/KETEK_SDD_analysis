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

def readfile(inputfile = "list_000008_01-00000000001.bin", ENum = -1, datadir="C:\Data/TMU2022C/list/", outdir="C:\Data/TMU2022C/pdf/", outname="SDD_ene.pdf"):
  #bdata = open('C:\Data\list_000008_01-00000000001.bin','rb')
  print('inputfile: ', datadir+inputfile)
  bdata = open(datadir+inputfile,'rb')
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
      ener=(pha-2500.)*(14.400-6.43)/(2500.-1110.)+14.400 #rough calib. ch-> keV
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

  pp=PdfPages(outdir+outname)
  fig,ax=plt.subplots(3,1)
  ax[0].hist(phas    , bins=1000, range=[500,10000],histtype='step',color='k', label='all'     )
  #ax[1].hist(phas_wG , bins=1700, range=[1000,2700],histtype='step',color='r', label='w/ Gate' )
  ax[1].hist(phas_wG , bins=1000, range=[500,10000],histtype='step',color='r', label='w/ Gate' )
  ax[2].hist(phas_woG, bins=1000, range=[500,10000],histtype='step',color='b', label='w/o Gate')
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
  fig.tight_layout()
  pp.savefig(fig)

  fig1,ax1=plt.subplots(3,1)
  #ax1[0].hist(phas    , bins=200, range=[2250,2750],histtype='step',color='k', label='all'     )
  ax1[0].hist(phas    , bins=1000, range=[0,10000],histtype='step',color='k', label='all'     )
  ax1[1].hist(phas_wG , bins=1000, range=[0,10000],histtype='step',color='r', label='w/ Gate' )
  ax1[2].hist(phas_woG, bins=1000, range=[0,10000],histtype='step',color='b', label='w/o Gate')
  ax1[2].set_xlabel('Pulse height [ch]')
  ax1[0].legend()
  ax1[1].legend()
  ax1[2].legend()
  #plt.show()
  fig1.tight_layout()
  #fig.savefig(outdir+outname)
  pp.savefig(fig1)

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
  muonicFe_xrays_sub=[37.23, 23.09, 46.32, 15.32, 29.32, 52.55, ]
  #Fe_xrays=[6.6
  fig2,ax2=plt.subplots(3,1)
  ax2[0].hist(eners    , bins=580, range=[2,60],histtype='step',color='k', label='all'     )
  #ax2[0].hist(eners    , bins=550, range=[4,15],histtype='step',color='k', label='all'     )
  ax2[1].hist(eners_wG , bins=580, range=[2,60],histtype='step',color='r', label='w/ Gate' )
  ax2[2].hist(eners_wG, bins=100, range=[4.9,6.9],histtype='step',color='r', label='w/ Gate (zoom)')
  #ax2[2].hist(eners_wG, bins=100, range=[2,4],histtype='step',color='r', label='w/ Gate (zoom)')
  #ax2[2].hist(eners_woG, bins=580, range=[2,60],histtype='step',color='b', label='w/o Gate')

  ax2[1].set_yscale('log')
  ax2[1].grid(which='minor', axis='x', linestyle='dashed')
  ax2[1].grid(which='major', axis='x')
  ax2[1].set_xticks(np.arange(0,60,1), minor=True)
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

  for muAr in muonicAr_xrays:
    ax2[1].axvline(muAr,color='m',linewidth=1)

  #for muCu in muonicCu_xrays:
  #  ax2[1].axvline(muCu,color='k',linewidth=1)
  ax2[2].set_xlabel('energy [keV] (do not trust me!)')
  ax2[0].legend()
  ax2[1].legend()
  ax2[2].legend()
  #plt.show()
  fig2.tight_layout()
  #fig.savefig(outdir+outname)
  pp.savefig(fig2)

  fig3,ax3=plt.subplots(3,1)
  ax3[0].hist(times    , bins=1000, range=[0,100],histtype='step',color='k', label='all'     )
  ax3[1].hist(times_wG , bins=1000, range=[0,100],histtype='step',color='r', label='w/ Gate' )
  ax3[2].hist(times_woG, bins=1000, range=[0,100],histtype='step',color='b', label='w/o Gate')
  ax3[2].set_xlabel('time [us]')
  ax3[0].legend()
  ax3[1].legend()
  ax3[2].legend()
  #plt.show()
  fig3.tight_layout()
  #fig.savefig(outdir+outname)
  pp.savefig(fig3)


  #fit
  ranges=[5.75,6.75]
  bins=50
  wrange=ranges[1]-ranges[0]
  wbin=wrange/bins
  integ_MnKalpha=len(eners_MnKalpha)
  integ_MnKbeta =len(eners_MnKbeta )
  print('wbin:',wbin,', integral Ka:',integ_MnKalpha,', integral Kb:',integ_MnKbeta)

  param_MnKalpha =norm.fit(eners_MnKalpha)
  param_MnKbeta  =norm.fit(eners_MnKbeta )
  print('param_MnKalpha:',param_MnKalpha)
  print('param_MnKbeta :',param_MnKbeta )
  x = np.linspace(5.75,6.75,100)
  #FWHM_MnKalpha=2*math.sqrt(2*math.log(2)*param_MnKalpha[1])
  #FWHM_MnKbeta =2*math.sqrt(2*math.log(2)*param_MnKbeta [1])
  FWHM_MnKalpha=2.35*param_MnKalpha[1]
  FWHM_MnKbeta =2.35*param_MnKbeta [1]
  print('FWHM Ka:',FWHM_MnKalpha,', Kb:',FWHM_MnKbeta)
  pdf_fitted_MnKalpha = norm.pdf(x,loc=param_MnKalpha[0], scale=param_MnKalpha[1])*integ_MnKalpha*wbin
  pdf_fitted_MnKbeta  = norm.pdf(x,loc=param_MnKbeta [0], scale=param_MnKbeta [1])*integ_MnKbeta *wbin
  fig4, ax4=plt.subplots(2,1)
  #plt.title('Mn Kalpha')
  ax4[0].plot(x, pdf_fitted_MnKalpha, 'r-')#, x,pdf, 'b-')
  ax4[0].plot(x, pdf_fitted_MnKbeta , 'r-')#, x,pdf, 'b-')
  ax4[0].hist(np.hstack([eners_MnKalpha,eners_MnKbeta]), bins=bins, range=ranges,histtype='step',color='r', label='Mn Ka, Kb', alpha=.3)
  ax4[0].legend()
  ax4[1].hist(np.hstack([eners_MnKalpha,eners_MnKbeta]), bins=bins, range=ranges,histtype='step',color='r', label='Mn Ka, Kb', alpha=.3)
  ax4[1].set_xlabel('energy [keV]')
  fig4.tight_layout()
  pp.savefig(fig4)

  pp.close()
  print(outdir+outname,'is saved.')

if __name__ == "__main__":
    #if len(sys.argv)==2:
    #  #readfile(sys.argv[1])
    #  runnum=int(sys.argv[1])
    #  input_binary=f'list_{runnum:06}_01-00000000001.bin'
    #  output_pdf  =f'SDD_PH_run{runnum:06}.pdf'
    #  print("input:",input_binary)
    #  print("output:",output_pdf)
    #  readfile(inputfile = input_binary, ENum = -1, datadir="C:\Data/TMU2022C/list/", outdir="C:\Data/TMU2022C/pdf/", outname=output_pdf)#test 3/16
    #else:
    #  readfile(inputfile = "list_000008_01-00000000001.bin", ENum = -1, datadir="C:\Data/", outdir="C:\Data/TMU2022C/pdf/", outname="SDD_ene_Co57.pdf")#test 3/16

    readfile(inputfile = "list_000001_01-00000000001.bin", ENum = -1, datadir="C:\Data/TMU2023A/list/", outdir="C:\Data/TMU2023A/pdf/", outname="SDD_ene_Fe55.pdf")

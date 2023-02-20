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
import h5py

PHAmsk=0x3fff
isGate=0x8000
CHmsk=0x000f

#def fit_gauss(data=[100, 101], f_range=[0, 100], d_range=[0,100]):
#  param = norm.fit(data)
#  return param
def ngaussian(x, area, mu, sigma):
   return area/(math.sqrt(2.*math.pi*sigma*sigma))*np.exp(-np.power((x - mu), 2.)/(2.*sigma*sigma))
def gaussian(x, peak, mu, sigma):
   return peak*np.exp(-np.power((x - mu), 2.)/(2.*sigma*sigma))

def getNearestValue(list, num):
  idx = np.abs(np.asarray(list) - num).argmin()
  return list[idx]


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

def readfile(inputfiles = ["list_000008_01-00000000001.bin"], ENum = -1, datadir="C:\Data/TMU2022C/list/", outdir="C:\Data/TMU2022C/pdf/", outname="SDD_ene.pdf", cal_source=""):
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
  
  cal_lines=get_cal_lines(cal_source)
  print("cal lines:",cal_lines)
  #Co57_ranges=[13.413+0.3,15.413-0.3]
  Co57_ranges=[14.2,14.65]
  for inputfile in inputfiles:
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

  muonicHe3_xray=8.13
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

                   # 3-2      4-2
  muonicC_xrays  =[13.95] # 18.83

                  #  2-1   3-2    4-2     3-1     4-1
  muonicBe_xrays=[33.378, 6.18, 8.343, 39.558, 41.721]
                  #  4-3   5-4   6-3     6-4      6-5    7-4
  muonicAl_xrays=[23.04, 10.66, 39.50, 16.46] #  5.79, 19.95
                  # 7-6    8-7    9-8   10-9 #  11-10 12-11  13-12 14-13
  muonicAg_xrays=[45.86, 29.75, 20.39, 14.58]#, 10.79, 8.21,  6.39, 5.07]
                  # 4-3    5-4    6-4  # 6-5
  muonicSi_xrays=[26.73, 12.37, 19.09] #6.72,
  bins_full_range=np.arange(0,48,0.1)
  #Fe_xrays=[6.46
  fig2,ax2=plt.subplots(3,1)
  ax2[0].hist(eners    , bins=bins_full_range,histtype='step',color='k', label='all'     )
  #ax2[0].hist(eners    , bins=550, range=[4,15],histtype='step',color='k', label='all'     )
  ax2[1].hist(eners_wG , bins=bins_full_range,histtype='step',color='r', label='w/ Gate' )
  ax2[2].hist(eners_wG, bins=50, range=[7.,9.],histtype='step',color='r', label='w/ Gate (zoom)')
  #ax2[2].hist(eners_wG, bins=100, range=[20.,30.],histtype='step',color='r', label='w/ Gate (zoom)') # for Ag target
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

  

  for cal_line in cal_lines:
    ax2[1].axvline(cal_line,color='b',linewidth=1, linestyle='dashed')
    #ax2[2].axvline(cal_line,color='b',linewidth=1, linestyle='dashed')
  ax2[2].axvline(muonicHe3_xray,color='b',linewidth=1, linestyle='-.')

  #for muCu in muonicCu_xrays:
  #  ax2[1].axvline(muCu,color='k',linewidth=1)
  ax2[2].set_xlabel('energy [keV]')
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


  ##fit
  fig4, ax4=plt.subplots(2,1)
  #Fe55
  if cal_source == 'Fe55':
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
    FWHM_MnKalpha=2*param_MnKalpha[1]*math.sqrt(2*math.log(2))
    FWHM_MnKbeta =2*param_MnKbeta [1]*math.sqrt(2*math.log(2))
    #FWHM_MnKalpha=2.35*param_MnKalpha[1]
    #FWHM_MnKbeta =2.35*param_MnKbeta [1]
    print('FWHM Ka:',FWHM_MnKalpha,', Kb:',FWHM_MnKbeta)
    pdf_fitted_MnKalpha = norm.pdf(x,loc=param_MnKalpha[0], scale=param_MnKalpha[1])*integ_MnKalpha*wbin
    pdf_fitted_MnKbeta  = norm.pdf(x,loc=param_MnKbeta [0], scale=param_MnKbeta [1])*integ_MnKbeta *wbin
    #plt.title('Mn Kalpha')
    ax4[0].plot(x, pdf_fitted_MnKalpha, 'r-')#, x,pdf, 'b-')
    ax4[0].plot(x, pdf_fitted_MnKbeta , 'r-')#, x,pdf, 'b-')
    ax4[0].hist(np.hstack([eners_MnKalpha,eners_MnKbeta]), bins=bins, range=ranges,histtype='step',color='r', label='Mn Ka, Kb', alpha=.3)
    ax4[0].legend()
    ax4[1].hist(np.hstack([eners_MnKalpha,eners_MnKbeta]), bins=bins, range=ranges,histtype='step',color='r', label='Mn Ka, Kb', alpha=.3)
    ax4[1].set_xlabel('energy [keV]')
    ####
    #Co57 (14.41300)
  elif cal_source == 'Co57':
    bins=20
    wrange=Co57_ranges[1]-Co57_ranges[0]
    wbin=wrange/bins
    integ_Co57=len(eners_Co57)
    print('wbin:',wbin,', integral Co57 peak:',integ_Co57)
    param_Co57 =norm.fit(eners_Co57)
    print('param_Co57:',param_Co57)
    x = np.linspace(Co57_ranges[0],Co57_ranges[1],100)
    FWHM_Co57=2*param_Co57[1]*math.sqrt(2*math.log(2))
    #FWHM_Co57=2.35*param_Co57[1]
    print('FWHM Co57:',FWHM_Co57)
    pdf_fitted_Co57 = norm.pdf(x,loc=param_Co57[0], scale=param_Co57[1])*integ_Co57*wbin
    fig4, ax4=plt.subplots(2,1)
    ax4[0].plot(x, pdf_fitted_Co57, 'r-')#, x,pdf, 'b-')
    ax4[0].hist(eners_Co57, bins=bins, range=Co57_ranges,histtype='step',color='r', label='Co57', alpha=.3)
    ax4[0].legend()
    ax4[1].hist(eners_Co57, bins=bins, range=Co57_ranges,histtype='step',color='r', label='Co57', alpha=.3)
    ax4[1].set_xlabel('energy [keV]')

    fig4.tight_layout()
    pp.savefig(fig4)

  fig5, ax5=plt.subplots(2,1)
  ebins1=np.arange( 2,30,0.1)
  ebins2=np.arange(30,48,0.1)
  #ax5[0].hist(eners    , bins=ebins1,histtype='step',color='k', label='all'     )
  #ax5[1].hist(eners    , bins=ebins2,histtype='step',color='k', label='all'     )
  ax5[0].hist(eners_wG    , bins=ebins1,histtype='step',color='r', label='w/ Gate'     )
  ax5[1].hist(eners_wG    , bins=ebins2,histtype='step',color='r', label='w/ Gate'     )
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
    ax5[0].axvline(cal_line,color='b',linewidth=1, linestyle='dashed')
    ax5[1].axvline(cal_line,color='b',linewidth=1, linestyle='dashed')
  for muX in muonicAg_xrays:
    if (muX < 30) & (muX > 5):
      ax5[0].axvline(muX,color='g',linewidth=1)
    elif (muX > 30) & (muX < 55):
      ax5[1].axvline(muX,color='g',linewidth=1)
  for muX in muonicAl_xrays:
    if (muX < 30) & (muX > 5):
      ax5[0].axvline(muX,color='m',linewidth=1)
    elif (muX > 30) & (muX < 55):
      ax5[1].axvline(muX,color='m',linewidth=1)
  for muX in muonicBe_xrays:
    if (muX < 30) & (muX > 5):
      ax5[0].axvline(muX,color='y',linewidth=1)
    elif (muX > 30) & (muX < 55):
      ax5[1].axvline(muX,color='y',linewidth=1)
  for muX in muonicSi_xrays:
    if (muX < 30) & (muX > 5):
      ax5[0].axvline(muX,color='c',linewidth=1)
    elif (muX > 30) & (muX < 55):
      ax5[1].axvline(muX,color='c',linewidth=1)
  for muX in muonicC_xrays:
    if (muX < 30) & (muX > 5):
      ax5[0].axvline(muX,color='#ff7f00',linewidth=1)
    elif (muX > 30) & (muX < 55):
      ax5[1].axvline(muX,color='#ff7f00',linewidth=1)
  pp.savefig(fig5)

  ####
  #muHe3 bg fit
  fig6, ax6=plt.subplots(3,1, sharex=True)
  muHe_bins=np.arange(7,9,4.*5.748e-03)
  np_eners_wG=np.array(eners_wG)
  muHe_bgcut = ((np_eners_wG>7.0) & (np_eners_wG<7.75)) | ((np_eners_wG>8.5) & (np_eners_wG<9.00))
  print('muHe_bgcut:',muHe_bgcut, ' type(muHe_bgcut)', type(muHe_bgcut), ' len(muHe_bgcut)', len(muHe_bgcut))
  ncount_bg, bins = np.histogram(np_eners_wG[muHe_bgcut],muHe_bins)
  ncount   , bins = np.histogram(np_eners_wG            ,muHe_bins)
  delta_E = (bins[:-1]+bins[1:])*0.5
  nonzero_nc_idx = tuple(list(np.where(ncount_bg>0)))
  #nonzero_nc_idx = list(np.where(ncount_bg>0))  #warning appeared
  coeff, cov = np.polyfit(delta_E[nonzero_nc_idx], ncount_bg[nonzero_nc_idx],1, w=1/np.sqrt(ncount_bg[nonzero_nc_idx]),cov=True)
  print('coeff:',coeff)
  print('delta_E:',delta_E)
  print('bins:',bins)

  ax6[0].errorbar(delta_E, ncount, yerr=np.sqrt(ncount) , drawstyle='steps-mid', color='r')
  ax6[0].plot(delta_E, coeff[0]*delta_E+coeff[1], color='m', linestyle='dashed')
  ax6[1].errorbar(delta_E, ncount-(coeff[0]*delta_E+coeff[1]), yerr=np.sqrt(ncount) , drawstyle='steps-mid', color='r')
  muBe_8keV=getNearestValue(muonicBe_xrays,8.)
  muBe_shift= -0.02
  #print(muBe_8keV)
  x=np.linspace(7,9,100)
  peak=max(ncount-(coeff[0]*delta_E+coeff[1]))
  print('peak: ',peak)
  print('mean counts (after sub. pol1bg, gauss): ',np.mean(ncount-(coeff[0]*delta_E+coeff[1])-gaussian(delta_E,peak,mu=muBe_8keV+muBe_shift,sigma=0.06)))
  ax6[1].plot(x,gaussian(x,peak,mu=muBe_8keV+muBe_shift,sigma=0.06),'r-',alpha=.3) #sigma 0.081@14.4keV, 0.060@5.9keV,  0.062@6.4keV
  #print(gaussian(delta_E,peak,mu=muBe_8keV,sigma=0.06))
  ax6[1].axvline(muBe_8keV,color='y',linewidth=1)
  ax6[2].errorbar(delta_E, ncount-(coeff[0]*delta_E+coeff[1])-gaussian(delta_E,peak,mu=muBe_8keV+muBe_shift,sigma=0.06), yerr=np.sqrt(ncount) , drawstyle='steps-mid', color='r')
  ax6[2].axhline(0,color='k',linewidth=1, linestyle='dashed')
  ax6[2].set_xlabel('Energy [keV]')
  ax6[2].axvline(muonicHe3_xray,color='b',linewidth=1, linestyle='-.')
  pp.savefig(fig6)
  

  pp.close()
  print(outdir+outname,'is saved.')

  ##hdf5 output (test)
  hdf_cut= (np_eners_wG) > 7 & (np_eners_wG < 10)
  outname_hdf5 =outname.replace('.pdf','.hdf5')
  print(outname_hdf5)
  with h5py.File(outdir+outname_hdf5, mode='w') as f:
    dataset = f.create_dataset(name='/SDD/eners_wG', data=np_eners_wG)

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

    #readfile(inputfiles = ["list_000001_01-00000000001.bin"], ENum = -1, datadir="C:\Data/TMU2023A/list/", outdir="C:\Data/TMU2023A/pdf/", outname="SDD_ene_Fe55.pdf", cal_source="Fe55")
    #readfile(inputfile = "list_000005_01-00000000001.bin", ENum = -1, datadir="C:\Data/TMU2023A/list/", outdir="C:\Data/TMU2023A/pdf/", outname="SDD_ene_Co57.pdf",  cal_source="Co57")
    #readfile(inputfile = "list_000006_01-00000000001.bin", ENum = -1, datadir="C:\Data/TMU2023A/list/", outdir="C:\Data/TMU2023A/pdf/", outname="SDD_ene_Eu152.pdf",  cal_source="Eu152")
    #readfile(inputfile = "list_000005_01-00000000001.bin", ENum = -1, datadir="C:\Data/TMU2023A/list/", outdir="C:\Data/TMU2023A/pdf/", outname="SDD_ene_Co57_testPd.pdf",  cal_source="Pd")
    #readfile(inputfile = "list_000006_01-00000000001.bin", ENum = -1, datadir="C:\Data/TMU2023A/list/", outdir="C:\Data/TMU2023A/pdf/", outname="SDD_ene_Eu152_testPd.pdf",  cal_source="Pd")
    #readfile(inputfiles = ["list_000004_01-00000000001.bin"], ENum = -1, datadir="C:\Data/TMU2023A/list/", outdir="C:\Data/TMU2023A/pdf/", outname="SDD_ene_Co57_20230204.pdf",  cal_source="Co57")

    #readfile(inputfiles = ["list_000008_01-00000000001.bin"], ENum = -1, datadir="C:\Data/TMU2023A/list/", outdir="C:\Data/TMU2023A/pdf/", outname="SDD_ene_MLF_run0008.pdf",  cal_source="Pd")
    #readfile(inputfiles = ["list_000009_01-00000000001.bin"], ENum = -1, datadir="C:\Data/TMU2023A/list/", outdir="C:\Data/TMU2023A/pdf/", outname="SDD_ene_MLF_run0009.pdf",  cal_source="Pd")

    #for run in range(41,48):
    #  readfile(inputfiles = ["list_{0:06d}_01-00000000001.bin".format(run)], ENum = -1, datadir="C:\Data/TMU2023A/list/", outdir="C:\Data/TMU2023A/pdf/", outname="SDD_ene_MLF_run{0:04d}.pdf".format(run),  cal_source="AgPd")

    run=118 #Ag 30MeV/c
    readfile(inputfiles = ["list_{0:06d}_01-00000000001.bin".format(run)], ENum = -1, datadir="C:\Data/TMU2023A/list/", outdir="C:\Data/TMU2023A/pdf/", outname="SDD_ene_MLF_run{0:04d}.pdf".format(run),  cal_source="AgPd")

    inputfiles=[]
    ##all D2 run
    #irun=21
    #frun=102

    #irun=10
    #irun=8
    #frun=20

    #Ag bg 27.5MeV/c
    irun=103
    frun=117
    for run in range(irun,frun+1):
      inputfiles.append("list_{0:06d}_01-00000000001.bin".format(run))
    #readfile(inputfiles =inputfiles, ENum = -1, datadir="C:\Data/TMU2023A/list/", outdir="C:\Data/TMU2023A/pdf/", outname="SDD_ene_MLF_run{0:04d}_{1:04d}.pdf".format(irun, frun),  cal_source="AgPd")

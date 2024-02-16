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
cmap=plt.get_cmap('tab10')

PHAmsk=0x3fff
isGate=0x8000
CHmsk=0x000f

alpha=0.3

#def fit_gauss(data=[100, 101], f_range=[0, 100], d_range=[0,100]):
#  param = norm.fit(data)
#  return param
def get_energy(adc):

  #2024/02/06 by muAr (rough)
  a=0.0115
  b=-0.067
  #2024/02/05
  #a=1.0190e-02
  #b=4.5443e+00
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
  elif source=="muAr":
    ene_list=[11.12, 20.48]#
  elif source=="muSi":
    #ene_list=[12.37, 26.73]#
    ene_list=[13.95, 24.9]# sonoba shinogi value
  elif source=="muLi":
    ene_list=[18.70, 22.167]#
  else :
    ene_list=[]
  return ene_list

def readfile(runs = [1], ENum = -1, datadir="C:\Data/TMU2022C/list/", outdir="C:\Data/TMU2022C/pdf/", outname="SDD_ene.pdf", cal_source="",labels=None, noise_E=9.0, rfactor_noi=1.):
  #bdata = open('C:\Data\list_000008_01-00000000001.bin','rb')
  counter=0

  colors=['k', 'r', 'b', 'g', 'y']

  if len(runs)>5:
    colors=[cmap(i) for i in range(len(runs))]
  
  cal_lines=get_cal_lines(cal_source)
  print("cal lines:",cal_lines)
  #Co57_ranges=[13.413+0.3,15.413-0.3]
  Co57_ranges=[14.2,14.65]
  fig,ax=plt.subplots(3,1)
  fig1,ax1=plt.subplots(3,1)
  fig2,ax2=plt.subplots(2,1)
  fig5, ax5=plt.subplots(2,1)
  delta_E=0.5

  if labels == None:
    labels=['{0:04d}'.format(run) for run in runs]

  for i, (run, label) in enumerate(zip(runs,labels)):
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

                   #
    muonicN_xrays=[19.0, 6.65, 25.7]
                   # 3-2, 4-3,  5-4
    muonicO_xrays=[24.86, 8.7, 4.03]
                   #  5-4,   6-5,   7-6,  8-7,  9-8, 10-9,11-10
    muonicFe_xrays=[42.78, 23.23, 14.01, 9.09, 6.23, 4.46, 3.3]
                   #  6-5,   7-6,   8-7,  9-8, 10-9,11-10    5-4, 
    muonicCu_xrays=[28.91, 17.43, 11.31, 7.75, 5.55, 4.10] #53.25, 
                   #  4-3,   5-4,   6-5,  7-6,  8-7,  9-8
    muonicAr_xrays=[44.26, 20.48, 11.12, 6.71, 4.35, 2.98]
                         #6-4,   7-5,   8-6,  9-7    9-6
    muonicAr_xrays_sub=[31.60, 17.83, 11.06, 7.34, 14.04]
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
                    #  3-2     4-3    5-4    6-5   7-6
    muonicSi_xrays =[ 76.42,  26.73, 12.37,  6.72, 4.05]
                       #   5-3    6-4    7-5   8-6    6-3    7-4    8-5   9-6
    muonicSi_xrays_sub=[ 39.10, 19.09, 10.77, 6.68, 45.82, 23.14, 13.40, 8.48]

    muonicTb_xrays =[ 56.98, 39.05, 27.92, 20.65, 15.71, 12.22, 9.7]
    muonicC_xrays  =[13.95]
                         # 5-2,   4-2
    muonicC_xrays_sub  =[21.09, 18.83]
                      # 6-5    7-6   8-7
    muonicCo_xrays  =[25.06, 15.11, 9.80]

    print('=====run{0:04d} {1}====='.format(run,label))
    #for muAr_xrays in muonicAr_xrays:
    #  cut=(eners_wG<muAr_xrays+0.3) & (eners_wG>muAr_xrays-0.3)
    noise_arr   =np.where( (np.array(eners_wG)>noise_E-rfactor_noi*delta_E) & (np.array(eners_wG)<noise_E+rfactor_noi*delta_E))
    noise_counts=len(noise_arr[0])/rfactor_noi
    print('noise counts: {}'.format(noise_counts))

    for cal_E in cal_lines:
        narr=np.where( (np.array(eners_wG)>cal_E-delta_E) & (np.array(eners_wG)<cal_E+delta_E))
        ncounts=len(narr[0])
        print('{}keV, nEv: {}, signal: {}, S/N: {}'.format(cal_E, ncounts, ncounts-noise_counts, ncounts/noise_counts))


    ebins1=np.arange( 2,30,0.1)
    ebins2=np.arange(30,48,0.1)
    ax[0].hist(phas    , bins=1000, range=[500,10000],histtype='step',color=colors[i], label='run{0:04d}'.format(run)     )
    #ax[1].hist(phas_wG , bins=1700, range=[1000,2700],histtype='step',color='r', label='w/ Gate' )
    ax[1].hist(phas_wG , bins=1000, range=[500,10000],histtype='step',color=colors[i])
    ax[2].hist(phas_woG, bins=1000, range=[500,10000],histtype='step',color=colors[i])
    ax[2].set_xlabel('Pulse height [ch]')
    #ax[0].set_xlim(500,10000)
    #ax[1].set_xlim(500,10000)
    #ax[2].set_xlim(500,10000)
    #ax[1].axvline(1110,color='b',linewidth=1)
    #ax[1].axvline(2500,color='b',linewidth=1)
    ax[1].set_xticks(np.arange(1000,2700,100), minor=True)
    ax[0].legend()
    #ax[1].legend()
    #ax[2].legend()
    #plt.show()

    #ax1[0].hist(phas    , bins=200, range=[2250,2750],histtype='step',color='k', label='all'     )
    ax1[0].hist(phas    , bins=1000, range=[0,10000],histtype='step',color=colors[i], label='run{0:04d}'.format(run)     )
    ax1[1].hist(phas_wG , bins=1000, range=[0,10000],histtype='step',color=colors[i] )
    ax1[2].hist(phas_woG, bins=1000, range=[0,10000],histtype='step',color=colors[i] )
    ax1[2].set_xlabel('Pulse height [ch]')
    ax1[0].legend()
    #ax1[1].legend()
    #ax1[2].legend()
    #plt.show()

    bins_full_range=np.arange(0,35,0.1)
    #Fe_xrays=[6.46
    #ax2[0].hist(eners    , bins=bins_full_range,histtype='step',color=colors[i], label='all'     )
    #ax2[0].hist(eners    , bins=550, range=[4,15],histtype='step',color='k', label='all'     )
    ax2[0].hist(eners_wG , bins=bins_full_range,histtype='step',color=colors[i], label=label )
    ax2[1].hist(eners_wG, bins=300, range=[5.,35.],histtype='step',color=colors[i], label='w/ Gate (zoom)')
    #ax2[2].hist(eners_wG, bins=100, range=[4.9,6.9],histtype='step',color='r', label='w/ Gate (zoom)')
    #ax2[2].hist(eners_wG, bins=100, range=[2,4],histtype='step',color='r', label='w/ Gate (zoom)')
    #ax2[2].hist(eners_woG, bins=580, range=[2,60],histtype='step',color='b', label='w/o Gate')

    ax2[0].set_yscale('log')
    #ax2[1].set_yscale('log')
    ax2[0].grid(which='minor', axis='x', linestyle='dashed')
    ax2[0].grid(which='major', axis='x')
    ax2[0].set_xticks(np.arange(0,35,1), minor=True)
    #ax2[0].grid(which='minor', axis='x', linestyle='dashed')
    #ax2[0].grid(which='major', axis='x')
    #ax2[0].set_xticks(np.arange(0,50,1), minor=True)
    #ax2[0].axvline(6.40,color='b',linewidth=1)
    #ax2[0].axvline(7.06,color='b',linewidth=1)
    #ax2[0].axvline(14.4,color='b',linewidth=1)


    #ax2[1].text(xpos+2*dx,ypos,'muN'  ,color='k')

    #ax5[0].hist(eners_wG    , bins=ebins1,histtype='step',color=colors[i], label='run{0:04d}'.format(run)     )
    #ax5[1].hist(eners_wG    , bins=ebins2,histtype='step',color=colors[i], label='run{0:04d}'.format(run)     )
    ax5[0].hist(eners_wG    , bins=ebins1,histtype='step',color=colors[i], label=label     )
    ax5[1].hist(eners_wG    , bins=ebins2,histtype='step',color=colors[i], label=label     )
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
  ax5[1].legend(ncol=2)

  for muO in muonicO_xrays:
    ax2[0].axvline(muO,color='b',linewidth=1,alpha=alpha)
    ax2[1].axvline(muO,color='b',linewidth=1,alpha=alpha)
  for muX in muonicN_xrays:
    ax2[0].axvline(muX,color='c',linewidth=1,alpha=alpha)
    ax2[1].axvline(muX,color='c',linewidth=1,alpha=alpha)
  #for muOsub in muonicO_xrays_sub:
  #  ax2[1].axvline(muOsub,color='b',linewidth=1,linestyle='dashed')
  #for muFe in muonicFe_xrays:
  #  ax2[1].axvline(muFe,color='m',linewidth=1)
  #for muFesub in muonicFe_xrays_sub:
  #  ax2[1].axvline(muFesub,color='m',linewidth=1, linestyle='dashed')
  #for muAr in muonicAr_xrays:
  #  ax2[1].axvline(muAr,color='m',linewidth=1)
  #for muArsub in muonicAr_xrays_sub:
  #  ax2[1].axvline(muArsub,color='m',linewidth=0.8, linestyle='dashed')
  for muX in muonicAl_xrays:
    ax2[0].axvline(muX,color='m',linewidth=1,alpha=alpha)
    ax2[1].axvline(muX,color='m',linewidth=1,alpha=alpha)


  for muX in muonicCo_xrays:
    ax2[0].axvline(muX,color='g',linewidth=1,alpha=alpha)
    ax2[1].axvline(muX,color='g',linewidth=1,alpha=alpha)
  for muCu in muonicCu_xrays:
    ax2[0].axvline(muCu,color='gold',linewidth=1,alpha=alpha)
    ax2[1].axvline(muCu,color='gold',linewidth=1,alpha=alpha)
  for muC in muonicC_xrays:
    ax2[0].axvline(muC,color='k',linewidth=1,alpha=alpha)
    ax2[1].axvline(muC,color='k',linewidth=1,alpha=alpha)
  for muC in muonicC_xrays_sub:
    ax2[0].axvline(muC,color='k',linewidth=1,alpha=alpha)
    ax2[1].axvline(muC,color='k',linewidth=1,alpha=alpha)
  ax2[1].set_xlabel('energy [keV] (do not trust me!)')
  ax2[0].legend()
  #ax2[1].legend()
  #ax2[2].legend()
  x_lim=ax2[1].get_xlim()
  y_lim=ax2[1].get_ylim()
  #xpos=x_lim[0]+5
  xpos=30
  ypos=y_lim[1]-50
  dx=5
  dy=0.1*y_lim[1]
  ax2[1].text(xpos,ypos     ,'muC'  ,color='k')
  ax2[1].text(xpos,ypos-  dy,'muAl' ,color='m')
  ax2[1].text(xpos,ypos-2*dy,'muCu' ,color='gold')
  ax2[1].text(xpos,ypos-3*dy,'muCo' ,color='g')
  ax2[1].text(xpos,ypos-4*dy,'muO'  ,color='b')
  ax2[1].text(xpos,ypos-5*dy,'muN'  ,color='c')

  for cal_line in cal_lines:
    #ax2[1].axvline(cal_line,color='b',linewidth=1)#, linestyle='dashed')
    ax5[0].axvspan(cal_line+delta_E,cal_line-delta_E,color='b',alpha=0.3)#, linestyle='dashed')
    ax5[1].axvspan(cal_line+delta_E,cal_line-delta_E,color='b',alpha=0.3)#, linestyle='dashed')
    #ax2[2].axvline(cal_line,color='b',linewidth=1)#, linestyle='dashed')
    ax5[0].axvline(cal_line,color='b',linewidth=1)#, linestyle='dashed')
    ax5[1].axvline(cal_line,color='b',linewidth=1)#, linestyle='dashed')
  ax5[0].axvspan(noise_E+rfactor_noi*delta_E,noise_E-rfactor_noi*delta_E,color='gold',alpha=0.3)#, linestyle='dashed')
  ax5[1].axvspan(noise_E+rfactor_noi*delta_E,noise_E-rfactor_noi*delta_E,color='gold',alpha=0.3)#, linestyle='dashed')
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
  #pp.savefig(fig)
  fig1.tight_layout()
  #fig.savefig(outdir+outname)
  #pp.savefig(fig1)
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
    runlists= [4,5]
    runs= [4,5]
    runlists= [7,8,9,10,11]
    runs= [7,8,9,10,11]
    #runs= [8]
    #labels=['20.5MeV']
    runs= [10,11,12,13,14,15]
    labels=['21.5MeV','20.0MeV','19.0MeV','21.0MeV','22.0MeV','21.5MeV']
    runs=[16]
    labels=['20.5MeV']
    runs=[17,18,19,20]
    labels=['19.0MeV','20.0MeV','21.0MeV','19.5MeV']
    runs=[21,22,23,24]
    labels=['19.5MeV','19.0MeV','18.5MeV','20.0MeV']
    runs=[64,65,66,67,68,69]
    labels=['31.0MeV','29.0MeV','32.0MeV','33.0MeV','34.0MeV','35.0MeV']
    #runs=[25,26]
    #labels=['w/o slit','w/ slit']
    #runs=[104, 105, 107, 108, 109]
    #labels=['22.0MeV', '24.0MeV', '24.0MeV (flip)','21.0MeV (flip)','22.0MeV (flip)']
    #runs=[107, 108, 109, 110, 111]
    #labels=['24.0MeV (flip)','21.0MeV (flip)','22.0MeV (flip)','23.0MeV (flip)','25.0MeV (flip)']
    runs=[107, 109, 112]
    labels=['24.0MeV (flip)','22.0MeV (flip)','26.0MeV (flip)']
    #runs=[107, 110, 111]
    #labels=['24.0MeV (flip)','23.0MeV (flip)','25.0MeV (flip)']
    runs=[107,113]
    labels=['Battery','Li only']

    for run in runs :
      str_runlist+="{0:04d}_".format(run)
    #readfile(runs = runs, ENum = -1, datadir="C:\Data/RKN2024A/list/", outdir="C:\Data/RKN2024A/pdf/", outname="SDD_ene_MLF_run{0:s}compare.pdf".format(str_runlist),  cal_source="muAr", noise_E=9.0,rfactor_noi=2,  labels=labels)

    #readfile(runs = runs, ENum = -1, datadir="C:\Data/RKN2024A/list/", outdir="C:\Data/RKN2024A/pdf/", outname="SDD_ene_MLF_run{0:s}compare.pdf".format(str_runlist),  cal_source="muSi", noise_E=9.0,rfactor_noi=1,  labels=labels)
    readfile(runs = runs, ENum = -1, datadir="C:\Data/RKN2024A/list/", outdir="C:\Data/RKN2024A/pdf/", outname="SDD_ene_MLF_run{0:s}compare.pdf".format(str_runlist),  cal_source="muLi", noise_E=16.0,rfactor_noi=1,  labels=labels)




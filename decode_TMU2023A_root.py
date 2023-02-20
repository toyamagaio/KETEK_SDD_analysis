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
import ROOT

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

def readfile(inputfiles = ["list_000008_01-00000000001.bin"], ENum = -1, datadir="C:\Data/TMU2022C/list/", outdir="C:\Data/TMU2022C/pdf/", hfile_name="SDD_ene.root", cal_source=""):
  #bdata = open('C:\Data\list_000008_01-00000000001.bin','rb')

  hfile = ROOT.TFile(outdir+hfile_name,'RECREATE')
  tree = ROOT.TTree("tree", "SDD Analysis 2023");

  counter=0
  g_counter=0
  t_trigID =np.empty(1, dtype=np.int64)
  t_pha    =np.empty(1, dtype=np.int64)
  t_gate   =np.empty(1, dtype='bool' )
  t_pileup =np.empty(1, dtype='bool' )
  t_energy =np.empty(1, dtype=np.float64)
  t_time   =np.empty(1, dtype=np.float64)

  tree.Branch("trigID",t_trigID,"trigID/L")
  tree.Branch("pha"   ,t_pha   ,"pha/L"   )
  tree.Branch("gate"  ,t_gate  ,"gate/O"  )
  tree.Branch("pileup",t_pileup,"pileup/O")
  tree.Branch("energy",t_energy,"energy/D")
  tree.Branch("time"  ,t_time  ,"time/D"  )
  tree.Reset()
  
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

        data = bdata.read(4)#Realtime2
        idata=int.from_bytes(data, byteorder='big', signed=False)
        time2=idata*0.01 #0.01: 10ns->us
        #print('time2:',time2)
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
        puf=((idata>>15)&1) #pileup flag
        t_trigID[0] =trigid
        t_pha   [0] =pha
        t_gate  [0] =gate# lambda gate: True if (gate>0) else False
        t_pileup[0] =puf # lambda puf : True if (puf >0) else False
        t_energy[0] =ener
        t_time  [0] =time2


      tree.Fill()
      if counter%1000 ==0:
        print('{0:d} event analyzed'.format(counter))
        #print('gate {0}, tgate[0]:{1}'.format(gate, t_gate[0]) )
      if len(data) == 0:
        break
      if counter >= ENum & ENum > 0:
        break
      counter += 1
      if t_gate[0] > 0:
        g_counter += 1
    bdata.close

  print(outdir+hfile_name,'is saved.')
  total_event=tree.GetEntries()
  print('total_event: {0:d}, gated: {1:d}'.format(total_event, g_counter))
  tree.Write()
  hfile.Close()


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

    datadir='/Users/toyama/Documents/SDD_TMU2023A/list/'
    outdir ='/Users/toyama/Documents/SDD_TMU2023A/root/'

    #run=118 #Ag 30MeV/c
    #run=87 #Ag 30MeV/c

    for run in range(1,123):
      inputfile="list_{0:06d}_01-00000000001.bin".format(run)
      if os.path.isfile(datadir+inputfile):
        print(datadir+inputfile,' exists.')
        readfile(inputfiles = [inputfile], ENum = -1, datadir=datadir, outdir=outdir, hfile_name="SDD_MLF202302_run{0:04d}.root".format(run))
      else :
        print(datadir+inputfile,' does not exist.')

    inputfiles=[]
    ##all D2 run
    #irun=21
    #frun=102

    #irun=10
    #irun=8
    #frun=20


import sys
import os
import numpy as np
import ROOT
import math
import uproot
import glob
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle
from ROOT import TH1D, TH2D, TCanvas
from ROOT import std

gROOT.SetBatch(1)

def mk_TH1D(inputfiles=[], datadir='./', pdfname='hoge.pdf', orootname='hoge.root', pdfdir='pdf/', outdir='oroot/'):
  
    #fname=datadir+fname
    #print(fname)
    #inputfiles = glob.glob(fname)
    #print(inputfiles)

    if len(inputfiles)==0:
      print('no file named', fname)
      return

    tmax=100
    emax=50
    de  =5.748e-03 #calib factor
    nch_bin=3. #n ch per bin for energy hist
    ebins=np.arange(0,emax,nch_bin*de)
    t_th=25

    ofile = ROOT.TFile(outdir+orootname,'RECREATE')
    h2_et_all         =TH2D("h2_et_all"          ,"",len(ebins),0,ebins[-1],500,0,tmax)
    h2_et_gate        =TH2D("h2_et_gate"         ,"",len(ebins),0,ebins[-1],500,0,tmax)
    h2_et_ngate       =TH2D("h2_et_ngate"        ,"",len(ebins),0,ebins[-1],500,0,tmax)
    h2_et_pu          =TH2D("h2_et_pu"           ,"",len(ebins),0,ebins[-1],500,0,tmax)
    h2_et_tcut        =TH2D("h2_et_tcut"         ,"",len(ebins),0,ebins[-1],500,0,tmax)

    h_a_all         =TH1D("h_a_all"          ,"",8000,0,8000)
    h_a_gate        =TH1D("h_a_gate"         ,"",8000,0,8000)
    h_a_ngate       =TH1D("h_a_ngate"        ,"",8000,0,8000)
    h_a_pu          =TH1D("h_a_pu"           ,"",8000,0,8000)
    h_a_tcut        =TH1D("h_a_tcut"         ,"",8000,0,8000)

    h_e_all         =TH1D("h_e_all"          ,"",len(ebins),0,ebins[-1])
    h_e_gate        =TH1D("h_e_gate"         ,"",len(ebins),0,ebins[-1])
    h_e_ngate       =TH1D("h_e_ngate"        ,"",len(ebins),0,ebins[-1])
    h_e_pu          =TH1D("h_e_pu"           ,"",len(ebins),0,ebins[-1])
    h_e_tcut        =TH1D("h_e_tcut"         ,"",len(ebins),0,ebins[-1])

    h_t_all         =TH1D("h_t_all"          ,"",100,0,tmax)
    h_t_gate        =TH1D("h_t_gate"         ,"",100,0,tmax)
    h_t_ngate       =TH1D("h_t_ngate"        ,"",100,0,tmax)
    h_t_pu          =TH1D("h_t_pu"           ,"",100,0,tmax)
    h_t_tcut        =TH1D("h_t_tcut"         ,"",100,0,tmax)


    #file=uproot.open(inputfiles[0])#read single rootfile
    #t=file["Timepix"]
    for a_file in inputfiles:
      if os.path.isfile(datadir+a_file):
        print(a_file,' open')
      else:
        print(a_file,' does not exist. skipped')
        continue
      t=uproot.open(datadir+a_file)['tree']
      print(t.keys())
      trigID  =t['trigID'].array()
      phas    =t['pha'   ].array()
      gates   =t['gate'  ].array()
      pileups =t['pileup'].array()
      energy  =t['energy'].array()
      times   =t['time'  ].array()

      ###Event loop###
      for i, (tID, pha, gate, pu, ene, time) in enumerate(zip(trigID, phas, gates, pileups, energy, times)):
        if (i % 5000==0):
          print(i, tID, ene, time)

        h2_et_all .Fill(ene, time)
        h_a_all   .Fill(pha)
        h_e_all   .Fill(ene)
        h_t_all   .Fill(time)
        if pu:
          h2_et_pu .Fill(ene, time)
          h_a_pu   .Fill(pha)
          h_e_pu   .Fill(ene)
          h_t_pu   .Fill(time)
          continue
 
        if gate:
          h2_et_gate .Fill(ene, time)
          h_a_gate   .Fill(pha)
          h_e_gate   .Fill(ene)
          h_t_gate   .Fill(time)
        else:
          h2_et_ngate .Fill(ene, time)
          h_a_ngate   .Fill(pha)
          h_e_ngate   .Fill(ene)
          h_t_ngate   .Fill(time)
        if time < t_th:
          h2_et_tcut .Fill(ene, time)
          h_a_tcut   .Fill(pha)
          h_e_tcut   .Fill(ene)
          h_t_tcut   .Fill(time)

    #print
    pdfname=pdfdir+pdfname
    cv1=TCanvas("cv1","cv1",800,800)
    cv1.Divide(2,2)
    cv1.cd(1);gPad.SetLogy(1)
    h_a_all.Draw()
    cv1.cd(2);gPad.SetLogy(1)
    #h_a_pu .Draw()
    h_a_tcut .Draw()
    cv1.cd(3);gPad.SetLogy(1)
    h_a_gate.Draw()
    cv1.cd(4);gPad.SetLogy(1)
    h_a_ngate.Draw()

    cv1.Print(pdfname+"[")
    cv1.Print(pdfname)

    cv2=TCanvas("cv2","cv2",800,800)
    cv2.Divide(2,2)
    cv2.cd(1);gPad.SetLogy(1)
    h_e_all.Draw()
    cv2.cd(2);gPad.SetLogy(1)
    #h_e_pu .Draw()
    h_e_tcut .Draw()
    cv2.cd(3);gPad.SetLogy(1)
    h_e_gate.Draw()
    cv2.cd(4);gPad.SetLogy(1)
    h_e_ngate.Draw()
    cv2.Print(pdfname)

    cv3=TCanvas("cv3","cv3",800,800)
    cv3.Divide(2,2)
    cv3.cd(1);gPad.SetLogy(1)
    h_t_all.Draw()
    cv3.cd(2);gPad.SetLogy(1)
    #h_t_pu .Draw()
    h_t_tcut .Draw()
    cv3.cd(3);gPad.SetLogy(1)
    h_t_gate.Draw()
    cv3.cd(4);gPad.SetLogy(1)
    h_t_ngate.Draw()
    cv3.Print(pdfname)
    cv1.Print(pdfname+"]")

    h2_et_all   .Write()
    h2_et_gate  .Write()
    h2_et_ngate .Write()
    h2_et_tcut  .Write()

    h_a_all   .Write()
    h_a_gate  .Write()
    h_a_ngate .Write()
    h_a_tcut  .Write()

    h_e_all   .Write()
    h_e_gate  .Write()
    h_e_ngate .Write()
    h_e_tcut  .Write()

    h_t_all   .Write()
    h_t_gate  .Write()
    h_t_ngate .Write()
    h_t_tcut  .Write()
    ofile.Close()

if __name__ == "__main__":

    datadir ='/Users/toyama/Documents/SDD_TMU2023A/root/'
    pdfdir  ='/Users/toyama/Documents/SDD_TMU2023A/pdf/'
    outdir  ='/Users/toyama/Documents/SDD_TMU2023A/root/myroot/'
    inputfiles=[]

    irun=41 #1st production run (D2 target)
    frun=102#The last production run (D2 target)


    pdfname  ='SDD_MLF202302_run{0:04d}_{1:04d}.pdf'.format(irun,frun)
    orootname='SDD_MLF202302_run{0:04d}_{1:04d}.root'.format(irun,frun)

    for run in range(irun,frun):
      inputfiles.append("SDD_MLF202302_run{0:04d}.root".format(run))
    

    mk_TH1D(inputfiles=inputfiles, datadir=datadir, pdfname=pdfname, orootname=orootname, pdfdir=pdfdir, outdir=outdir)

    ####
 
    inputfiles=[]
    irun=103#1st Ag bg run w/ same beam condition of D2 target runs
    frun=121#The last 

    pdfname  ='SDD_MLF202302_run{0:04d}_{1:04d}.pdf'.format(irun,frun)
    orootname='SDD_MLF202302_run{0:04d}_{1:04d}.root'.format(irun,frun)

    for run in range(irun,frun):
      inputfiles.append("SDD_MLF202302_run{0:04d}.root".format(run))
    

    mk_TH1D(inputfiles=inputfiles, datadir=datadir, pdfname=pdfname, orootname=orootname, pdfdir=pdfdir, outdir=outdir)
 

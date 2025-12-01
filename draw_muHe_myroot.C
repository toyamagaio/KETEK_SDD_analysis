#include "Setting.cc"
const int NCanvas=6;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gaus_pol1bg(double *x, double *par) {
  //par[0]=area (gaus)
  //par[1]=mean (gaus)
  //par[2]=sigma (gaus)
  //par[3]=pol const
  //par[4]=pol tilt
  double val;
  double bg;
  double ga;

  ga=par[0]*TMath::Gaus(x[0],par[1],par[2],1);
  bg =  par[3]+ x[0]*par[4] ;
  val = ga + bg;
  return val;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gaus_pol0bg(double *x, double *par) {
  //par[0]=area (gaus)
  //par[1]=mean (gaus)
  //par[2]=sigma (gaus)
  //par[3]=pol const
  //par[4]=pol tilt
  double val;
  double bg;
  double ga;

  ga=par[0]*TMath::Gaus(x[0],par[1],par[2],1);
  bg =  par[3];
  val = ga + bg;
  return val;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gaus_pol3bg(double *x, double *par) {
  //par[0]=area (gaus)
  //par[1]=mean (gaus)
  //par[2]=sigma (gaus)
  //par[3]=pol tilt
  //par[4]=pol tilt
  //par[5]=pol tilt
  //par[6]=pol const
  double val;
  double bg;
  double ga;

  ga=par[0]*TMath::Gaus(x[0],par[1],par[2],1);
  bg =  par[3]*pow(x[0]-par[4],3.) + par[5]*(x[0]-par[4])+par[6];
  val = ga + bg;
  return val;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void draw_muHe_myroot(){
  int colors[]={1,2,4,6};
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);


  Setting *set= new Setting();
  TFile *file[2];
  string rootfiles[2];
  rootfiles[0]="./root/myroot/SDD_MLF202302_run0041_0102.root";
  rootfiles[1]="./root/myroot/SDD_MLF202302_run0103_0121.root";
  string pdfname="SDD_muHe.pdf";
  TH1D *h_e_tcut[2];
  TH1D *h_e_tcut_cl[2][4];
  TH1D *h_e_tcut_cl_bgresi[2][4];
  TH1D *h_e_tcut_cl_resi[2][4];

  double range[4][2];
  range[0][0]=6.8;
  range[0][1]=9.8;
  double nbg_min=6.9;
  double nbg_max=7.6;
  double nbg[2];

  TLine *tline = new TLine();
  set->SetTLine(tline, 14, 2, 2);


  double muBe_3_2=6.18;
  double muAl_5_4=10.66;
  range[1][0]=muBe_3_2-1.5;
  range[1][1]=muBe_3_2+1.5;
  range[2][0]=muAl_5_4-1.5;
  range[2][1]=muAl_5_4+1.5;
  range[3][0]=0.;
  range[3][1]=47.;

  //Get histograms from root file made by mk_TH1D.py
  for(int n=0;n<2;n++){
    file[n]=new TFile(rootfiles[n].c_str(),"readonly");

    h_e_tcut[n]=(TH1D*)file[n]->Get("h_e_tcut" );
    set->SetTH1(h_e_tcut[n],"","Energy [keV]","",colors[n],1000,0);   
    //h_e_tcut[n]->SetMinimum(0);
    h_e_tcut[n]->SetStats(0);
    nbg[n]=h_e_tcut[n]->Integral(h_e_tcut[n]->FindBin(nbg_min),h_e_tcut[n]->FindBin(nbg_max));
    cout<<"nbg"<<n<<": "<<nbg[n]<<"  in "<<nbg_min<<" -- "<<nbg_max<<endl;
  }
  h_e_tcut[1]->Scale(nbg[0]/nbg[1]);
  h_e_tcut[1]->SetLineWidth(1);

  double wbin=h_e_tcut[0]->GetBinWidth(1);
  cout<<"wbin: "<<wbin<<" keV"<<endl;

  for(int i=0;i<4;i++){
    h_e_tcut_cl[0][i]=(TH1D*)h_e_tcut[0]->Clone(Form("h_e_tcut_%d",i+1));
    h_e_tcut_cl[1][i]=(TH1D*)h_e_tcut[1]->Clone(Form("h_e_tcut_%d",i+1));
    h_e_tcut_cl[0][i]->GetXaxis()->SetRangeUser(range[i][0],range[i][1]);
    h_e_tcut_cl[1][i]->GetXaxis()->SetRangeUser(range[i][0],range[i][1]);
    h_e_tcut_cl[0][i]->SetMinimum(0);
    h_e_tcut_cl[1][i]->SetMinimum(0);

    h_e_tcut_cl_bgresi[0][i]=(TH1D*)h_e_tcut[0]->Clone(Form("h_e_tcut_%d_%d_bgresi",1,i+1));
    h_e_tcut_cl_bgresi[1][i]=(TH1D*)h_e_tcut[1]->Clone(Form("h_e_tcut_%d_%d_bgresi",2,i+1));
    h_e_tcut_cl_bgresi[0][i]->GetXaxis()->SetRangeUser(range[i][0],range[i][1]);
    h_e_tcut_cl_bgresi[1][i]->GetXaxis()->SetRangeUser(range[i][0],range[i][1]);

    h_e_tcut_cl_resi[0][i]=(TH1D*)h_e_tcut[0]->Clone(Form("h_e_tcut_%d_%d_resi",1,i+1));
    h_e_tcut_cl_resi[1][i]=(TH1D*)h_e_tcut[1]->Clone(Form("h_e_tcut_%d_%d_resi",2,i+1));
    h_e_tcut_cl_resi[0][i]->GetXaxis()->SetRangeUser(range[i][0],range[i][1]);
    h_e_tcut_cl_resi[1][i]->GetXaxis()->SetRangeUser(range[i][0],range[i][1]);
  }
  h_e_tcut_cl[0][3]->SetMinimum(0.8);
  h_e_tcut_cl[1][3]->SetMinimum(0.8);
  //h_e_tcut_cl[1]->GetXaxis()->SetRangeUser(muBe_3_2-1.5,muBe_3_2+1.5);
  //h_e_tcut_cl[2]->GetXaxis()->SetRangeUser(muAl_5_4-1.5,muAl_5_4+1.5);

  for(int i=0;i<3;i++){
    int maxbin=h_e_tcut_cl[0][i]->GetMaximumBin();
    double maxcont = h_e_tcut_cl[0][i]->GetBinContent(maxbin);
    double maxx    = h_e_tcut_cl[0][i]->GetBinCenter(maxbin);
    double maxbg   = h_e_tcut[1]->GetBinContent(h_e_tcut[1]->FindBin(maxx));
    if(maxbg > maxcont)
      h_e_tcut_cl[0][i]->GetYaxis()->SetRangeUser(0,1.2*maxbg);
    else
      h_e_tcut_cl[1][i]->GetYaxis()->SetRangeUser(0,1.2*maxcont);
  }


  ///////////
  //Fitting//
  ///////////
  
  TF1 *f_sum_Ag[3], *f_gau_Ag[3], *f_bg_Ag[3];
  TF1 *f_sum_D2[3], *f_gau_D2[3], *f_bg_D2[3];
  double par_Ag[3][5], par_D2[3][5];
  double err_Ag[3][5], err_D2[3][5];
  for(int i=0;i<3;i++){
    f_sum_Ag[i]= new TF1(Form("f_sum_Ag_%d",i+1),gaus_pol1bg ,range[i][0], range[i][1], 5);
    f_gau_Ag[i]= new TF1(Form("f_gau_Ag_%d",i+1),"gausn(0)"  ,range[i][0], range[i][1]);
    f_bg_Ag [i]= new TF1(Form("f_bg_Ag_%d" ,i+1),"[0]+[1]*x" ,range[i][0], range[i][1]);
    f_sum_D2[i]= new TF1(Form("f_sum_D2_%d",i+1),gaus_pol1bg ,range[i][0], range[i][1], 5);
    f_gau_D2[i]= new TF1(Form("f_gau_D2_%d",i+1),"gausn(0)"  ,range[i][0], range[i][1]);
    f_bg_D2 [i]= new TF1(Form("f_bg_D2_%d" ,i+1),"[0]+[1]*x" ,range[i][0], range[i][1]);

    set->SetTF1(f_sum_Ag[i],4,1,1);
    set->SetTF1(f_gau_Ag[i],4,1,1);
    set->SetTF1(f_bg_Ag [i],4,1,1);
    set->SetTF1(f_sum_D2[i],2,1,1);
    set->SetTF1(f_gau_D2[i],2,1,1);
    set->SetTF1(f_bg_D2 [i],2,1,1);

    h_e_tcut[0] -> Fit(f_bg_D2[i],"0QR","",range[i][0],range[i][0]+1);
    h_e_tcut[1] -> Fit(f_bg_Ag[i],"0QR","",range[i][0],range[i][0]+1);

    int nbins =h_e_tcut_cl_bgresi[1][i]->GetNbinsX();
    for(int j=1;j<nbins;j++){
      double x=h_e_tcut_cl_bgresi[0][i]->GetBinCenter(j);
      double bg_D2=f_bg_D2[i]->Eval(x);
      double bg_Ag=f_bg_Ag[i]->Eval(x);
      h_e_tcut_cl_bgresi[0][i]->SetBinContent(j,h_e_tcut_cl_bgresi[0][i]->GetBinContent(j)-bg_D2);
      h_e_tcut_cl_bgresi[1][i]->SetBinContent(j,h_e_tcut_cl_bgresi[1][i]->GetBinContent(j)-bg_Ag);
    }
    h_e_tcut_cl_bgresi[0][i]->Fit(f_gau_D2[i],"0QR","",range[i][0]+1.,range[i][1]-1.);
    h_e_tcut_cl_bgresi[1][i]->Fit(f_gau_Ag[i],"0QR","",range[i][0]+1.,range[i][1]-1.);
    f_bg_D2[i] ->GetParameters(&par_D2[i][3]);
    f_bg_Ag[i] ->GetParameters(&par_Ag[i][3]);
    f_gau_D2[i]->GetParameters(&par_D2[i][0]);
    f_gau_Ag[i]->GetParameters(&par_Ag[i][0]);

    f_sum_D2[i]->SetParameters(&par_D2[i][0]);
    f_sum_Ag[i]->SetParameters(&par_Ag[i][0]);

    h_e_tcut[0]->Fit(f_sum_D2[i],"0QR","",range[i][0],range[i][1]);
    h_e_tcut[1]->Fit(f_sum_Ag[i],"0QR","",range[i][0],range[i][1]);

    f_sum_D2[i]->GetParameters(&par_D2[i][0]);
    f_sum_Ag[i]->GetParameters(&par_Ag[i][0]);
    if(i==0){
      f_sum_D2[i]->FixParameter(1, par_Ag[i][1]);
      f_sum_D2[i]->FixParameter(2, par_Ag[i][2]);
      f_sum_D2[i]->FixParameter(3, par_Ag[i][3]);
      f_sum_D2[i]->FixParameter(4, par_Ag[i][4]);
      h_e_tcut[0]->Fit(f_sum_D2[i],"0QR","",range[i][0],range[i][1]);
      f_sum_D2[i]->GetParameters(&par_D2[i][0]);
    }


    for(int j=0;j<5;j++){
      err_D2[i][j]=f_sum_D2[i]->GetParError(j);
      err_Ag[i][j]=f_sum_Ag[i]->GetParError(j);
    }

    nbins =h_e_tcut_cl_resi[1][i]->GetNbinsX();
    for(int j=1;j<nbins;j++){
      double x=h_e_tcut_cl_resi[0][i]->GetBinCenter(j);
      double ft_D2=f_sum_D2[i]->Eval(x);
      double ft_Ag=f_sum_Ag[i]->Eval(x);
      h_e_tcut_cl_resi[0][i]->SetBinContent(j,h_e_tcut_cl_resi[0][i]->GetBinContent(j)-ft_D2);
      h_e_tcut_cl_resi[1][i]->SetBinContent(j,h_e_tcut_cl_resi[1][i]->GetBinContent(j)-ft_Ag);
    }
 
    cout<<"Range:"<<range[i][0]<<" -- "<<range[i][1]<<endl;
    cout<<"Gauss center D2: "<<par_D2[i][1]<<" +/- "<<err_D2[i][1]<<endl;
    cout<<"Gauss center Ag: "<<par_Ag[i][1]<<" +/- "<<err_Ag[i][1]<<endl;
    cout<<"Gauss width D2: " <<par_D2[i][2]<<" +/- "<<err_D2[i][2]<<endl;
    cout<<"Gauss width Ag: " <<par_Ag[i][2]<<" +/- "<<err_Ag[i][2]<<endl;
    cout<<"==========="<<endl;
  }
  

  TCanvas *c[NCanvas];
  for(int i=0;i<NCanvas;i++){
    c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1), 1600,900);
  }

  c[0]->cd();
  //h_e_tcut_cl[0][0]->Draw("hist");
  h_e_tcut_cl[1][0]->Draw("hist");
  h_e_tcut_cl[0][0]->Draw("samehist");
  //f_gau_Ag[1]->Draw("same");

  c[1]->Divide(2,2);
  c[1]->cd(1);
  h_e_tcut_cl[1][0]->Draw("hist");
  h_e_tcut_cl[0][0]->Draw("samehist");
  c[1]->cd(2);
  h_e_tcut_cl[0][0]->Draw();//f_sum_D2[0]->Draw("same");
  c[1]->cd(3);
  h_e_tcut_cl[1][0]->Draw("hist");//f_sum_Ag[0]->Draw("same");
  c[1]->cd(4);
  h_e_tcut_cl[0][0]->Draw();//f_sum_D2[0]->Draw("same");
  h_e_tcut_cl[1][0]->Draw("samehist");//f_sum_Ag[0]->Draw("same");

  c[2]->Divide(2,2);
  c[2]->cd(1);
  h_e_tcut_cl[1][1]->Draw("hist");
  h_e_tcut_cl[0][1]->Draw("samehist");
  c[2]->cd(2);
  h_e_tcut_cl[0][1]->Draw();//f_sum_D2[1]->Draw("same");
  c[2]->cd(3);
  h_e_tcut_cl[1][1]->Draw();//f_sum_Ag[1]->Draw("same");
  c[2]->cd(4);
  h_e_tcut_cl[0][1]->Draw();//f_sum_D2[1]->Draw("same");
  h_e_tcut_cl[1][1]->Draw("samehist");//f_sum_Ag[1]->Draw("same");
  //c[2]->cd(5);
  //h_e_tcut_cl_bgresi[0][1]->Draw();f_gau_D2[1]->Draw("same");
  //c[2]->cd(6);
  //h_e_tcut_cl_bgresi[1][1]->Draw();f_gau_Ag[1]->Draw("same");

  c[3]->Divide(2,2);
  c[3]->cd(1);
  h_e_tcut_cl[1][2]->Draw("hist");
  h_e_tcut_cl[0][2]->Draw("samehist");
  c[3]->cd(2);
  h_e_tcut_cl[0][2]->Draw();//f_sum_D2[2]->Draw("same");
  c[3]->cd(3);
  h_e_tcut_cl[1][2]->Draw();//f_sum_Ag[2]->Draw("same");
  c[3]->cd(4);
  h_e_tcut_cl[0][2]->Draw();//f_sum_D2[2]->Draw("same");
  h_e_tcut_cl[1][2]->Draw("samehist");//f_sum_Ag[2]->Draw("same");
  //c[3]->cd(5);
  //h_e_tcut_cl_bgresi[0][2]->Draw();f_gau_D2[2]->Draw("same");
  //c[3]->cd(6);
  //h_e_tcut_cl_bgresi[1][2]->Draw();f_gau_Ag[2]->Draw("same");

  c[4]->Divide(2,2);
  c[4]->cd(1);
  h_e_tcut_cl_resi[0][0]->Draw("hist");
  h_e_tcut_cl_resi[1][0]->Draw("samehist");
  tline->DrawLine(range[0][0],0,range[0][1],0);
  c[4]->cd(2);
  h_e_tcut_cl_resi[0][0]->Draw("hist");
  tline->DrawLine(range[0][0],0,range[0][1],0);
  c[4]->cd(3);
  h_e_tcut_cl_resi[1][0]->Draw("hist");
  tline->DrawLine(range[0][0],0,range[0][1],0);

  //c[5]->Divide(1,2);
  c[5]->cd();gPad->SetLogy(1);//gPad->SetGrid(0);
  h_e_tcut_cl[1][3]->Draw("hist");
  h_e_tcut_cl[0][3]->Draw("samehist");
  gPad->Update();
  double yax_max = (double)h_e_tcut_cl[1][3]->GetYaxis()->GetXmax();
  cout<<"y max display: "<<yax_max<<endl;
  yax_max=0.48e5;
  double yax_min=0.8;

  //int npeak=30;
  //TSpectrum *s=new TSpectrum(npeak);
  //s->Search(h_e_tcut_cl[1][3], 2,"new");
  //npeak=s->GetNPeaks();
  //TList *functions = h_e_tcut_cl[1][3]->GetListOfFunctions();
  //TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
  //cout<<"npeak:"<<npeak<<endl;
  ////double *peakposx=s->GetPositionX();
  //double *peakposx=pm->GetX();
  //for(int i=0;i<npeak;i++){
  //  //double a=peakposx[i];
  //  //int bin=1+int(a+0.5);
  //  //double pos=h_e_tcut_cl[1][3]->GetBinCenter(bin);
  //  cout<<"peak pos: "<<peakposx[i]<<endl;
  //}


  double muBe[]     ={33.378, 6.18, 2.16};//2-1, 3-2, 4-3
  double muBe_Ly[]  ={39.558,41.21,43.27};//3-1, 4-1, 5-1
  double muBe_sub[] ={8.343 ,9.344};//4-2, 5-2
  double muAg[]     ={45.86, 29.75, 20.39, 14.58, 10.79, 8.21, 6.39};//7-6,8-7,...,13-12
  double muAg_Ly[] ={34.97, 25.37, 18.99};//10-8, 11-9, 12-10
  double muAg_sub[]={-100 ,-100};//
  double muAl[]    ={23.04, 10.66, 5.79};//4-3, 5-4, 6-5
  double muAl_Ly[] ={33.71 ,39.50, 42.99};//5-3, 6-3, 7-3
  double muAl_sub[]={16.46,-100.0,19.95,22.22};//6-4 , 7-4, 8-4
  double muSi[]    ={26.73, 12.37, 6.72};//4-3, 5-4, 6-5
  double muSi_Ly[] ={39.10, 45.82, 49.87};//5-3, 6-3, 7-3
  double muSi_sub[]={19.09, 23.14, 25.77};//6-4 , 7-4, 8-4
  double muC[]     ={13.95,4.88,2.26};//3-2, 4-3, 5-4
  double muC_Ly[]  ={18.83,21.09,22.32};//3-2, 4-3, 5-4
  double muC_sub[] ={-100,-100.};//3-2, 4-3, 5-4

  double Pd_Ka=21.177; //elec X-ray (muAg as psued Pd)
  double Pd_Kb=23.830; //elec X-ray (muAg as psued Pd)
  double Ag_Ka=21.990;//Ka1=22.163, Ka2=21.99
  double Ag_Kb=24.942;//Kb1=24.942
  //double 
  TLine *tline_Be[3],*tline_Ag[3],*tline_Al[3], *tline_Si[3], *tline_C[3];
  TLine *tline_X=new TLine();
  tline_X->DrawLine(Pd_Ka,yax_min,Pd_Ka,yax_max);
  tline_X->DrawLine(Pd_Kb,yax_min,Pd_Kb,yax_max);
  tline_X->DrawLine(Ag_Ka,yax_min,Ag_Ka,yax_max);
  tline_X->DrawLine(Ag_Kb,yax_min,Ag_Kb,yax_max);
  for(int i=0;i<3;i++){
    tline_Be[i]= new TLine();
    tline_Ag[i]= new TLine();
    tline_Al[i]= new TLine();
    tline_Si[i]= new TLine();
    tline_C[i] = new TLine();
    set->SetTLine(tline_Be[i], kBlue+2  , i+1, 1);
    set->SetTLine(tline_Ag[i], kGreen+3 , i+1, 1);
    set->SetTLine(tline_Al[i], kYellow+3, i+1, 1);
    set->SetTLine(tline_Si[i], kViolet+1, i+1, 1);
    set->SetTLine(tline_C[i] , kOrange+2, i+1, 1);
  }
  for(int i=0;i<3;i++){
    tline_Ag[0]->DrawLine(muAg[i],yax_min,muAg[i],yax_max);
    tline_Al[0]->DrawLine(muAl[i],yax_min,muAl[i],yax_max);
    tline_Si[0]->DrawLine(muSi[i],yax_min,muSi[i],yax_max);
  }
  tline_Be[0]->DrawLine(muBe[0],yax_min,muBe[0],yax_max);
  tline_Be[0]->DrawLine(muBe[1],yax_min,muBe[1],yax_max);
  tline_C [0]->DrawLine(muC [0],yax_min,muC [0],yax_max);
  tline_Ag[0]->DrawLine(muAg[3],yax_min,muAg[3],yax_max);
  //tline_Ag[0]->DrawLine(muAg[4],yax_min,muAg[4],yax_max);
  //for(int i=0;i<3;i++){
  //  tline_Be[1]->DrawLine(muBe_Ly[i],yax_min,muBe_Ly[i],yax_max);
  //  tline_Ag[1]->DrawLine(muAg_Ly[i],yax_min,muAg_Ly[i],yax_max);
  //  tline_Al[1]->DrawLine(muAl_Ly[i],yax_min,muAl_Ly[i],yax_max);
  //  tline_Si[1]->DrawLine(muSi_Ly[i],yax_min,muSi_Ly[i],yax_max);
  //  tline_C [1]->DrawLine(muC_Ly [i],yax_min,muC_Ly [i],yax_max);
  //}
  tline_Be[1]->DrawLine(muBe_Ly[0],yax_min,muBe_Ly[0],yax_max);
  tline_Ag[1]->DrawLine(muAg_Ly[1],yax_min,muAg_Ly[1],yax_max);
  tline_Ag[1]->DrawLine(muAg_Ly[2],yax_min,muAg_Ly[2],yax_max);
  for(int i=0;i<2;i++){
    tline_Be[2]->DrawLine(muBe_sub[i],yax_min,muBe_sub[i],yax_max);
    tline_Ag[2]->DrawLine(muAg_sub[i],yax_min,muAg_sub[i],yax_max);
    tline_Al[2]->DrawLine(muAl_sub[i],yax_min,muAl_sub[i],yax_max);
    //tline_Si[2]->DrawLine(muSi_sub[i],yax_min,muSi_sub[i],yax_max);
  }


  c[0]->Print(Form("%s[",pdfname.c_str()));
  for(int i=0;i<NCanvas;i++)
    c[i]->Print(Form("%s",pdfname.c_str()));
  c[NCanvas-1]->Print(Form("%s]",pdfname.c_str()));
  
}

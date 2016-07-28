//##############################
//### this code is modified ####
//### dliu13@mail.ustc.edu.cn###
//##############################


#include "TFile.h"
//#include "TFolder.h"
#include "TTree.h"
#include "TEventList.h"
#include "TLeaf.h"
#include "TLorentzVector.h"
#include "TH2D.h"
#include "TH1.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

using namespace std;
int Ana(const char *filename=0, const char* outdir=0, TFile* fileout=0);
double GetEnergy(int run);
//double FitSpectrum(TTree *&dataraw, double beame, const char* namesfx=0); // cb signal, gaus back
//double FitSpectrum1(TTree *&dataraw, double beame, const char* namesfx=0); // cb signal, poly bg
//double FitSpectrum2(TTree *&dataraw, double beame, const char* namesfx=0);
  double FitSpectrum3(TTree *dataraw, double beame, double &err, const char* namesfx=0);
  double FitSpectrum_rg(TTree *dataraw, double beame, const char* namesfx=0);
  double FitSpectrum_mcshape(TTree *&dataraw, TH1D* hsig, TH1D* hbck, double beame, double &err, const char* namesfx, TFile* fileout, double* par_pool=0, int rw = 10);
  double FitSpectrum_smcs(TTree *dataraw, TH1D* hsig, double beame, const char* namesfx=0, double *par_pool=0);
  double FitSpectrum_smcsAbckgaus(TTree *dataraw, TH1D* hsig, double beame, const char* namesfx=0, double *par_pool=0);
  double FitSpectrum_bmcs(TTree *dataraw, TH1D* hbck, double beame, const char* namesfx=0, double *par_pool=0);
  double FitSpectrum_bmcsAsiggaus(TTree *dataraw, TH1D* hbck, double beame, const char* namesfx=0, double *par_pool=0);
//double FitSpectrum3(TEventList *dataraw, double beame, const char* namesfx=0);
double FitSpectrum4(TTree *&dataraw, double beame, double &err, const char* namesfx=0, TFile* fileout=0); // fit with MC shape
double FitSpectrum_mcshape_fitpar(TTree *&dataraw, TH1D* hsig, TH1D* hbck, double beame, double &err, const char* namesfx, TFile* fileout);
const char* getPureName(const char* name);
const char* getPureName2(const char* name);
double getEne(const char* name);


int main(int argc, char** argv)
{
  //if (argc==1) Ana(argv[1]);
  //if (argc>1)  Ana(argv[1],argv[2]);
  TFile *fileout = new TFile("output/output_noF_1.root","recreate");
  Ana("lalala",0,fileout);
  //for (int i=1; i<argc; i++){
  //  Ana(argv[i],0,fileout);
  //}
  fileout->Write();
  delete fileout;
  return 0;
}

int Ana(const char *filename, const char* outdir, TFile *fileout)
{
  //TFile *file = new TFile("mc_KPI_22324.root");
  TFile *file = new TFile("output/output_openang.root");
  std::cout<<"File name is "<<file->GetName()<<std::endl;
  if (file==0) return -1;
//TTree *tree = (TTree*)file->Get("TwoProng");
//if (tree==0) return -2;

  double energys[22] = {
  2.0000, 2.0500, 2.1000, 2.125, 2.1500, 2.1750,
  2.2000, 2.2324, 2.3094, 2.3864, 2.3960,
  2.5000, 2.6444, 2.6464, 2.7000, 2.8000,
  2.9000, 2.9500, 2.9810, 3.0000, 3.0200,
  3.0800
  };
  double ene[22];
  double pcut[22];
  double epcut[22];
  double thecut[22];
  double m_pcut;
  double m_epcut;
  double m_thecut;
  ifstream cutin("cutpar");
  char line[1000];
  if (cutin.is_open()) cutin.getline(line,1000);
  int iene=0;
  // set E/p and p cut, mark: set cut
  while ( !cutin.eof()){
    cutin.getline(line,1000);
    istringstream iss;
    iss.str(line);
    iss>>ene[iene]>>pcut[iene]>>epcut[iene];
    iene++;
    if (iene==22) break;
  }

  double pi = TMath::Pi();
  double mka = 0.493677;
  //double mka = 0.13957;
  //double mka = 0.1057;

  //TFile *fileout = new TFile(name,"recreate");
  const char *pureName = getPureName2(filename);
  std::cout<<"Pure Name: "<< pureName <<std::endl;
  char name1[1000];
  //sprintf(name1,"output/%s.root",pureName);
  sprintf(name1,"output/%s.root",pureName);
  TFile *dir = new TFile(name1,"recreate");
  
  double totp, costheta1, costheta2;
  double theta;
  double p1,p2;
  TTree *elist = new TTree("elist","elist");
  elist->Branch("p1",&p1,"p1/D");
//TTree *elist2 = new TTree("elist2","elist2");
//elist2->Branch("p1",&p1,"p1/D");
//TTree *elist3 = new TTree("elist3","elist3");
//elist3->Branch("p1",&p1,"p1/D");
//TTree *elist4 = new TTree("elist4","elist4");
//elist4->Branch("p1",&p1,"p1/D");
   
 
  fileout->cd();

  double pexp[22]={
  
  };
  double psigma[22];
  for (int i=0;i<22;i++){
    pexp[i] = sqrt(pow(energys[i]/2,2)-pow(mka,2));
    psigma[i] = (pcut[i]-pexp[i])/3;
    //cout << "Ecm = "<< energys[i]<<", pexp = "<< pexp[i]<<", sigma = "<< (pcut[i]-pexp[i])/3<<endl;
  }

  int count0=0,count1=0,count2=0,count3=0;
  int count[10]={0};
  double par_pool[100];

  for (int iene=0; iene<22;iene++){
    char name[100];
    sprintf(name,"KK_%d",(int)(energys[iene]*10000));
    TDirectory* dir_data = (TDirectory*)file->Get(name);
    TTree* var_data = (TTree*) dir_data->Get("vars");
    var_data->SetBranchAddress("p1",&p1);
    var_data->SetBranchAddress("p2",&p2);
    var_data->SetBranchAddress("theta",&theta);
    m_pcut = pcut[iene];

    elist->Reset();
    //elist2->Reset();
    double entries = var_data->GetEntries();
    for (int ien=0;ien<entries;ien++){
      var_data->GetEntry(ien);
      if (theta>179 && p2<m_pcut && p1>0) elist->Fill();
    }

    
    char evt_flt[1000];
    double p_low, p_up;
    p_low = pexp[iene]-0.1;
    p_up  = pexp[iene]+0.2;
    if (energys[iene]<2.6) {
      p_low = pexp[iene]-0.05;
      p_up  = pexp[iene]+0.15;
    }
    sprintf(evt_flt,"theta>179 && p1>%f && p1<%f && p2<%f ",p_low,p_up,pcut[iene]);
    TH1D* hp_mckk  = new TH1D("hp_mckk","",100,p_low,p_up); //###########
    TH1D* hp_mcmu  = new TH1D("hp_mcmu","",100,p_low,p_up);
    sprintf(name,"mcKK_ecal_%03d_dst",iene+1);
    TDirectory* dir_mckk = (TDirectory*)file->Get(name);
    TTree* var_mckk = (TTree*) dir_mckk->Get("vars");
    sprintf(name,"mcdimu_%03d_%.4f",iene+1,energys[iene]);
    if (iene==3) sprintf(name,"mcdimu_%.4f",energys[iene]);
    if (iene>3) sprintf(name,"mcdimu_%03d_%.4f",iene,energys[iene]);
    TDirectory* dir_mcmu = (TDirectory*)file->Get(name);
    TTree* var_mcmu = (TTree*) dir_mcmu->Get("vars");
    var_mckk->Draw("p1>>hp_mckk",evt_flt);
    var_mcmu->Draw("p1>>hp_mcmu",evt_flt);


  //sprintf(name,"mcKK_ecal_%03d_dst",iene+1);
  //TDirectory* dir_mckk = (TDirectory*)file->Get(name);
  //TTree* var_mckk = (TTree*) dir_mckk->Get("vars");
  //var_mckk->SetBranchAddress("p1",&p1);
  //var_mckk->SetBranchAddress("p2",&p2);
  //var_mckk->SetBranchAddress("theta",&theta);
  //elist3->Reset();
  //elist4->Reset();
  //entries = var_mckk->GetEntries();
  //for (int ien=0;ien<entries;ien++){
  //  var_mckk->GetEntry(ien);
  //  if (theta>179 && p2<m_pcut) elist3->Fill();
  //  if (theta>179.5 && p2<m_pcut) elist4->Fill();
  //}

    ofstream cutflow("cutflow2",std::ios::app);
    //ofstream cutflow("cutflow_cmp665andp01",std::ios::app);
    cutflow<<energys[iene]<<" change fitting shape"<<std::endl;
    cutflow<<"Initial size      :"<<var_data->GetEntries()<<std::endl;
    cutflow<<"After theta cut   :"<<var_data->Draw("","theta<179")<<std::endl;
    cutflow<<"p2-exp<3 sigma    :"<<elist->GetEntries()<<std::endl;
    
    double err;
    sprintf(name,"KK_%d",(int)(energys[iene]*10000));
    //double sigp1=FitSpectrum3(elist, energys[iene],err, name);
    // 4 parameters in pool, smear par for signal, mean, sigma; smear par for background, mean, sigma
    double sigp1=FitSpectrum_mcshape(elist, hp_mckk, hp_mcmu, energys[iene],err, name, fileout, par_pool, 1);
    cutflow<<"p1 fit signal pre :"<<sigp1<<" +/-"<<err <<std::endl;
    //double sigp1=FitSpectrum_mcshape_fitpar(elist, hp_mckk, hp_mcmu, energys[iene],err, name, fileout);
    //double sigp1=FitSpectrum4(elist, energys[iene], name,fileout);
  //double err;
    //sprintf(name,"KK_%d",(int)(energys[iene]*10000));
    //double sigp1=FitSpectrum4(elist, energys[iene],err, name, fileout);
    //cutflow<<"p1 fit signal pre :"<<sigp1<<" +/-"<<err <<std::endl;
  //double sigp1_rg=FitSpectrum_rg(elist, energys[iene], name);
  //cutflow<<"p1 fit signal crg :"<<sigp1_rg<< " relative changes "<<(sigp1_rg-sigp1)/sigp1 <<std::endl;
    
  //double par_pool2[100];
  //par_pool2[0]  = par_pool[0];
  //par_pool2[1]  = par_pool[1];
  //par_pool2[2]  = par_pool[2];
  //par_pool2[3]  = par_pool[3];
  //par_pool2[10] = par_pool[10];
  //par_pool2[11] = par_pool[11];
  //par_pool2[12] = par_pool[12];
  //par_pool2[13] = par_pool[13];
  //sprintf(name,"KK_repeat_%d",(int)(energys[iene]*10000));
  //double sigp1_rep=FitSpectrum_mcshape(elist, hp_mckk, hp_mcmu, energys[iene],err, name, fileout, par_pool, 10);
  //cutflow<<"p1 fit signal repe:"<<sigp1_rep << " relative changes "<<(sigp1_rep-sigp1)/sigp1<<std::endl;
  //par_pool2[2]  = par_pool[2]+par_pool[3];
  //sprintf(name,"KK_sigsmear_%d",(int)(energys[iene]*10000));
  //// 2 parameters in pool, smear par for signal, mean, sigma; 
  //double sigp1_smc=FitSpectrum_mcshape(elist, hp_mckk, hp_mcmu, energys[iene],err, name, fileout, par_pool2, 10);
  //cutflow<<"p1 fit signal smcs:"<<sigp1_smc << " relative changes "<<(sigp1_smc-sigp1)/sigp1 << " relative changes2 "<<(sigp1_smc-sigp1_rep)/sigp1_rep<<std::endl;
  //par_pool2[2]  = par_pool[2];
  //par_pool2[12] = par_pool[12]+par_pool[13];
  //sprintf(name,"KK_bcksmear_%d",(int)(energys[iene]*10000));
  //// 2 parameters in pool,  smear par for background, mean, sigma
  //double sigp1_bmc=FitSpectrum_mcshape(elist, hp_mckk, hp_mcmu, energys[iene],err, name, fileout, par_pool2, 10);
  //cutflow<<"p1 fit signal bmcs:"<<sigp1_bmc  << " relative changes "<<(sigp1_bmc-sigp1)/sigp1 << " relative changes2 "<<(sigp1_bmc-sigp1_rep)/sigp1_rep<<std::endl;
    
    // 2 parameters in pool, smear par for signal, mean, sigma; 
    double sigp1_smc=FitSpectrum_smcs(elist, hp_mckk, energys[iene], name);
  //double sigp1_smc=FitSpectrum_smcsAbckgaus(elist, hp_mckk, energys[iene], name);
    cutflow<<"p1 fit signal smcs:"<<sigp1_smc << " relative changes "<<(sigp1_smc-sigp1)/sigp1<<std::endl;
    // 2 parameters in pool,  smear par for background, mean, sigma
    double sigp1_bmc=FitSpectrum_bmcs(elist, hp_mcmu, energys[iene], name);
    //double sigp1_bmc=FitSpectrum_bmcsAsiggaus(elist, hp_mcmu, energys[iene], name);
    cutflow<<"p1 fit signal bmcs:"<<sigp1_bmc << " relative changes "<<(sigp1_bmc-sigp1)/sigp1<<std::endl;
    
    delete hp_mckk;
    delete hp_mcmu;

  //sprintf(name,"KK_%d_179.5",(int)(energys[iene]*10000));
  //double sigp2=FitSpectrum3(elist2, energys[iene],err, name);
  //cutflow<<"p1 fit signal pre :"<<sigp2<<" +/-"<<err <<std::endl;
  //sprintf(name,"mcKK_%d_179",(int)(energys[iene]*10000));
  //double sigp3=FitSpectrum3(elist3, energys[iene],err, name);
  //cutflow<<"p1 fit signal pre :"<<sigp3<<" +/-"<<err <<std::endl;
  //sprintf(name,"mcKK_%d_179.5",(int)(energys[iene]*10000));
  //double sigp4=FitSpectrum3(elist4, energys[iene],err, name);
  //cutflow<<"p1 fit signal pre :"<<sigp4<<" +/-"<<err <<std::endl;
  //cutflow<<"uncertainty "<< (sigp1/sigp3 - sigp2/sigp4)/(sigp1/sigp3)<<endl;
  }
  // finish program here
  return 0;
 
}

#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include "TPaveText.h"
#include "TMath.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooCBShape.h"
#include "RooGenericPdf.h"
#include "RooChebychev.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooMsgService.h"
#include "TGaxis.h"
#include "TRandom.h"
#include <ctime>
using namespace RooFit;


double exp1(double *x, double *par)
{
    return TMath::Exp(par[0]*(x[0]-par[1]))+par[2];
}




double FitSpectrum3(TTree *dataraw, double beame,double &err, const char* namesfx)
{
   TF1 fsigma("fsigma","-0.00705519+0.00600716*x",2.0,3.1);
   double sigma_ini = fsigma.Eval(beame);
   TCanvas *c1 = new TCanvas("","",800,600);
   
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double mmu = 0.1057;
   double Ecm = beame;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   double pmu = sqrt(pow(beame,2)-pow(mmu,2));
   std::cout<<"Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;
   double beamlow=peakvalue-0.1;
   double beamup=peakvalue+0.2;
   if (Ecm<2.6) {
     beamlow=peakvalue-0.05;
     beamup=peakvalue+0.15;
   }
   // save histogram of pmu
  //TH1D *hp = new TH1D("hp","hp",100,beamlow,beamup);
  //dataraw->Draw("x1>>hp");
  //char name[100];
  //sprintf(name,"hp_mcKK_%.4f",beame*2);
  //hp->Write(name);
  //delete hp;
  //return 0;
   //
   //
   // try to use roofit
   //
   RooRealVar x("p1","momentum",peakvalue,beamlow,beamup,"GeV");
   // signal
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.003,peakvalue+0.003);
   RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.013);
   //RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.014);
   RooRealVar alpha1("alpha1","#alpha",1.5,1.0,5.0);
   RooRealVar nnn1("n1","n",2,1,10);
   RooCBShape cbshape("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);
   RooRealVar sigma2("sigma2","width of gaussian",0.007,0.005,0.02);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma2);
   RooRealVar frac1("frac1","frac1",0.9,0.6,1.0);

   RooAddPdf sig("sig","signal",RooArgList(cbshape,gaus),RooArgList(frac1));
   
   // background, mainly di-mu
   RooRealVar meanb("meanb","mean of gaussian",pmu,pmu-0.002,pmu+0.003);
   //RooRealVar sigmab("sigmab","width of gaussian",sigma_ini-0.001, sigma_ini-0.001,sigma_ini+0.001);
   RooRealVar sigmab("sigmab","width of gaussian",sigma_ini, sigma_ini-0.001,sigma_ini+0.001);
   //RooRealVar sigmab("sigmab","width of gaussian",0.005,0.002,0.02);
 //  RooGaussian gausb("gausb","gauss(x,m,s)",x,meanb,sigmab);
   //RooRealVar alphab("alphab","#alpha",0.8,0.5,1.5);
   RooRealVar alphab("alphab","#alpha",1.2,1.0,1.5);
   //RooRealVar alphab("alphab","#alpha",0.5,0.1,1.5);
   //RooRealVar nnnb("nnnb","n",0.5,0.4,2.0);
   RooRealVar nnnb("nnnb","n",0.5,0.1,2.0);
   //RooRealVar nnnb("nnnb","n",beame);
   RooCBShape cbshapeb("cbshapeb","crystal ball",x,meanb,sigmab,alphab,nnnb);
   //RooRealVar meane("meane","mean of gaussian",beame,beame-0.005,beame+0.003);
   RooRealVar sigmae("sigmae","width of gaussian",0.01,0.005,0.08);
   RooGaussian gause("gause","gauss(x,m,s)",x,meanb,sigmae);
   RooRealVar frac2("frac2","frac2",0.9,0.6,1.0);
   RooAddPdf bck("bck","signal",RooArgList(cbshapeb,gause),RooArgList(frac2));
   
   RooRealVar signal("signal"," ",1000,0,1000000000);//event number
   RooRealVar background("background"," ",20000,0,100000000);
 
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title(""));
   c1->SetMargin(0.13,0.1,0.13,0.1);
   TGaxis::SetMaxDigits(3);
   xframe->SetTitle("");
   xframe->GetXaxis()->SetTitle("p (GeV/c)");
   xframe->GetXaxis()->SetNdivisions(505);
   xframe->GetXaxis()->SetTitleSize(0.05);
   xframe->GetXaxis()->SetLabelSize(0.05);
   xframe->GetXaxis()->SetTitleOffset(1.1);
   xframe->GetYaxis()->SetNdivisions(505);
   xframe->GetYaxis()->SetTitleSize(0.05);
   xframe->GetYaxis()->SetTitleOffset(1.25);
   xframe->GetYaxis()->SetLabelSize(0.05);
   xframe->GetYaxis()->SetNoExponent(0);
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(sig,bck),RooArgList(signal,background));
   Npar = 10;
   //x.setRange("sigragi",peakvalue-0.1,beame+0.012);
   //sum->fitTo(*dataset);
   sum->fitTo(*dataset);
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bck),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt;
   if (Ecm<2.6) pt = new TPaveText(0.5,0.55,0.85,0.90,"BRNDC");
   else pt = new TPaveText(0.15,0.55,0.50,0.90,"BRNDC");
   pt->SetBorderSize(0);
   //pt->SetFillStyle(4000);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.045);
   sprintf(tmpchr,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"signal = %.1f #pm %.1f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
 //sprintf(tmpchr,"p_spectrum_%s",namesfx);
 //c1->SetName(tmpchr);
 //c1->Write();
   sprintf(tmpchr,"p_spectrum_%s.pdf",namesfx);
   c1->Print(tmpchr);

   ofstream ofsig("signal.txt",ios::app);
   ofsig << Ecm <<"\t"<< signal.getVal() <<"\t"<<signal.getError()<<endl;
 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   ofstream fitpar("fitpar",std::ios::app);
 //fitpar<<" ene = "<< beame*2 <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
 //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
 //fitpar<<"\t sigNo = " << signal.getVal() << "\t sigNoE = " << signal.getError();
 //fitpar<<"\t bckNo = " << background.getVal() << "\t bckNoE = " << background.getError();
 //fitpar<< "\t chi2 = " << xframe->chiSquare(Npar);
 //fitpar<< std::endl;
   fitpar<< beame*2 <<"\t ";
   
   fitpar<< mean.getVal() <<"\t ";
   fitpar<< sigma.getVal() <<"\t ";
   fitpar<< sigma2.getVal() <<"\t ";
   fitpar<< alpha1.getVal()<<"\t ";
   fitpar<< nnn1.getVal()<<"\t ";
   fitpar<< frac1.getVal()<<"\t ";
   
   fitpar<< meanb.getVal()<<"\t ";
   fitpar<< sigmab.getVal()<<"\t ";
   fitpar<< sigmae.getVal()<<"\t ";
   fitpar<< alphab.getVal()<<"\t ";
   fitpar<< nnnb.getVal()<<"\t ";
   fitpar<< frac2.getVal()<<"\t ";
   
   fitpar<< endl;
   //fitpar<<"\t e mean = " << meane.getVal() << "\t e sigma = " << sigmae.getVal()<<std::endl;
   
   double intel3 = mean.getVal()-3*sigma.getVal();
   double inteu3 = mean.getVal()+3*sigma.getVal();
   double intel5 = mean.getVal()-5*sigma.getVal();
   double inteu5 = mean.getVal()+5*sigma.getVal();
   //RooRealVar xint("xint","xint",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame+0.01);
   x.setRange("sigrag3",intel3, inteu3);
   x.setRange("sigrag5",intel5, inteu5);
   RooAbsReal* intsigi = sig.createIntegral( x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intbcki = bck.createIntegral(x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intsig3 = sig.createIntegral( x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intbck3 = bck.createIntegral(x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intsig5 = sig.createIntegral( x,  NormSet(x),  Range("sigrag5"));
   RooAbsReal* intbck5 = bck.createIntegral(x,  NormSet(x),  Range("sigrag5"));
   double signalN3 =   signal.getVal()*intsig3->getVal();
   double backN3 = background.getVal()*intbck3->getVal();
   double signalN5 =   signal.getVal()*intsig5->getVal();
   double backN5 = background.getVal()*intbck5->getVal();
 //double signalN3nom =   signal.getVal()*intsig3->getVal()/intsigi->getVal();
 //double backN3nom = background.getVal()*intbck3->getVal()/intbcki->getVal();
 //double signalN5nom =   signal.getVal()*intsig5->getVal()/intsigi->getVal();
 //double backN5nom = background.getVal()*intbck5->getVal()/intbcki->getVal();
   ofstream angSNR("angSNR.dat",std::ios::app);
   angSNR << namesfx <<"\t"<<signalN3 << "\t" << backN3 <<"\t"<<signalN5 << "\t" << backN5 << std::endl;
 //epSNR << "\t"<<signalN3nom << "\t" << backN3nom <<"\t"<<signalN5nom << "\t" << backN5nom<< std::endl;
   std::cout<< "Total signal int is "<< intsigi->getVal() <<" Total bck int is "<< intbcki->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig3->getVal() <<" Total bck int is "<< intbck3->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig5->getVal() <<" Total bck int is "<< intbck5->getVal()<<std::endl;

   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   err = signal.getError();
   return signal.getVal();
}

/*
double FitSpectrum3(TTree *dataraw, double beame, double &err,const char* namesfx)
{ 
   TF1 fsigma("fsigma","-0.00705519+0.00600716*x",2.0,3.1);
   double sigma_ini = fsigma.Eval(beame);
   TCanvas *c1 = new TCanvas("","",800,600);
   
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double mmu = 0.1057;
   double Ecm = beame;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   double pmu = sqrt(pow(beame,2)-pow(mmu,2));
   std::cout<<"Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;

   double beamlow=peakvalue-0.1;
   double beamup=peakvalue+0.2;
   if (Ecm<2.6) {
     beamlow=peakvalue-0.05;
     beamup=peakvalue+0.15;
   }
 
   //
   // try to use roofit
   //
   RooRealVar x("p1","momentum",peakvalue,beamlow,beamup,"GeV/c");
   // signal
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.003,peakvalue+0.003);
   RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.013);
   //RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.014);
   RooRealVar alpha1("alpha1","#alpha",1.5,1.0,5.0);
   RooRealVar nnn1("n1","n",2,1,10);
   RooCBShape cbshape("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);
   RooRealVar sigma2("sigma2","width of gaussian",0.007,0.005,0.02);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma2);
   RooRealVar frac1("frac1","frac1",0.9,0.1,1.0);

   RooAddPdf sig("sig","signal",RooArgList(cbshape,gaus),RooArgList(frac1));
   
   // background, mainly di-mu
   RooRealVar meanb("meanb","mean of gaussian",pmu,pmu-0.002,pmu+0.003);
   //RooRealVar sigmab("sigmab","width of gaussian",sigma_ini-0.001, sigma_ini-0.001,sigma_ini+0.001);
   RooRealVar sigmab("sigmab","width of gaussian",sigma_ini, sigma_ini-0.001,sigma_ini+0.001);
   //RooRealVar sigmab("sigmab","width of gaussian",0.005,0.002,0.02);
 //  RooGaussian gausb("gausb","gauss(x,m,s)",x,meanb,sigmab);
   RooRealVar alphab("alphab","#alpha",0.5,0.5,5);
   //RooRealVar alphab("alphab","#alpha",1.2,1.0,1.5);
   //RooRealVar alphab("alphab","#alpha",0.5,0.1,1.5);
   //RooRealVar nnnb("nnnb","n",0.5,0.4,2.0);
   RooRealVar nnnb("nnnb","n",0.6,0.1,1.5);
   //RooRealVar nnnb("nnnb","n",beame);
   RooCBShape cbshapeb("cbshapeb","crystal ball",x,meanb,sigmab,alphab,nnnb);
   //RooRealVar meane("meane","mean of gaussian",beame,beame-0.005,beame+0.003);
   RooRealVar sigmae("sigmae","width of gaussian",0.01,0.005,0.08);
   RooGaussian gause("gause","gauss(x,m,s)",x,meanb,sigmae);
   RooRealVar frac2("frac2","frac2",0.9,0.6,1.0);
   RooAddPdf bck("bck","signal",RooArgList(cbshapeb,gause),RooArgList(frac2));
   
   RooRealVar signal("signal"," ",1000,0,1000000000);//event number
   RooRealVar background("background"," ",20000,0,100000000);
 
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title(""));
   c1->SetMargin(0.13,0.1,0.13,0.1);
   TGaxis::SetMaxDigits(3);
   xframe->SetTitle("");
   xframe->GetXaxis()->SetTitle("p (GeV/c)");
   xframe->GetXaxis()->SetNdivisions(505);
   xframe->GetXaxis()->SetTitleSize(0.05);
   xframe->GetXaxis()->SetLabelSize(0.05);
   xframe->GetXaxis()->SetTitleOffset(1.1);
   xframe->GetYaxis()->SetNdivisions(505);
   xframe->GetYaxis()->SetTitleSize(0.05);
   xframe->GetYaxis()->SetTitleOffset(1.25);
   xframe->GetYaxis()->SetLabelSize(0.05);
   xframe->GetYaxis()->SetNoExponent(0);
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(sig,bck),RooArgList(signal,background));
   Npar = 10;
   //x.setRange("sigragi",peakvalue-0.1,beame+0.012);
   //sum->fitTo(*dataset);
   sum->fitTo(*dataset);
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bck),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt;
   if (Ecm<2.65) pt = new TPaveText(0.5,0.55,0.85,0.90,"BRNDC");
   else pt = new TPaveText(0.15,0.55,0.50,0.90,"BRNDC");
   pt->SetBorderSize(0);
   //pt->SetFillStyle(4000);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.045);
   sprintf(tmpchr,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"signal = %.1f #pm %.1f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
   //sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,xframe->chiSquare(Npar));
   sprintf(tmpchr,"#chi^{2} = %5.3f",xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
   sprintf(tmpchr,"p_spectrum_%s",namesfx);
   c1->SetName(tmpchr);
   c1->Write();
   sprintf(tmpchr,"p_spectrum_%s.pdf",namesfx);
   c1->Print(tmpchr);

 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   ofstream fitpar("fitpar",std::ios::app);
 //fitpar<<" ene = "<< beame*2 <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
 //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
 //fitpar<<"\t sigNo = " << signal.getVal() << "\t sigNoE = " << signal.getError();
 //fitpar<<"\t bckNo = " << background.getVal() << "\t bckNoE = " << background.getError();
 //fitpar<< "\t chi2 = " << xframe->chiSquare(Npar);
 //fitpar<< std::endl;
   fitpar<< beame*2 <<"\t ";
   
   fitpar<< mean.getVal() <<"\t ";
   fitpar<< sigma.getVal() <<"\t ";
   fitpar<< sigma2.getVal() <<"\t ";
   fitpar<< alpha1.getVal()<<"\t ";
   fitpar<< nnn1.getVal()<<"\t ";
   fitpar<< frac1.getVal()<<"\t ";
   
   fitpar<< meanb.getVal()<<"\t ";
   fitpar<< sigmab.getVal()<<"\t ";
   fitpar<< sigmae.getVal()<<"\t ";
   fitpar<< alphab.getVal()<<"\t ";
   fitpar<< nnnb.getVal()<<"\t ";
   fitpar<< frac2.getVal()<<"\t ";
   
   fitpar<< endl;
   //fitpar<<"\t e mean = " << meane.getVal() << "\t e sigma = " << sigmae.getVal()<<std::endl;
   
   double intel3 = mean.getVal()-3*sigma.getVal();
   double inteu3 = mean.getVal()+3*sigma.getVal();
   double intel5 = mean.getVal()-5*sigma.getVal();
   double inteu5 = mean.getVal()+5*sigma.getVal();
   //RooRealVar xint("xint","xint",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame+0.01);
   x.setRange("sigrag3",intel3, inteu3);
   x.setRange("sigrag5",intel5, inteu5);
   RooAbsReal* intsigi = sig.createIntegral( x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intbcki = bck.createIntegral(x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intsig3 = sig.createIntegral( x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intbck3 = bck.createIntegral(x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intsig5 = sig.createIntegral( x,  NormSet(x),  Range("sigrag5"));
   RooAbsReal* intbck5 = bck.createIntegral(x,  NormSet(x),  Range("sigrag5"));
   double signalN3 =   signal.getVal()*intsig3->getVal();
   double backN3 = background.getVal()*intbck3->getVal();
   double signalN5 =   signal.getVal()*intsig5->getVal();
   double backN5 = background.getVal()*intbck5->getVal();
 //double signalN3nom =   signal.getVal()*intsig3->getVal()/intsigi->getVal();
 //double backN3nom = background.getVal()*intbck3->getVal()/intbcki->getVal();
 //double signalN5nom =   signal.getVal()*intsig5->getVal()/intsigi->getVal();
 //double backN5nom = background.getVal()*intbck5->getVal()/intbcki->getVal();
   ofstream angSNR("angSNR.dat",std::ios::app);
   angSNR << namesfx <<"\t"<<signalN3 << "\t" << backN3 <<"\t"<<signalN5 << "\t" << backN5 << std::endl;
 //epSNR << "\t"<<signalN3nom << "\t" << backN3nom <<"\t"<<signalN5nom << "\t" << backN5nom<< std::endl;
   std::cout<< "Total signal int is "<< intsigi->getVal() <<" Total bck int is "<< intbcki->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig3->getVal() <<" Total bck int is "<< intbck3->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig5->getVal() <<" Total bck int is "<< intbck5->getVal()<<std::endl;

   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   err = signal.getError();
   return signal.getVal();
}
*/

// modify parameter limits , frac
/*
double FitSpectrum3(TTree *dataraw, double beame, double &err, const char* namesfx)
{
   TF1 fsigma("fsigma","-0.00705519+0.00600716*x",2.0,3.1);
   double sigma_ini = fsigma.Eval(beame);
   TCanvas *c1 = new TCanvas("","",800,600);
   
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double mmu = 0.1057;
   double Ecm = beame;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   double pmu = sqrt(pow(beame,2)-pow(mmu,2));
   std::cout<<"Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;

   double beamlow=peakvalue-0.1;
   double beamup=peakvalue+0.2;
   if (Ecm<2.6) {
     beamlow=peakvalue-0.05;
     beamup=peakvalue+0.15;
   }
 
   //
   // try to use roofit
   //
   RooRealVar x("p1","momentum",peakvalue,beamlow,beamup,"GeV/c");
   // signal
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.003,peakvalue+0.003);
   RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.013);
   //RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.014);
   RooRealVar alpha1("alpha1","#alpha",1.5,1.0,5.0);
   RooRealVar nnn1("n1","n",2,1,10);
   RooCBShape cbshape("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);
   RooRealVar sigma2("sigma2","width of gaussian",0.007,0.005,0.02);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma2);
   RooRealVar frac1("frac1","frac1",0.6,0.4,0.8);

   RooAddPdf sig("sig","signal",RooArgList(cbshape,gaus),RooArgList(frac1));
   
   // background, mainly di-mu
   RooRealVar meanb("meanb","mean of gaussian",pmu,pmu-0.002,pmu+0.003);
   //RooRealVar sigmab("sigmab","width of gaussian",sigma_ini-0.001, sigma_ini-0.001,sigma_ini+0.001);
   RooRealVar sigmab("sigmab","width of gaussian",sigma_ini, sigma_ini-0.001,sigma_ini+0.001);
   //RooRealVar sigmab("sigmab","width of gaussian",0.005,0.002,0.02);
 //  RooGaussian gausb("gausb","gauss(x,m,s)",x,meanb,sigmab);
   RooRealVar alphab("alphab","#alpha",0.5,0.5,5);
   //RooRealVar alphab("alphab","#alpha",1.2,1.0,1.5);
   //RooRealVar alphab("alphab","#alpha",0.5,0.1,1.5);
   //RooRealVar nnnb("nnnb","n",0.5,0.4,2.0);
   RooRealVar nnnb("nnnb","n",0.6,0.1,1.5);
   //RooRealVar nnnb("nnnb","n",beame);
   RooCBShape cbshapeb("cbshapeb","crystal ball",x,meanb,sigmab,alphab,nnnb);
   //RooRealVar meane("meane","mean of gaussian",beame,beame-0.005,beame+0.003);
   RooRealVar sigmae("sigmae","width of gaussian",0.01,0.005,0.08);
   RooGaussian gause("gause","gauss(x,m,s)",x,meanb,sigmae);
   RooRealVar frac2("frac2","frac2",0.8,0.7,0.9);
   RooAddPdf bck("bck","signal",RooArgList(cbshapeb,gause),RooArgList(frac2));
   
   RooRealVar signal("signal"," ",1000,0,1000000000);//event number
   RooRealVar background("background"," ",20000,0,100000000);
 
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title(""));
   c1->SetMargin(0.13,0.1,0.13,0.1);
   TGaxis::SetMaxDigits(3);
   xframe->SetTitle("");
   xframe->GetXaxis()->SetTitle("p (GeV/c)");
   xframe->GetXaxis()->SetNdivisions(505);
   xframe->GetXaxis()->SetTitleSize(0.05);
   xframe->GetXaxis()->SetLabelSize(0.05);
   xframe->GetXaxis()->SetTitleOffset(1.1);
   xframe->GetYaxis()->SetNdivisions(505);
   xframe->GetYaxis()->SetTitleSize(0.05);
   xframe->GetYaxis()->SetTitleOffset(1.25);
   xframe->GetYaxis()->SetLabelSize(0.05);
   xframe->GetYaxis()->SetNoExponent(0);
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(sig,bck),RooArgList(signal,background));
   Npar = 10;
   //x.setRange("sigragi",peakvalue-0.1,beame+0.012);
   sum->fitTo(*dataset);
   
   
//vector<double> v_Nsig; 
//vector<double> v_Esig; 
// 
//int max_time = 10;
//int i_time = 0;
//do
//{
//  sum->fitTo(*dataset);
//  v_Nsig.push_back(signal.getVal());
//  v_Esig.push_back(signal.getError());
//  //sum.fitTo(data);
//  i_time++;
//  cout<<i_time<<"th sigma "<< sigma.getVal()<<endl;
//  if (i_time==max_time) break;
//  else {
//    time_t now_time = time(NULL);
//    gRandom->SetSeed(now_time);
//    mean.setVal(gRandom->Uniform(peakvalue-0.003, peakvalue+0.003));
//    sigma.setVal(gRandom->Uniform(0.004, 0.013));
//    alpha1.setVal(gRandom->Uniform(1.0, 5.0));
//    nnn1.setVal(gRandom->Uniform(1.0, 10.0));
//    sigma2.setVal(gRandom->Uniform(0.005, 0.02));
//    frac1.setVal(gRandom->Uniform(0.4, 0.8));
//    
//    meanb.setVal(gRandom->Uniform(pmu-0.002, pmu+0.003));
//    sigmab.setVal(gRandom->Uniform(sigma_ini-0.001, sigma_ini+0.001));
//    alphab.setVal(gRandom->Uniform(0.5, 5));
//    nnnb.setVal(gRandom->Uniform(0.1, 1.5));
//    sigmae.setVal(gRandom->Uniform(0.005, 0.03));
//    frac2.setVal(gRandom->Uniform(0.7, 0.9));
//    
//    signal.setVal(gRandom->Uniform(100, 2000));
//    background.setVal(gRandom->Uniform(100, 2000));

//  }
//} while (i_time<max_time);
// 
//double Nsig_mean=0, Esig_mean=0;
//for (int i=0; i<v_Nsig.size(); i++){
//  Nsig_mean += v_Nsig.at(i);
//  Esig_mean += v_Esig.at(i);
//} 
//Nsig_mean = Nsig_mean/v_Nsig.size();
//Esig_mean = Esig_mean/v_Nsig.size();

//TH1D *hNsig = new TH1D("hNsig","Nsig",100,Nsig_mean-5*Esig_mean, Nsig_mean+5*Esig_mean);
//TH1D *hEsig = new TH1D("hEsig","Esig",100, 0, 2*Esig_mean);
//for (int i=0; i<v_Nsig.size(); i++){
//  hNsig->Fill(v_Nsig.at(i));
//} 
//hNsig->Draw("E");
//hNsig->Fit("gaus");
//Nsig_mean = hNsig->GetFunction("gaus")->GetParameter(1);
//sprintf(tmpchr,"hNsig_%s.pdf",namesfx);
//c1->Print(tmpchr);
//c1->Clear();
//delete hNsig;
//delete hEsig;

   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bck),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt;
   if (Ecm<2.65) pt = new TPaveText(0.5,0.55,0.85,0.90,"BRNDC");
   else pt = new TPaveText(0.15,0.55,0.50,0.90,"BRNDC");
   pt->SetBorderSize(0);
   //pt->SetFillStyle(4000);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.045);
   sprintf(tmpchr,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"signal = %.1f #pm %.1f",signal.getVal(),signal.getError());
   //sprintf(tmpchr,"signal = %.1f #pm %.1f", Nsig_mean, Esig_mean);
   pt->AddText(tmpchr);
   //sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,xframe->chiSquare(Npar));
   sprintf(tmpchr,"#chi^{2} = %5.3f",xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
 //sprintf(tmpchr,"p_spectrum_%s",namesfx);
 //c1->SetName(tmpchr);
 //c1->Write();
   sprintf(tmpchr,"p_spectrum_%s.pdf",namesfx);
   c1->Print(tmpchr);

   err = signal.getError();
   //err = Esig_mean;

 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   ofstream fitpar("fitpar",std::ios::app);
 //fitpar<<" ene = "<< beame*2 <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
 //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
 //fitpar<<"\t sigNo = " << signal.getVal() << "\t sigNoE = " << signal.getError();
 //fitpar<<"\t bckNo = " << background.getVal() << "\t bckNoE = " << background.getError();
 //fitpar<< "\t chi2 = " << xframe->chiSquare(Npar);
 //fitpar<< std::endl;
   fitpar<< beame*2 <<"\t ";
   
   fitpar<< mean.getVal() <<"\t ";
   fitpar<< sigma.getVal() <<"\t ";
   fitpar<< sigma2.getVal() <<"\t ";
   fitpar<< alpha1.getVal()<<"\t ";
   fitpar<< nnn1.getVal()<<"\t ";
   fitpar<< frac1.getVal()<<"\t ";
   
   fitpar<< meanb.getVal()<<"\t ";
   fitpar<< sigmab.getVal()<<"\t ";
   fitpar<< sigmae.getVal()<<"\t ";
   fitpar<< alphab.getVal()<<"\t ";
   fitpar<< nnnb.getVal()<<"\t ";
   fitpar<< frac2.getVal()<<"\t ";
   
   fitpar<< endl;
   //fitpar<<"\t e mean = " << meane.getVal() << "\t e sigma = " << sigmae.getVal()<<std::endl;
   
   double intel3 = mean.getVal()-3*sigma.getVal();
   double inteu3 = mean.getVal()+3*sigma.getVal();
   double intel5 = mean.getVal()-5*sigma.getVal();
   double inteu5 = mean.getVal()+5*sigma.getVal();
   //RooRealVar xint("xint","xint",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame+0.01);
   x.setRange("sigrag3",intel3, inteu3);
   x.setRange("sigrag5",intel5, inteu5);
   RooAbsReal* intsigi = sig.createIntegral( x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intbcki = bck.createIntegral(x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intsig3 = sig.createIntegral( x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intbck3 = bck.createIntegral(x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intsig5 = sig.createIntegral( x,  NormSet(x),  Range("sigrag5"));
   RooAbsReal* intbck5 = bck.createIntegral(x,  NormSet(x),  Range("sigrag5"));
   double signalN3 =   signal.getVal()*intsig3->getVal();
   double backN3 = background.getVal()*intbck3->getVal();
   double signalN5 =   signal.getVal()*intsig5->getVal();
   double backN5 = background.getVal()*intbck5->getVal();
 //double signalN3nom =   signal.getVal()*intsig3->getVal()/intsigi->getVal();
 //double backN3nom = background.getVal()*intbck3->getVal()/intbcki->getVal();
 //double signalN5nom =   signal.getVal()*intsig5->getVal()/intsigi->getVal();
 //double backN5nom = background.getVal()*intbck5->getVal()/intbcki->getVal();
   ofstream angSNR("angSNR.dat",std::ios::app);
   angSNR << namesfx <<"\t"<<signalN3 << "\t" << backN3 <<"\t"<<signalN5 << "\t" << backN5 << std::endl;
 //epSNR << "\t"<<signalN3nom << "\t" << backN3nom <<"\t"<<signalN5nom << "\t" << backN5nom<< std::endl;
   std::cout<< "Total signal int is "<< intsigi->getVal() <<" Total bck int is "<< intbcki->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig3->getVal() <<" Total bck int is "<< intbck3->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig5->getVal() <<" Total bck int is "<< intbck5->getVal()<<std::endl;

   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete c1;
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
   //return Nsig_mean;
}
*/

double FitSpectrum_rg(TTree *dataraw, double beame, const char* namesfx)
{ 
   TF1 fsigma("fsigma","-0.00705519+0.00600716*x",2.0,3.1);
   double sigma_ini = fsigma.Eval(beame);
   TCanvas *c1 = new TCanvas("","",800,600);
   
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double mmu = 0.1057;
   double Ecm = beame;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   double pmu = sqrt(pow(beame,2)-pow(mmu,2));
   std::cout<<"Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;
   double beamlow=peakvalue-0.1-0.02;
   double beamup=peakvalue+0.2+0.02;
   if (Ecm<2.6) {
     beamlow=peakvalue-0.05-0.02;
     beamup=peakvalue+0.15+0.02;
   }
 //double beamlow=peakvalue-0.1;
 //double beamup=peakvalue+0.2;
 //if (Ecm<2.6) {
 //  beamlow=peakvalue-0.05;
 //  beamup=peakvalue+0.15;
 //}
 
   // save histogram of pmu
  //TH1D *hp = new TH1D("hp","hp",100,beamlow,beamup);
  //dataraw->Draw("x1>>hp");
  //char name[100];
  //sprintf(name,"hp_mcKK_%.4f",beame*2);
  //hp->Write(name);
  //delete hp;
  //return 0;
   //
   //
   // try to use roofit
   //
   RooRealVar x("p1","momentum",peakvalue,beamlow,beamup,"GeV/c");
   // signal
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.003,peakvalue+0.003);
   RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.013);
   //RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.014);
   RooRealVar alpha1("alpha1","#alpha",1.5,1.0,5.0);
   RooRealVar nnn1("n1","n",2,1,10);
   RooCBShape cbshape("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);
   RooRealVar sigma2("sigma2","width of gaussian",0.007,0.005,0.02);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma2);
   RooRealVar frac1("frac1","frac1",0.6,0.4,0.8);

   RooAddPdf sig("sig","signal",RooArgList(cbshape,gaus),RooArgList(frac1));
   
   // background, mainly di-mu
   RooRealVar meanb("meanb","mean of gaussian",pmu,pmu-0.002,pmu+0.003);
   //RooRealVar sigmab("sigmab","width of gaussian",sigma_ini-0.001, sigma_ini-0.001,sigma_ini+0.001);
   RooRealVar sigmab("sigmab","width of gaussian",sigma_ini, sigma_ini-0.001,sigma_ini+0.001);
   //RooRealVar sigmab("sigmab","width of gaussian",0.005,0.002,0.02);
 //  RooGaussian gausb("gausb","gauss(x,m,s)",x,meanb,sigmab);
   RooRealVar alphab("alphab","#alpha",0.5,0.5,5);
   //RooRealVar alphab("alphab","#alpha",1.2,1.0,1.5);
   //RooRealVar alphab("alphab","#alpha",0.5,0.1,1.5);
   //RooRealVar nnnb("nnnb","n",0.5,0.4,2.0);
   RooRealVar nnnb("nnnb","n",0.6,0.1,1.5);
   //RooRealVar nnnb("nnnb","n",beame);
   RooCBShape cbshapeb("cbshapeb","crystal ball",x,meanb,sigmab,alphab,nnnb);
   //RooRealVar meane("meane","mean of gaussian",beame,beame-0.005,beame+0.003);
   RooRealVar sigmae("sigmae","width of gaussian",0.01,0.005,0.08);
   RooGaussian gause("gause","gauss(x,m,s)",x,meanb,sigmae);
   RooRealVar frac2("frac2","frac2",0.8,0.7,0.9);
   RooAddPdf bck("bck","signal",RooArgList(cbshapeb,gause),RooArgList(frac2));
   
   RooRealVar signal("signal"," ",1000,0,1000000000);//event number
   RooRealVar background("background"," ",20000,0,100000000);
 
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title(""));
   c1->SetMargin(0.13,0.1,0.13,0.1);
   TGaxis::SetMaxDigits(3);
   xframe->SetTitle("");
   xframe->GetXaxis()->SetTitle("p (GeV/c)");
   xframe->GetXaxis()->SetNdivisions(505);
   xframe->GetXaxis()->SetTitleSize(0.05);
   xframe->GetXaxis()->SetLabelSize(0.05);
   xframe->GetXaxis()->SetTitleOffset(1.1);
   xframe->GetYaxis()->SetNdivisions(505);
   xframe->GetYaxis()->SetTitleSize(0.05);
   xframe->GetYaxis()->SetTitleOffset(1.25);
   xframe->GetYaxis()->SetLabelSize(0.05);
   xframe->GetYaxis()->SetNoExponent(0);
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(sig,bck),RooArgList(signal,background));
   Npar = 10;
   //x.setRange("sigragi",peakvalue-0.1,beame+0.012);
   //sum->fitTo(*dataset);
   sum->fitTo(*dataset);
 
//vector<double> v_Nsig; 
//vector<double> v_Esig; 
// 
//int max_time = 10;
//int i_time = 0;
//do
//{
//  sum->fitTo(*dataset);
//  v_Nsig.push_back(signal.getVal());
//  v_Esig.push_back(signal.getError());
//  //sum.fitTo(data);
//  i_time++;
//  cout<<i_time<<"th sigma "<< sigma.getVal()<<endl;
//  if (i_time==max_time) break;
//  else {
//    time_t now_time = time(NULL);
//    gRandom->SetSeed(now_time);
//    mean.setVal(gRandom->Uniform(peakvalue-0.003, peakvalue+0.003));
//    sigma.setVal(gRandom->Uniform(0.004, 0.013));
//    alpha1.setVal(gRandom->Uniform(1.0, 5.0));
//    nnn1.setVal(gRandom->Uniform(1.0, 10.0));
//    sigma2.setVal(gRandom->Uniform(0.005, 0.02));
//    frac1.setVal(gRandom->Uniform(0.4, 0.8));
//    
//    meanb.setVal(gRandom->Uniform(pmu-0.002, pmu+0.003));
//    sigmab.setVal(gRandom->Uniform(sigma_ini-0.001, sigma_ini+0.001));
//    alphab.setVal(gRandom->Uniform(0.5, 5));
//    nnnb.setVal(gRandom->Uniform(0.1, 1.5));
//    sigmae.setVal(gRandom->Uniform(0.005, 0.03));
//    frac2.setVal(gRandom->Uniform(0.7, 0.9));
//    
//    signal.setVal(gRandom->Uniform(100, 2000));
//    background.setVal(gRandom->Uniform(100, 2000));

//  }
//} while (i_time<max_time);
// 
//double Nsig_mean=0, Esig_mean=0;
//for (int i=0; i<v_Nsig.size(); i++){
//  Nsig_mean += v_Nsig.at(i);
//  Esig_mean += v_Esig.at(i);
//} 
//Nsig_mean = Nsig_mean/v_Nsig.size();
//Esig_mean = Esig_mean/v_Nsig.size();

//TH1D *hNsig = new TH1D("hNsig","Nsig",100,Nsig_mean-5*Esig_mean, Nsig_mean+5*Esig_mean);
//TH1D *hEsig = new TH1D("hEsig","Esig",100, 0, 2*Esig_mean);
//for (int i=0; i<v_Nsig.size(); i++){
//  hNsig->Fill(v_Nsig.at(i));
//} 
//hNsig->Draw("E");
//hNsig->Fit("gaus");
//Nsig_mean = hNsig->GetFunction("gaus")->GetParameter(1);
//sprintf(tmpchr,"hNsig_%s_crg.pdf",namesfx);
//c1->Print(tmpchr);
//c1->Clear();
//delete hNsig;
//delete hEsig;
 
   
   
   
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bck),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt;
   if (Ecm<2.65) pt = new TPaveText(0.5,0.55,0.85,0.90,"BRNDC");
   else pt = new TPaveText(0.15,0.55,0.50,0.90,"BRNDC");
   pt->SetBorderSize(0);
   //pt->SetFillStyle(4000);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.045);
   sprintf(tmpchr,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"signal = %.1f #pm %.1f",signal.getVal(),signal.getError());
   //sprintf(tmpchr,"signal = %.1f #pm %.1f",Nsig_mean,Esig_mean);
   pt->AddText(tmpchr);
   //sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,xframe->chiSquare(Npar));
   sprintf(tmpchr,"#chi^{2} = %5.3f",xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
 //sprintf(tmpchr,"p_spectrum_%s",namesfx);
 //c1->SetName(tmpchr);
 //c1->Write();
   sprintf(tmpchr,"p_spectrum_%s_crg.pdf",namesfx);
   c1->Print(tmpchr);

 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   ofstream fitpar("fitpar",std::ios::app);
 //fitpar<<" ene = "<< beame*2 <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
 //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
 //fitpar<<"\t sigNo = " << signal.getVal() << "\t sigNoE = " << signal.getError();
 //fitpar<<"\t bckNo = " << background.getVal() << "\t bckNoE = " << background.getError();
 //fitpar<< "\t chi2 = " << xframe->chiSquare(Npar);
 //fitpar<< std::endl;
   fitpar<< beame*2 <<"\t ";
   
   fitpar<< mean.getVal() <<"\t ";
   fitpar<< sigma.getVal() <<"\t ";
   fitpar<< sigma2.getVal() <<"\t ";
   fitpar<< alpha1.getVal()<<"\t ";
   fitpar<< nnn1.getVal()<<"\t ";
   fitpar<< frac1.getVal()<<"\t ";
   
   fitpar<< meanb.getVal()<<"\t ";
   fitpar<< sigmab.getVal()<<"\t ";
   fitpar<< sigmae.getVal()<<"\t ";
   fitpar<< alphab.getVal()<<"\t ";
   fitpar<< nnnb.getVal()<<"\t ";
   fitpar<< frac2.getVal()<<"\t ";
   
   fitpar<< endl;
   //fitpar<<"\t e mean = " << meane.getVal() << "\t e sigma = " << sigmae.getVal()<<std::endl;
   
   double intel3 = mean.getVal()-3*sigma.getVal();
   double inteu3 = mean.getVal()+3*sigma.getVal();
   double intel5 = mean.getVal()-5*sigma.getVal();
   double inteu5 = mean.getVal()+5*sigma.getVal();
   //RooRealVar xint("xint","xint",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame+0.01);
   x.setRange("sigrag3",intel3, inteu3);
   x.setRange("sigrag5",intel5, inteu5);
   RooAbsReal* intsigi = sig.createIntegral( x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intbcki = bck.createIntegral(x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intsig3 = sig.createIntegral( x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intbck3 = bck.createIntegral(x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intsig5 = sig.createIntegral( x,  NormSet(x),  Range("sigrag5"));
   RooAbsReal* intbck5 = bck.createIntegral(x,  NormSet(x),  Range("sigrag5"));
   double signalN3 =   signal.getVal()*intsig3->getVal();
   double backN3 = background.getVal()*intbck3->getVal();
   double signalN5 =   signal.getVal()*intsig5->getVal();
   double backN5 = background.getVal()*intbck5->getVal();
 //double signalN3nom =   signal.getVal()*intsig3->getVal()/intsigi->getVal();
 //double backN3nom = background.getVal()*intbck3->getVal()/intbcki->getVal();
 //double signalN5nom =   signal.getVal()*intsig5->getVal()/intsigi->getVal();
 //double backN5nom = background.getVal()*intbck5->getVal()/intbcki->getVal();
   ofstream angSNR("angSNR.dat",std::ios::app);
   angSNR << namesfx <<"\t"<<signalN3 << "\t" << backN3 <<"\t"<<signalN5 << "\t" << backN5 << std::endl;
 //epSNR << "\t"<<signalN3nom << "\t" << backN3nom <<"\t"<<signalN5nom << "\t" << backN5nom<< std::endl;
   std::cout<< "Total signal int is "<< intsigi->getVal() <<" Total bck int is "<< intbcki->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig3->getVal() <<" Total bck int is "<< intbck3->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig5->getVal() <<" Total bck int is "<< intbck5->getVal()<<std::endl;

   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete c1;
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
   //return Nsig_mean;
}



/*
double FitSpectrum3(TTree *&dataraw, double beame, const char* namesfx)
{
   TF1 fsigma("fsigma","-0.00705519+0.00600716*x",2.0,3.1);
   double sigma_ini = fsigma.Eval(beame);
   TCanvas *c1 = new TCanvas("","",800,600);
   
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double mmu = 0.1057;
   double Ecm = beame;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   double pmu = sqrt(pow(beame,2)-pow(mmu,2));
   std::cout<<"Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;
   double beamlow=peakvalue-0.1;
   double beamup=peakvalue+0.2;
   if (Ecm<2.6) {
     beamlow=peakvalue-0.05;
     beamup=peakvalue+0.15;
   }
   // save histogram of pmu or pK
  //TH1D *hp = new TH1D("hp","hp",100,beamlow,beamup);
  //dataraw->Draw("x1>>hp");
  //char name[100];
  //sprintf(name,"hp_mckk_%.4f",beame*2);
  //hp->Write(name);
  //delete hp;
  //return 0;
   //
   //
   // try to use roofit
   //
   RooRealVar x("x1","momentum",peakvalue,beamlow,beamup,"GeV/c");
   // signal
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.003,peakvalue+0.003);
   RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.013);
   //RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.014);
   RooRealVar alpha1("alpha1","#alpha",2,1.0,5.0);
   RooRealVar nnn1("n1","n",2,1,10);
   RooCBShape cbshape("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);
   RooRealVar sigma2("sigma2","width of gaussian",0.007,0.005,0.02);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma2);
   RooRealVar frac1("frac1","frac1",0.8,0.1,1.0);

   RooAddPdf sig("sig","signal",RooArgList(cbshape,gaus),RooArgList(frac1));
   
   // background, mainly di-mu
   RooRealVar meanb("meanb","mean of gaussian",pmu,pmu-0.002,pmu+0.003);
   //RooRealVar sigmab("sigmab","width of gaussian",sigma_ini-0.001, sigma_ini-0.001,sigma_ini+0.001);
   RooRealVar sigmab("sigmab","width of gaussian",sigma_ini-0.0005, sigma_ini-0.002,sigma_ini+0.002);
   //RooRealVar sigmab("sigmab","width of gaussian",0.005,0.002,0.02);
 //  RooGaussian gausb("gausb","gauss(x,m,s)",x,meanb,sigmab);
   RooRealVar alphab("alphab","#alpha",1.0,0.1,5);
   //RooRealVar alphab("alphab","#alpha",1.2,1.0,1.5);
   //RooRealVar alphab("alphab","#alpha",0.5,0.1,1.5);
   //RooRealVar nnnb("nnnb","n",0.5,0.4,2.0);
   RooRealVar nnnb("nnnb","n",0.6,0.1,1.5);
   //RooRealVar nnnb("nnnb","n",beame);
   RooCBShape cbshapeb("cbshapeb","crystal ball",x,meanb,sigmab,alphab,nnnb);
   //RooRealVar meane("meane","mean of gaussian",beame,beame-0.005,beame+0.003);
   RooRealVar sigmae("sigmae","width of gaussian",0.01,0.005,0.08);
   RooGaussian gause("gause","gauss(x,m,s)",x,meanb,sigmae);
   RooRealVar frac2("frac2","frac2",0.9,0.6,1.0);
   RooAddPdf bck("bck","signal",RooArgList(cbshapeb,gause),RooArgList(frac2));
   
   RooRealVar signal("signal"," ",1000,0,1000000);//event number
   RooRealVar background("background"," ",2000,0,1000000);
 
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title(""));
   c1->SetMargin(0.13,0.1,0.13,0.1);
   TGaxis::SetMaxDigits(3);
   xframe->SetTitle("");
   xframe->GetXaxis()->SetTitle("p (GeV/c)");
   xframe->GetXaxis()->SetNdivisions(505);
   xframe->GetXaxis()->SetTitleSize(0.05);
   xframe->GetXaxis()->SetLabelSize(0.05);
   xframe->GetXaxis()->SetTitleOffset(1.1);
   xframe->GetYaxis()->SetNdivisions(505);
   xframe->GetYaxis()->SetTitleSize(0.05);
   xframe->GetYaxis()->SetTitleOffset(1.25);
   xframe->GetYaxis()->SetLabelSize(0.05);
   xframe->GetYaxis()->SetNoExponent(0);
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(sig,bck),RooArgList(signal,background));
   Npar = 10;
   //x.setRange("sigragi",peakvalue-0.1,beame+0.012);
   //sum->fitTo(*dataset);
   sum->fitTo(*dataset);
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bck),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt;
   if (Ecm<2.65) pt = new TPaveText(0.5,0.55,0.85,0.90,"BRNDC");
   else pt = new TPaveText(0.15,0.55,0.50,0.90,"BRNDC");
   pt->SetBorderSize(0);
   //pt->SetFillStyle(4000);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.045);
   sprintf(tmpchr,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"signal = %.1f #pm %.1f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"bckgrd = %.1f #pm %.1f",background.getVal(),background.getError());
   pt->AddText(tmpchr);
   //sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,xframe->chiSquare(Npar));
   sprintf(tmpchr,"#chi^{2} = %5.3f",xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
   sprintf(tmpchr,"p_spectrum_%s",namesfx);
   c1->SetName(tmpchr);
   c1->Write();
   sprintf(tmpchr,"p_spectrum_%s.pdf",namesfx);
   c1->Print(tmpchr);

 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   ofstream fitpar("fitpar",std::ios::app);
 //fitpar<<" ene = "<< beame*2 <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
 //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
 //fitpar<<"\t sigNo = " << signal.getVal() << "\t sigNoE = " << signal.getError();
 //fitpar<<"\t bckNo = " << background.getVal() << "\t bckNoE = " << background.getError();
 //fitpar<< "\t chi2 = " << xframe->chiSquare(Npar);
 //fitpar<< std::endl;
   fitpar<< beame*2 <<"\t ";
   
   fitpar<< mean.getVal() <<"\t ";
   fitpar<< sigma.getVal() <<"\t ";
   fitpar<< sigma2.getVal() <<"\t ";
   fitpar<< alpha1.getVal()<<"\t ";
   fitpar<< nnn1.getVal()<<"\t ";
   fitpar<< frac1.getVal()<<"\t ";
   
   fitpar<< meanb.getVal()<<"\t ";
   fitpar<< sigmab.getVal()<<"\t ";
   fitpar<< sigmae.getVal()<<"\t ";
   fitpar<< alphab.getVal()<<"\t ";
   fitpar<< nnnb.getVal()<<"\t ";
   fitpar<< frac2.getVal()<<"\t ";
   
   fitpar<< endl;
   //fitpar<<"\t e mean = " << meane.getVal() << "\t e sigma = " << sigmae.getVal()<<std::endl;
   
   double intel3 = mean.getVal()-3*sigma.getVal();
   double inteu3 = mean.getVal()+3*sigma.getVal();
   double intel5 = mean.getVal()-5*sigma.getVal();
   double inteu5 = mean.getVal()+5*sigma.getVal();
   //RooRealVar xint("xint","xint",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame+0.01);
   x.setRange("sigrag3",intel3, inteu3);
   x.setRange("sigrag5",intel5, inteu5);
   RooAbsReal* intsigi = sig.createIntegral( x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intbcki = bck.createIntegral(x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intsig3 = sig.createIntegral( x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intbck3 = bck.createIntegral(x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intsig5 = sig.createIntegral( x,  NormSet(x),  Range("sigrag5"));
   RooAbsReal* intbck5 = bck.createIntegral(x,  NormSet(x),  Range("sigrag5"));
   double signalN3 =   signal.getVal()*intsig3->getVal();
   double backN3 = background.getVal()*intbck3->getVal();
   double signalN5 =   signal.getVal()*intsig5->getVal();
   double backN5 = background.getVal()*intbck5->getVal();
 //double signalN3nom =   signal.getVal()*intsig3->getVal()/intsigi->getVal();
 //double backN3nom = background.getVal()*intbck3->getVal()/intbcki->getVal();
 //double signalN5nom =   signal.getVal()*intsig5->getVal()/intsigi->getVal();
 //double backN5nom = background.getVal()*intbck5->getVal()/intbcki->getVal();
   ofstream angSNR("SNR.dat",std::ios::app);
   angSNR << namesfx <<"\t"<<signalN3 << "\t" << backN3 <<"\t"<<signalN5 << "\t" << backN5 << "\t"<<mean.getVal() << "\t"<<sigma.getVal()<< std::endl;
 //epSNR << "\t"<<signalN3nom << "\t" << backN3nom <<"\t"<<signalN5nom << "\t" << backN5nom<< std::endl;
   std::cout<< "Total signal int is "<< intsigi->getVal() <<" Total bck int is "<< intbcki->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig3->getVal() <<" Total bck int is "<< intbck3->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig5->getVal() <<" Total bck int is "<< intbck5->getVal()<<std::endl;

   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
}
*/

/*
// prove dimu back can discribe with the function 
double FitSpectrum3(TTree *&dataraw, double beame, const char* namesfx)
{
   TF1 fsigma("fsigma","-0.00705519+0.00600716*x",2.0,3.1);
   double sigma_ini = fsigma.Eval(beame);
   TCanvas *c1 = new TCanvas("","",800,600);
   
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double mmu = 0.1057;
   double Ecm = beame;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   double pmu = sqrt(pow(beame,2)-pow(mmu,2));
   std::cout<<"Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;
   double beamlow=peakvalue-0.1;
   double beamup=peakvalue+0.2;
   if (Ecm<2.6) {
     beamlow=peakvalue-0.05;
     beamup=peakvalue+0.15;
   }
   // save histogram of pmu or pK
  //TH1D *hp = new TH1D("hp","hp",100,beamlow,beamup);
  //dataraw->Draw("x1>>hp");
  //char name[100];
  //sprintf(name,"hp_mckk_%.4f",beame*2);
  //hp->Write(name);
  //delete hp;
  //return 0;
   //
   //
   // try to use roofit
   //
   RooRealVar x("x1","momentum",peakvalue,beamlow,beamup,"GeV/c");
   // signal
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.003,peakvalue+0.003);
   RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.013);
   //RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.014);
   RooRealVar alpha1("alpha1","#alpha",2,1.0,5.0);
   RooRealVar nnn1("n1","n",2,1,10);
   RooCBShape cbshape("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);
   RooRealVar sigma2("sigma2","width of gaussian",0.007,0.005,0.02);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma2);
   RooRealVar frac1("frac1","frac1",0.8,0.1,1.0);

   RooAddPdf sig("sig","signal",RooArgList(cbshape,gaus),RooArgList(frac1));
   
   // background, mainly di-mu
   RooRealVar meanb("meanb","mean of gaussian",pmu-0.001,pmu-0.003,pmu+0.003);
   //RooRealVar sigmab("sigmab","width of gaussian",sigma_ini-0.001, sigma_ini-0.001,sigma_ini+0.001);
   RooRealVar sigmab("sigmab","width of gaussian",sigma_ini-0.001, sigma_ini-0.002,sigma_ini+0.002);
   //RooRealVar sigmab("sigmab","width of gaussian",0.005,0.002,0.02);
 //  RooGaussian gausb("gausb","gauss(x,m,s)",x,meanb,sigmab);
   RooRealVar alphab("alphab","#alpha",1.0,0.1,5);
   //RooRealVar alphab("alphab","#alpha",1.2,1.0,1.5);
   //RooRealVar alphab("alphab","#alpha",0.5,0.1,1.5);
   //RooRealVar nnnb("nnnb","n",0.5,0.4,2.0);
   RooRealVar nnnb("nnnb","n",0.6,0.1,1.5);
   //RooRealVar nnnb("nnnb","n",beame);
   RooCBShape cbshapeb("cbshapeb","crystal ball",x,meanb,sigmab,alphab,nnnb);
   //RooRealVar meane("meane","mean of gaussian",beame,beame-0.005,beame+0.003);
   RooRealVar sigmae("sigmae","width of gaussian",0.01,0.005,0.08);
   RooGaussian gause("gause","gauss(x,m,s)",x,meanb,sigmae);
   RooRealVar frac2("frac2","frac2",0.9,0.6,1.0);
   RooAddPdf bck("bck","signal",RooArgList(cbshapeb,gause),RooArgList(frac2));
   
   RooRealVar signal("signal"," ",1000,0,1000000);//event number
   RooRealVar background("background"," ",2000,0,1000000);
 
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title(""));
   c1->SetMargin(0.13,0.1,0.13,0.1);
   TGaxis::SetMaxDigits(3);
   xframe->SetTitle("");
   xframe->GetXaxis()->SetTitle("p (GeV/c)");
   xframe->GetXaxis()->SetNdivisions(505);
   xframe->GetXaxis()->SetTitleSize(0.05);
   xframe->GetXaxis()->SetLabelSize(0.05);
   xframe->GetXaxis()->SetTitleOffset(1.1);
   xframe->GetYaxis()->SetNdivisions(505);
   xframe->GetYaxis()->SetTitleSize(0.05);
   xframe->GetYaxis()->SetTitleOffset(1.25);
   xframe->GetYaxis()->SetLabelSize(0.05);
   xframe->GetYaxis()->SetNoExponent(0);
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   //sum = new RooAddPdf("sum","sum",RooArgList(sig,bck),RooArgList(signal,background));
   sum = new RooAddPdf("sum","sum",RooArgList(bck),RooArgList(background));
   Npar = 10;
   //x.setRange("sigragi",peakvalue-0.1,beame+0.012);
   //sum->fitTo(*dataset);
   sum->fitTo(*dataset);
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bck),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt;
   if (Ecm<2.65) pt = new TPaveText(0.5,0.55,0.85,0.90,"BRNDC");
   else pt = new TPaveText(0.15,0.55,0.50,0.90,"BRNDC");
   pt->SetBorderSize(0);
   //pt->SetFillStyle(4000);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.045);
 //sprintf(tmpchr,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError());
 //pt->AddText(tmpchr);
 //sprintf(tmpchr,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError());
 //pt->AddText(tmpchr);
 //sprintf(tmpchr,"signal = %.1f #pm %.1f",signal.getVal(),signal.getError());
 //pt->AddText(tmpchr);
   //sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,xframe->chiSquare(Npar));
   sprintf(tmpchr,"#chi^{2} = %5.3f",xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
   sprintf(tmpchr,"p_spectrum_%s",namesfx);
   c1->SetName(tmpchr);
   c1->Write();
   sprintf(tmpchr,"p_spectrum_%s.pdf",namesfx);
   c1->Print(tmpchr);

 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   ofstream fitpar("fitpar",std::ios::app);
 //fitpar<<" ene = "<< beame*2 <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
 //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
 //fitpar<<"\t sigNo = " << signal.getVal() << "\t sigNoE = " << signal.getError();
 //fitpar<<"\t bckNo = " << background.getVal() << "\t bckNoE = " << background.getError();
 //fitpar<< "\t chi2 = " << xframe->chiSquare(Npar);
 //fitpar<< std::endl;
   fitpar<< beame*2 <<"\t ";
   
   fitpar<< mean.getVal() <<"\t ";
   fitpar<< sigma.getVal() <<"\t ";
   fitpar<< sigma2.getVal() <<"\t ";
   fitpar<< alpha1.getVal()<<"\t ";
   fitpar<< nnn1.getVal()<<"\t ";
   fitpar<< frac1.getVal()<<"\t ";
   
   fitpar<< meanb.getVal()<<"\t ";
   fitpar<< sigmab.getVal()<<"\t ";
   fitpar<< sigmae.getVal()<<"\t ";
   fitpar<< alphab.getVal()<<"\t ";
   fitpar<< nnnb.getVal()<<"\t ";
   fitpar<< frac2.getVal()<<"\t ";
   
   fitpar<< endl;
   //fitpar<<"\t e mean = " << meane.getVal() << "\t e sigma = " << sigmae.getVal()<<std::endl;
   
   double intel3 = mean.getVal()-3*sigma.getVal();
   double inteu3 = mean.getVal()+3*sigma.getVal();
   double intel5 = mean.getVal()-5*sigma.getVal();
   double inteu5 = mean.getVal()+5*sigma.getVal();
   //RooRealVar xint("xint","xint",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame+0.01);
   x.setRange("sigrag3",intel3, inteu3);
   x.setRange("sigrag5",intel5, inteu5);
   RooAbsReal* intsigi = sig.createIntegral( x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intbcki = bck.createIntegral(x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intsig3 = sig.createIntegral( x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intbck3 = bck.createIntegral(x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intsig5 = sig.createIntegral( x,  NormSet(x),  Range("sigrag5"));
   RooAbsReal* intbck5 = bck.createIntegral(x,  NormSet(x),  Range("sigrag5"));
   double signalN3 =   signal.getVal()*intsig3->getVal();
   double backN3 = background.getVal()*intbck3->getVal();
   double signalN5 =   signal.getVal()*intsig5->getVal();
   double backN5 = background.getVal()*intbck5->getVal();
 //double signalN3nom =   signal.getVal()*intsig3->getVal()/intsigi->getVal();
 //double backN3nom = background.getVal()*intbck3->getVal()/intbcki->getVal();
 //double signalN5nom =   signal.getVal()*intsig5->getVal()/intsigi->getVal();
 //double backN5nom = background.getVal()*intbck5->getVal()/intbcki->getVal();
   ofstream angSNR("angSNR.dat",std::ios::app);
   angSNR << namesfx <<"\t"<<signalN3 << "\t" << backN3 <<"\t"<<signalN5 << "\t" << backN5 << std::endl;
 //epSNR << "\t"<<signalN3nom << "\t" << backN3nom <<"\t"<<signalN5nom << "\t" << backN5nom<< std::endl;
   std::cout<< "Total signal int is "<< intsigi->getVal() <<" Total bck int is "<< intbcki->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig3->getVal() <<" Total bck int is "<< intbck3->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig5->getVal() <<" Total bck int is "<< intbck5->getVal()<<std::endl;

   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
}
*/


#include "RooHistPdf.h"
#include "RooFFTConvPdf.h"
double FitSpectrum4(TTree *&dataraw, double beame, double &err, const char* namesfx, TFile* fileout)
{
   TCanvas *c1 = new TCanvas("","",800,600);
   
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double mmu = 0.1057;
   double Ecm = beame;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   double pmu = sqrt(pow(beame,2)-pow(mmu,2));
   std::cout<<"Ecm is "<<Ecm<<", Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;
   //double beamlow=peakvalue-0.2;
   //double beamup=peakvalue+0.2;
   //double beamlow=peakvalue-0.10;
   //double beamup=peakvalue+0.10;
   double beamlow=peakvalue-0.1;
   double beamup=peakvalue+0.2;
   if (Ecm<2.6) {
     beamlow=peakvalue-0.05;
     beamup=peakvalue+0.15;
   }
   
   char hname[1000];
   TFile *file=new TFile("output/output_bck.root");
   sprintf(hname,"hp_mcKK2_%.4f",Ecm);
   TH1D *hk = (TH1D*)file->Get(hname);
   sprintf(hname,"hp_mcmumu_%.4f",Ecm);
   TH1D *hmu = (TH1D*)file->Get(hname);
   cout<<"aaaaaaaaa"<<endl;
   cout<<hk<<"\t"<<hmu<<endl;
 //sprintf(hname,"hkrebin_%.4f",Ecm);
 //TH1D *hkrebin = (TH1D*)hk->Rebin(2,hname);
 //sprintf(hname,"hmurebin_%.4f",Ecm);
 //TH1D *hmurebin = (TH1D*)hmu->Rebin(2,hname);
 //cout<<"hist dim "<< hkrebin->GetDimension()<<endl;
 //cout<<"aaaaaaaaa hkrebin pointer "<<hkrebin<<", hmurebin "<<hmurebin<<endl;
   //return 0;
     if (hk==0 || hmu==0) {cout<<"Error : spare pointer!!!!"<< endl; exit(-1);}
 //  hkrebin->Draw();
 //sprintf(hname,"hktmp_%.4f.pdf",Ecm);
 //  c1->Print(hname);
 //  hmurebin->Draw();
 //sprintf(hname,"hmutmp_%.4f.pdf",Ecm);
 //  c1->Print(hname);
   //
   if (fileout!=0) fileout->cd();
   //
   // try to use roofit
   //
   cout<<"aaaaaaaaa"<<endl;
   RooRealVar x("p1","momentum",peakvalue,beamlow,beamup,"GeV");
    // signal
 
   cout<<"aaaaaaaaa"<<endl;
   RooDataHist mckhist("mcKhist","mc K spec",x,hk);
   RooDataHist mcmuhist("mcmuhist","mc mu spec",x,hmu);
   cout<<"aaaaaaaaa"<<endl;
   RooHistPdf kpdf("kpdf","signal",x,mckhist,4); 
   RooHistPdf mupdf("mupdf","background",x,mcmuhist,4);
   cout<<"aaaaaaaaa"<<endl;

   RooRealVar mean_k("mean_k","mean of gaussian",0.0007,-0.005, 0.005);
   //RooRealVar mean_k("mean_k","mean of gaussian",0.001,-0.005, 0.005);
   RooRealVar sigma_k("sigma_k","width of gaussian",0.0024,0.0005,0.005);
   RooGaussian gaus_k("gaus_k","gauss(x,m,s)",x,mean_k,sigma_k);

   RooRealVar mean_u("mean","mean of gaussian",0.0007,-0.005, 0.005);
   RooRealVar sigma_u("sigma","width of gaussian",0.0024,0.0005,0.003);
   RooGaussian gaus_u("gaus_u","gauss(x,m,s)",x,mean_u,sigma_u);
   
   // signal
   RooFFTConvPdf sig("sig","signal",x,kpdf,gaus_k); 
   //RooHistPdf sig("sig","signal",x,mckhist,2); 
   // background, mainly di-mu
   RooFFTConvPdf bck("bck","background",x,mupdf,gaus_u); 
   //RooHistPdf bck("bck","background",x,mcmuhist,2);
   
   RooRealVar signal("signal"," ",1000,0,100000);//event number
   RooRealVar background("background"," ",200,0,100000);
 
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title("fit p"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(sig,bck),RooArgList(signal,background));
   //sum = new RooAddPdf("sum","sum",RooArgList(sig),RooArgList(signal));
   Npar = 6;
   //x.setRange("sigragi",peakvalue-0.1,beame+0.012);
   //sum->fitTo(*dataset,Range("sigragi"));
   //sum->fitTo(*dataset,Minimizer("Minuit","simplex"),Hesse(false));
   //sum->fitTo(*dataset,Minimizer("GSLSimAn"));
   //sum->fitTo(*dataset,Minimizer("GSLMultiMin"),Hesse(false) );
   //sum->chi2FitTo(*dataset);
   sum->fitTo(*dataset);
   err = signal.getError();
   dataset->plotOn(xframe,Binning(nBins));
   //sum->fitTo(mckhist);
   //mckhist.plotOn(xframe);
   sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bck),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt = new TPaveText(0.15,0.55,0.45,0.90,"BRNDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(4000);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.035);
   sprintf(tmpchr,"mean_{K} = %1.6f #pm %1.6f",mean_k.getVal(),mean_k.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma_{K} = %1.6f #pm %1.6f",sigma_k.getVal(),sigma_k.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"mean_{u} = %1.6f #pm %1.6f",mean_u.getVal(),mean_u.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma_{u} = %1.6f #pm %1.6f",sigma_u.getVal(),sigma_u.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"N_{sig} = %.2f #pm %.2f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"N_{bck} = %.2f #pm %.2f",background.getVal(),background.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
 //sprintf(tmpchr,"p_spectrum_%s",namesfx);
 //c1->SetName(tmpchr);
 //c1->Write();
   sprintf(tmpchr,"p_spectrum_%s.pdf",namesfx);
   c1->Print(tmpchr);
   
   ofstream ofsig("signal.txt",ios::app);
   ofsig << Ecm <<"\t"<< signal.getVal() <<"\t"<<signal.getError()<<endl;
   
   file->Close();
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
}

double FitSpectrum_mcshape(TTree *&dataraw, TH1D* hsig, TH1D* hbck, double beame, double &err, const char* namesfx, TFile* fileout, double *par_pool, int rw)
{
   TCanvas *c1 = new TCanvas("","",800,600);
   
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double mmu = 0.1057;
   double Ecm = beame;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   double pmu = sqrt(pow(beame,2)-pow(mmu,2));
   std::cout<<"Ecm is "<<Ecm<<", Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;
   double beamlow=peakvalue-0.1;
   double beamup=peakvalue+0.2;
   if (Ecm<2.6) {
     beamlow=peakvalue-0.05;
     beamup=peakvalue+0.15;
   }
   
   if (fileout!=0) fileout->cd();
   //
   // try to use roofit
   //
   cout<<"aaaaaaaaa"<<endl;
   RooRealVar x("p1","p",beamlow,beamup,"GeV/c");
    // signal
 
   cout<<"aaaaaaaaa"<<endl;
   RooDataHist mckhist("mcKhist","mc K spec",x,hsig);
   RooDataHist mcmuhist("mcmuhist","mc mu spec",x,hbck);
   cout<<"aaaaaaaaa"<<endl;
   RooHistPdf kpdf("kpdf","signal",x,mckhist,4); 
   RooHistPdf mupdf("mupdf","background",x,mcmuhist,4);
   cout<<"aaaaaaaaa"<<endl;

   RooRealVar mean_k("mean_k","mean of gaussian",0.0015,-0.005, 0.005);
   RooRealVar sigma_k("sigma_k","width of gaussian",0.002,0.0005,0.005);
   RooGaussian gaus_k("gaus_k","gauss(x,m,s)",x,mean_k,sigma_k);

   RooRealVar mean_u("mean_u","mean of gaussian",0.002,-0.005, 0.005);
   RooRealVar sigma_u("sigma_u","width of gaussian",0.0024,0.0005,0.005);
   RooGaussian gaus_u("gaus_u","gauss(x,m,s)",x,mean_u,sigma_u);
   
   if (par_pool!=0 && rw == 10){
     mean_k.setRange(par_pool[0], par_pool[0]);
     mean_k.setVal(par_pool[0]);
     sigma_k.setRange(par_pool[2], par_pool[2]);
     sigma_k.setVal(par_pool[2]);
     mean_u.setRange(par_pool[10], par_pool[10]);
     mean_u.setVal(par_pool[10]);
     sigma_u.setRange(par_pool[12], par_pool[12]);
     sigma_u.setVal(par_pool[12]);
   }

 //double sigma_f = -0.0279011 + 0.0234759*Ecm -0.0044394*Ecm*Ecm;
 //sigma_k.setRange(sigma_f, sigma_f);
 //sigma_k.setVal(sigma_f);
 //double sigmasp_f = 0.003;
 //sigma_u.setRange(sigmasp_f, sigmasp_f);
 //sigma_u.setVal(sigmasp_f);
   
   // signal
   RooFFTConvPdf sig("sig","signal",x,kpdf,gaus_k); 
   //RooHistPdf sig("sig","signal",x,mckhist,2); 
   // background, mainly di-mu
   RooFFTConvPdf bck("bck","background",x,mupdf,gaus_u); 
   //RooHistPdf bck("bck","background",x,mcmuhist,2);
   
   RooRealVar signal("signal"," ",1000,0,100000);//event number
   RooRealVar background("background"," ",200,0,100000);
 
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title(""));
   c1->SetMargin(0.15,0.1,0.15,0.1);
   TGaxis::SetMaxDigits(3);
   xframe->SetTitle("");
   xframe->GetXaxis()->SetTitle("p (GeV/c)");
   xframe->GetXaxis()->SetNdivisions(505);
   xframe->GetXaxis()->SetTitleSize(0.05);
   xframe->GetXaxis()->SetLabelSize(0.05);
   xframe->GetXaxis()->SetTitleOffset(1.1);
   xframe->GetYaxis()->SetNdivisions(505);
   xframe->GetYaxis()->SetTitleSize(0.05);
   xframe->GetYaxis()->SetTitleOffset(1.25);
   xframe->GetYaxis()->SetLabelSize(0.05);
   xframe->GetYaxis()->SetNoExponent(0);
 
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(sig,bck),RooArgList(signal,background));
   //sum = new RooAddPdf("sum","sum",RooArgList(sig),RooArgList(signal));
   Npar = 4;
   sum->fitTo(*dataset);
   err = signal.getError();
   dataset->plotOn(xframe,Binning(nBins));
   sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bck),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt;
   if (Ecm<2.6) pt = new TPaveText(0.55,0.65,0.88,0.85,"BRNDC");
   if (Ecm>2.6) pt = new TPaveText(0.18,0.65,0.45,0.85,"BRNDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(4000);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.035);
 //sprintf(tmpchr,"mean_{K} = %1.6f #pm %1.6f",mean_k.getVal(),mean_k.getError());
 //pt->AddText(tmpchr);
 //sprintf(tmpchr,"#sigma_{K} = %1.6f #pm %1.6f",sigma_k.getVal(),sigma_k.getError());
 //pt->AddText(tmpchr);
 //sprintf(tmpchr,"mean_{u} = %1.6f #pm %1.6f",mean_u.getVal(),mean_u.getError());
 //pt->AddText(tmpchr);
 //sprintf(tmpchr,"#sigma_{u} = %1.6f #pm %1.6f",sigma_u.getVal(),sigma_u.getError());
 //pt->AddText(tmpchr);
   sprintf(tmpchr,"N_{sig} = %.1f #pm %.1f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"N_{bck} = %.1f #pm %.1f",background.getVal(),background.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
   sprintf(tmpchr,"p_spectrum_%s.pdf",namesfx);
   c1->Print(tmpchr);
   
   ofstream ofpar("mcpar.txt",ios::app);
   ofpar << Ecm <<"\t"<< mean_k.getVal() <<"\t"<<mean_k.getError();
   ofpar        <<"\t"<<sigma_k.getVal() <<"\t"<<sigma_k.getError();
   ofpar        <<"\t"<< mean_u.getVal() <<"\t"<<mean_u.getError();
   ofpar        <<"\t"<<sigma_u.getVal() <<"\t"<<sigma_u.getError();
   ofpar << endl;
   
   ofstream ofsig("signal.txt",ios::app);
   ofsig << Ecm <<"\t"<< signal.getVal() <<"\t"<<signal.getError()<<endl;

   if (par_pool!=0 && rw == 1) {
     par_pool[0] = mean_k.getVal();
     par_pool[1] = mean_k.getError();
     par_pool[2] =sigma_k.getVal();
     par_pool[3] =sigma_k.getError();
     par_pool[10] = mean_u.getVal();
     par_pool[11] = mean_u.getError();
     par_pool[12] =sigma_u.getVal(); 
     par_pool[13] =sigma_u.getError(); 
   }

   delete xframe;
   delete dataset;
   delete sum;
   delete c1; 
   return signal.getVal();
}

double FitSpectrum_smcs(TTree *dataraw, TH1D* hsig, double beame, const char* namesfx, double *par_pool)
{
   TF1 fsigma("fsigma","-0.00705519+0.00600716*x",2.0,3.1);
   double sigma_ini = fsigma.Eval(beame);
   TCanvas *c1 = new TCanvas("","",800,600);
   
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double mmu = 0.1057;
   double Ecm = beame;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   double pmu = sqrt(pow(beame,2)-pow(mmu,2));
   std::cout<<"Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;
   double beamlow=peakvalue-0.1;
   double beamup=peakvalue+0.2;
   if (Ecm<2.6) {
     beamlow=peakvalue-0.05;
     beamup=peakvalue+0.15;
   }
   
   char hname[1000];
 //TFile *file=new TFile("output/output_bck.root");
 //sprintf(hname,"hp_mcKK2_%.4f",Ecm);
 //TH1D *hk = (TH1D*)file->Get(hname);
 //sprintf(hname,"hp_mcmumu_%.4f",Ecm);
 //TH1D *hmu = (TH1D*)file->Get(hname);
 //cout<<"aaaaaaaaa"<<endl;
 
   // try to use roofit
   //
   RooRealVar x("p1","momentum",peakvalue,beamlow,beamup,"GeV/c");
   // signal
 //RooDataHist mckhist("mcKhist","mc K spec",x,hk);
   RooDataHist mckhist("mcKhist","mc K spec",x,hsig);
   RooHistPdf kpdf("kpdf","signal",x,mckhist,4); 
   RooRealVar mean("mean","mean of gaussian",0.0015,-0.005, 0.005);
   RooRealVar sigma("sigma","width of gaussian",0.002,0.0005,0.003);
   if (par_pool!=0){
     mean.setRange(par_pool[0],par_pool[0]);
     mean.setVal(par_pool[0]);
     sigma.setRange(par_pool[2],par_pool[2]);
     sigma.setVal(par_pool[2]);
   }
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooFFTConvPdf sig("sig","signal",x,kpdf,gaus); 


   // background, mainly di-mu
   RooRealVar meanb("meanb","mean of gaussian",pmu,pmu-0.005,pmu+0.005);
   RooRealVar sigmab("sigmab","width of gaussian",sigma_ini, sigma_ini-0.001,sigma_ini+0.001);
   RooRealVar alphab("alphab","#alpha",0.5,0.5,5);
   RooRealVar nnnb("nnnb","n",0.6,0.1,1.5);
   RooCBShape cbshapeb("cbshapeb","crystal ball",x,meanb,sigmab,alphab,nnnb);
   RooRealVar sigmae("sigmae","width of gaussian",0.01,0.005,0.08);
   RooGaussian gause("gause","gauss(x,m,s)",x,meanb,sigmae);
   RooRealVar frac2("frac2","frac2",0.8,0.7,0.9);
   RooAddPdf bck("bck","signal",RooArgList(cbshapeb,gause),RooArgList(frac2));
   
   RooRealVar signal("signal"," ",1000,0,100000);//event number
   RooRealVar background("background"," ",2000,0,100000);
   
   // fix some parameter
 //double sigma_f = -0.0279011 + 0.0234759*Ecm -0.0044394*Ecm*Ecm;
 //sigma.setRange(sigma_f, sigma_f);
 //sigma.setVal(sigma_f);
 //double sigmab_f = -0.00286479 + 0.00447862*Ecm;
 //sigmab.setRange(sigmab_f, sigmab_f);
 //sigmab.setVal(sigmab_f);
 //double sigmae_f = -0.063683 + 0.0337589*Ecm;
 //sigmae.setRange(sigmae_f, sigmae_f);
 //sigmae.setVal(sigmae_f);
  
 
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title(""));
   c1->SetMargin(0.13,0.1,0.13,0.1);
   TGaxis::SetMaxDigits(3);
   xframe->SetTitle("");
   xframe->GetXaxis()->SetTitle("p (GeV/c)");
   xframe->GetXaxis()->SetNdivisions(505);
   xframe->GetXaxis()->SetTitleSize(0.05);
   xframe->GetXaxis()->SetLabelSize(0.05);
   xframe->GetXaxis()->SetTitleOffset(1.1);
   xframe->GetYaxis()->SetNdivisions(505);
   xframe->GetYaxis()->SetTitleSize(0.05);
   xframe->GetYaxis()->SetTitleOffset(1.25);
   xframe->GetYaxis()->SetLabelSize(0.05);
   xframe->GetYaxis()->SetNoExponent(0);
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(sig,bck),RooArgList(signal,background));
   Npar = 10;
   //x.setRange("sigragi",peakvalue-0.1,beame+0.012);
   //sum->fitTo(*dataset);
   sum->fitTo(*dataset);
  
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bck),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt;
   if (Ecm<2.65) pt = new TPaveText(0.5,0.55,0.85,0.90,"BRNDC");
   else pt = new TPaveText(0.15,0.55,0.50,0.90,"BRNDC");
   pt->SetBorderSize(0);
   //pt->SetFillStyle(4000);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.045);
 //sprintf(tmpchr,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError());
 //pt->AddText(tmpchr);
 //sprintf(tmpchr,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError());
 //pt->AddText(tmpchr);
   sprintf(tmpchr,"signal = %.1f #pm %.1f",signal.getVal(),signal.getError());
   //sprintf(tmpchr,"signal = %.1f #pm %.1f",Nsig_mean,Esig_mean);
   pt->AddText(tmpchr);
   //sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,xframe->chiSquare(Npar));
   sprintf(tmpchr,"#chi^{2} = %5.3f",xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
 //sprintf(tmpchr,"p_spectrum_%s",namesfx);
 //c1->SetName(tmpchr);
 //c1->Write();
   sprintf(tmpchr,"p_spectrum_%s_smcs.pdf",namesfx);
   c1->Print(tmpchr);
   
   ofstream ofsig("signal.txt",ios::app);
   ofsig << Ecm <<"\t"<< signal.getVal() <<"\t"<<signal.getError()<<endl;

 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   ofstream fitpar("fitpar",std::ios::app);
 //fitpar<<" ene = "<< beame*2 <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
 //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
 //fitpar<<"\t sigNo = " << signal.getVal() << "\t sigNoE = " << signal.getError();
 //fitpar<<"\t bckNo = " << background.getVal() << "\t bckNoE = " << background.getError();
 //fitpar<< "\t chi2 = " << xframe->chiSquare(Npar);
 //fitpar<< std::endl;
   fitpar<<setprecision(4);
   fitpar<< beame*2 <<"\t signal with mcshape \n";
   fitpar<< setprecision(6);
   fitpar<< meanb.getVal()<<"\t";
   fitpar<< sigmab.getVal()<<"\t";
 
   fitpar<< endl;
   
   ofstream fitpar_sigma("fitpar_sigsmr_sigma.txt", std::ios::app);
   fitpar_sigma<<setprecision(4);
   fitpar_sigma<< beame*2;
   fitpar_sigma<< "\t"<<sigma.getVal()<<"\t"<<sigma.getError()<<endl;
   
   ofstream fitpar_wid("fitpar_bckwid_sigma.txt", std::ios::app);
   fitpar_wid<<setprecision(4);
   fitpar_wid<< beame*2;
   fitpar_wid<< "\t"<<sigmae.getVal()<<"\t"<<sigmae.getError();
   fitpar_wid<< "\t"<<sigmab.getVal()<<"\t"<<sigmab.getError();
   fitpar_wid<< "\t"<<frac2.getVal()<<"\t"<<frac2.getError()<<endl;

   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete c1;
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
   //return Nsig_mean;
}

double FitSpectrum_smcsAbckgaus(TTree *dataraw, TH1D* hsig, double beame, const char* namesfx, double *par_pool)
{
   TF1 fsigma("fsigma","-0.00705519+0.00600716*x",2.0,3.1);
   double sigma_ini = fsigma.Eval(beame);
   TCanvas *c1 = new TCanvas("","",800,600);
   
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double mmu = 0.1057;
   double Ecm = beame;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   double pmu = sqrt(pow(beame,2)-pow(mmu,2));
   std::cout<<"Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;
   double beamlow=peakvalue-0.1;
   double beamup=peakvalue+0.2;
   if (Ecm<2.6) {
     beamlow=peakvalue-0.05;
     beamup=peakvalue+0.15;
   }
   
   char hname[1000];
 //TFile *file=new TFile("output/output_bck.root");
 //sprintf(hname,"hp_mcKK2_%.4f",Ecm);
 //TH1D *hk = (TH1D*)file->Get(hname);
 //sprintf(hname,"hp_mcmumu_%.4f",Ecm);
 //TH1D *hmu = (TH1D*)file->Get(hname);
 //cout<<"aaaaaaaaa"<<endl;
 
   // try to use roofit
   //
   RooRealVar x("p1","momentum",peakvalue,beamlow,beamup,"GeV/c");
   // signal
 //RooDataHist mckhist("mcKhist","mc K spec",x,hk);
   RooDataHist mckhist("mcKhist","mc K spec",x,hsig);
   RooHistPdf kpdf("kpdf","signal",x,mckhist,4); 
   RooRealVar mean("mean","mean of gaussian",0.0015,-0.005, 0.005);
   RooRealVar sigma("sigma","width of gaussian",0.002,0.0005,0.003);
   if (par_pool!=0){
     mean.setRange(par_pool[0],par_pool[0]);
     mean.setVal(par_pool[0]);
     sigma.setRange(par_pool[2],par_pool[2]);
     sigma.setVal(par_pool[2]);
   }
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooFFTConvPdf sig("sig","signal",x,kpdf,gaus); 


   // background, mainly di-mu
   RooRealVar meanb("meanb","mean of gaussian",pmu,pmu-0.005,pmu+0.005);
   RooRealVar sigmab("sigmab","width of gaussian",sigma_ini, sigma_ini-0.001,sigma_ini+0.001);
   RooRealVar meane("meane","mean of gaussian",pmu,pmu-0.03,pmu+0.005);
   RooRealVar sigmae("sigmae","width of gaussian",0.01,0.005,0.1);
   RooGaussian gausb("gausb","gauss(x,m,s)",x,meanb,sigmab);
   RooGaussian gause("gause","gauss(x,m,s)",x,meane,sigmae);
   RooRealVar frac2("frac2","frac2",0.8,0.7,1.0);
   
   RooAddPdf bck("bck","signal",RooArgList(gausb,gause),RooArgList(frac2));
   
   RooRealVar signal("signal"," ",1000,0,100000);//event number
   RooRealVar background("background"," ",2000,0,100000);
   
   // fix some parameter
 //double sigma_f = -0.0279011 + 0.0234759*Ecm -0.0044394*Ecm*Ecm;
 //sigma.setRange(sigma_f, sigma_f);
 //sigma.setVal(sigma_f);
 //double sigmab_f = -0.00286479 + 0.00447862*Ecm;
 //sigmab.setRange(sigmab_f, sigmab_f);
 //sigmab.setVal(sigmab_f);
 //double sigmae_f = -0.063683 + 0.0337589*Ecm;
 //sigmae.setRange(sigmae_f, sigmae_f);
 //sigmae.setVal(sigmae_f);
  
 
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title(""));
   c1->SetMargin(0.13,0.1,0.13,0.1);
   TGaxis::SetMaxDigits(3);
   xframe->SetTitle("");
   xframe->GetXaxis()->SetTitle("p (GeV/c)");
   xframe->GetXaxis()->SetNdivisions(505);
   xframe->GetXaxis()->SetTitleSize(0.05);
   xframe->GetXaxis()->SetLabelSize(0.05);
   xframe->GetXaxis()->SetTitleOffset(1.1);
   xframe->GetYaxis()->SetNdivisions(505);
   xframe->GetYaxis()->SetTitleSize(0.05);
   xframe->GetYaxis()->SetTitleOffset(1.25);
   xframe->GetYaxis()->SetLabelSize(0.05);
   xframe->GetYaxis()->SetNoExponent(0);
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(sig,bck),RooArgList(signal,background));
   Npar = 10;
   //x.setRange("sigragi",peakvalue-0.1,beame+0.012);
   //sum->fitTo(*dataset);
   sum->fitTo(*dataset);
  
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bck),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt;
   if (Ecm<2.65) pt = new TPaveText(0.5,0.55,0.85,0.90,"BRNDC");
   else pt = new TPaveText(0.15,0.55,0.50,0.90,"BRNDC");
   pt->SetBorderSize(0);
   //pt->SetFillStyle(4000);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.045);
 //sprintf(tmpchr,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError());
 //pt->AddText(tmpchr);
 //sprintf(tmpchr,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError());
 //pt->AddText(tmpchr);
   sprintf(tmpchr,"signal = %.1f #pm %.1f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"background = %.1f #pm %.1f",background.getVal(),background.getError());
   pt->AddText(tmpchr);
   //sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,xframe->chiSquare(Npar));
   sprintf(tmpchr,"#chi^{2} = %5.3f",xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
 //sprintf(tmpchr,"p_spectrum_%s",namesfx);
 //c1->SetName(tmpchr);
 //c1->Write();
   sprintf(tmpchr,"p_spectrum_%s_smcs.pdf",namesfx);
   c1->Print(tmpchr);
   
   ofstream ofsig("signal.txt",ios::app);
   ofsig << Ecm <<"\t"<< signal.getVal() <<"\t"<<signal.getError()<<endl;

 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   ofstream fitpar("fitpar",std::ios::app);
 //fitpar<<" ene = "<< beame*2 <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
 //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
 //fitpar<<"\t sigNo = " << signal.getVal() << "\t sigNoE = " << signal.getError();
 //fitpar<<"\t bckNo = " << background.getVal() << "\t bckNoE = " << background.getError();
 //fitpar<< "\t chi2 = " << xframe->chiSquare(Npar);
 //fitpar<< std::endl;
   fitpar<<setprecision(4);
   fitpar<< beame*2 <<"\t signal with mcshape \n";
   fitpar<< setprecision(6);
   fitpar<< meanb.getVal()<<"\t";
   fitpar<< sigmab.getVal()<<"\t";
 
   fitpar<< endl;
   
   ofstream fitpar_sigma("fitpar_sigsmr_sigma.txt", std::ios::app);
   fitpar_sigma<<setprecision(4);
   fitpar_sigma<< beame*2;
   fitpar_sigma<< "\t"<<sigma.getVal()<<"\t"<<sigma.getError()<<endl;
   
   ofstream fitpar_wid("fitpar_bckwid_sigma.txt", std::ios::app);
   fitpar_wid<<setprecision(4);
   fitpar_wid<< beame*2;
   fitpar_wid<< "\t"<<sigmae.getVal()<<"\t"<<sigmae.getError();
   fitpar_wid<< "\t"<<sigmab.getVal()<<"\t"<<sigmab.getError()<<endl;

   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete c1;
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
   //return Nsig_mean;
}

double FitSpectrum_bmcs(TTree *dataraw, TH1D* hbck, double beame, const char* namesfx, double *par_pool)
{
   TF1 fsigma("fsigma","-0.00705519+0.00600716*x",2.0,3.1);
   double sigma_ini = fsigma.Eval(beame);
   TCanvas *c1 = new TCanvas("","",800,600);
   
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double mmu = 0.1057;
   double Ecm = beame;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   double pmu = sqrt(pow(beame,2)-pow(mmu,2));
   std::cout<<"Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;

   double beamlow=peakvalue-0.1;
   double beamup=peakvalue+0.2;
   if (Ecm<2.6) {
     beamlow=peakvalue-0.05;
     beamup=peakvalue+0.15;
   }
   
   char hname[1000];
 //TFile *file=new TFile("output/output_bck.root");
 //sprintf(hname,"hp_mcKK2_%.4f",Ecm);
 //TH1D *hk = (TH1D*)file->Get(hname);
 //sprintf(hname,"hp_mcmumu_%.4f",Ecm);
 //TH1D *hmu = (TH1D*)file->Get(hname);
// TH1D *hmurebin = (TH1D*)hmu->Rebin(2,hname);
 
   //
   // try to use roofit
   //
   RooRealVar x("p1","momentum",peakvalue,beamlow,beamup,"GeV/c");
   // signal
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.005,peakvalue+0.005);
   RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.013);
   //RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.014);
   RooRealVar alpha1("alpha1","#alpha",1.5,1.0,5.0);
   RooRealVar nnn1("n1","n",2,1,10);
   RooCBShape cbshape("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);
   RooRealVar sigma2("sigma2","width of gaussian",0.007,0.005,0.02);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma2);
   RooRealVar frac1("frac1","frac1",0.7,0.6,0.9);

   RooAddPdf sig("sig","signal",RooArgList(cbshape,gaus),RooArgList(frac1));
   
   // background, mainly di-mu
   cout<<"aaaaaaaaa"<<endl;
   //RooDataHist mcmuhist("mcmuhist","mc mu spec",x,hmu);
   RooDataHist mcmuhist("mcmuhist","mc mu spec",x,hbck);
   RooHistPdf mupdf("mupdf","background",x,mcmuhist,4);
   //RooHistPdf bck("bck","background",x,mcmuhist,5);
   RooRealVar mean_sp("mean_sp","mean of gaussian",0.002,-0.005, 0.005);
   RooRealVar sigma_sp("sigma_sp","width of gaussian",0.0024,0.0005,0.003);
   //RooRealVar sigma_sp("sigma_sp","width of gaussian",0.0024,0.0005,0.005);
   RooGaussian gaus_sp("gaus_sp","gauss(x,m,s)",x,mean_sp,sigma_sp);
   if (par_pool!=0){
     mean_sp.setRange(par_pool[0],par_pool[0]);
     mean_sp.setVal(par_pool[0]);
     sigma_sp.setRange(par_pool[2],par_pool[2]);
     sigma_sp.setVal(par_pool[2]);
   }


   // fix some parameters
 //double sigmasp_f = 0.003;
 //sigma_sp.setRange(sigmasp_f, sigmasp_f);
 //sigma_sp.setVal(sigmasp_f);
 //double sigma_f = -0.00497225 + 0.0046221*Ecm;
 //sigma.setRange(sigma_f, sigma_f);
 //sigma.setVal(sigma_f);
 //double sigma2_f = 0.0107904;
 //sigma2.setRange(sigma2_f, sigma2_f);
 //sigma2.setVal(sigma2_f);

   // background, mainly di-mu
   RooFFTConvPdf bck("bck","background",x,mupdf,gaus_sp); 
   //RooHistPdf bck("bck","background",x,mcmuhist,2);
 
   
   RooRealVar signal("signal"," ",1000,0,1000000000);//event number
   RooRealVar background("background"," ",20000,0,100000000);
 
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title(""));
   c1->SetMargin(0.13,0.1,0.13,0.1);
   TGaxis::SetMaxDigits(3);
   xframe->SetTitle("");
   xframe->GetXaxis()->SetTitle("p (GeV/c)");
   xframe->GetXaxis()->SetNdivisions(505);
   xframe->GetXaxis()->SetTitleSize(0.05);
   xframe->GetXaxis()->SetLabelSize(0.05);
   xframe->GetXaxis()->SetTitleOffset(1.1);
   xframe->GetYaxis()->SetNdivisions(505);
   xframe->GetYaxis()->SetTitleSize(0.05);
   xframe->GetYaxis()->SetTitleOffset(1.25);
   xframe->GetYaxis()->SetLabelSize(0.05);
   xframe->GetYaxis()->SetNoExponent(0);
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(sig,bck),RooArgList(signal,background));
   Npar = 10;
   //x.setRange("sigragi",peakvalue-0.1,beame+0.012);
   //sum->fitTo(*dataset);
   sum->fitTo(*dataset);
   
   
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bck),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt;
   if (Ecm<2.65) pt = new TPaveText(0.5,0.55,0.85,0.90,"BRNDC");
   else pt = new TPaveText(0.15,0.55,0.50,0.90,"BRNDC");
   pt->SetBorderSize(0);
   //pt->SetFillStyle(4000);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.045);
   sprintf(tmpchr,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"signal = %.1f #pm %.1f",signal.getVal(),signal.getError());
   //sprintf(tmpchr,"signal = %.1f #pm %.1f",Nsig_mean,Esig_mean);
   pt->AddText(tmpchr);
   //sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,xframe->chiSquare(Npar));
   sprintf(tmpchr,"#chi^{2} = %5.3f",xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
// sprintf(tmpchr,"p_spectrum_%s",namesfx);
// c1->SetName(tmpchr);
// c1->Write();
   sprintf(tmpchr,"p_spectrum_%s_bmcs.pdf",namesfx);
   c1->Print(tmpchr);

   ofstream ofsig("signal.txt",ios::app);
   ofsig << Ecm <<"\t"<< signal.getVal() <<"\t"<<signal.getError()<<endl;
 
 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   ofstream fitpar("fitpar",std::ios::app);
 //fitpar<<" ene = "<< beame*2 <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
 //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
 //fitpar<<"\t sigNo = " << signal.getVal() << "\t sigNoE = " << signal.getError();
 //fitpar<<"\t bckNo = " << background.getVal() << "\t bckNoE = " << background.getError();
 //fitpar<< "\t chi2 = " << xframe->chiSquare(Npar);
 //fitpar<< std::endl;
   fitpar<<setprecision(4);
   fitpar<< beame*2 <<"\t bacground with mcshape \n";
   
   fitpar<<setprecision(6);
   fitpar<< mean.getVal() <<"\t ";
   fitpar<< sigma.getVal() <<"\t ";
 //fitpar<< sigma2.getVal() <<"\t ";
 //fitpar<< alpha1.getVal()<<"\t ";
 //fitpar<< nnn1.getVal()<<"\t ";
 //fitpar<< frac1.getVal()<<"\t ";
   
  
   fitpar<< endl;
   
   ofstream fitpar_wid("fitpar_sigwid_sigma.txt", std::ios::app);
   fitpar_wid<<setprecision(4);
   fitpar_wid<< beame*2;
   fitpar_wid<< "\t"<<sigma2.getVal()<<"\t"<<sigma2.getError();
   fitpar_wid<< "\t"<<sigma.getVal()<<"\t"<<sigma.getError();
   fitpar_wid<< "\t"<<frac1.getVal()<<"\t"<<frac1.getError()<<endl;
   
   ofstream fitpar_sigma("fitpar_bcksmr_sigma.txt", std::ios::app);
   fitpar_sigma<<setprecision(4);
   fitpar_sigma<< beame*2;
   fitpar_sigma<< "\t"<<sigma_sp.getVal()<<"\t"<<sigma_sp.getError()<<endl;



   double intel3 = mean.getVal()-3*sigma.getVal();
   double inteu3 = mean.getVal()+3*sigma.getVal();
   double intel5 = mean.getVal()-5*sigma.getVal();
   double inteu5 = mean.getVal()+5*sigma.getVal();
   //RooRealVar xint("xint","xint",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame+0.01);
   x.setRange("sigrag3",intel3, inteu3);
   x.setRange("sigrag5",intel5, inteu5);
   RooAbsReal* intsigi = sig.createIntegral( x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intbcki = bck.createIntegral(x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intsig3 = sig.createIntegral( x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intbck3 = bck.createIntegral(x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intsig5 = sig.createIntegral( x,  NormSet(x),  Range("sigrag5"));
   RooAbsReal* intbck5 = bck.createIntegral(x,  NormSet(x),  Range("sigrag5"));
   double signalN3 =   signal.getVal()*intsig3->getVal();
   double backN3 = background.getVal()*intbck3->getVal();
   double signalN5 =   signal.getVal()*intsig5->getVal();
   double backN5 = background.getVal()*intbck5->getVal();
 //double signalN3nom =   signal.getVal()*intsig3->getVal()/intsigi->getVal();
 //double backN3nom = background.getVal()*intbck3->getVal()/intbcki->getVal();
 //double signalN5nom =   signal.getVal()*intsig5->getVal()/intsigi->getVal();
 //double backN5nom = background.getVal()*intbck5->getVal()/intbcki->getVal();
   ofstream angSNR("angSNR.dat",std::ios::app);
   angSNR << namesfx <<"\t"<<signalN3 << "\t" << backN3 <<"\t"<<signalN5 << "\t" << backN5 << std::endl;
 //epSNR << "\t"<<signalN3nom << "\t" << backN3nom <<"\t"<<signalN5nom << "\t" << backN5nom<< std::endl;
   std::cout<< "Total signal int is "<< intsigi->getVal() <<" Total bck int is "<< intbcki->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig3->getVal() <<" Total bck int is "<< intbck3->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig5->getVal() <<" Total bck int is "<< intbck5->getVal()<<std::endl;

   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete c1;
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
   //return Nsig_mean;
}

double FitSpectrum_bmcsAsiggaus(TTree *dataraw, TH1D* hbck, double beame, const char* namesfx, double *par_pool)
{
   TF1 fsigma("fsigma","-0.00705519+0.00600716*x",2.0,3.1);
   double sigma_ini = fsigma.Eval(beame);
   TCanvas *c1 = new TCanvas("","",800,600);
   
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double mmu = 0.1057;
   double Ecm = beame;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   double pmu = sqrt(pow(beame,2)-pow(mmu,2));
   std::cout<<"Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;

   double beamlow=peakvalue-0.1;
   double beamup=peakvalue+0.2;
   if (Ecm<2.6) {
     beamlow=peakvalue-0.05;
     beamup=peakvalue+0.15;
   }
   
   char hname[1000];
 //TFile *file=new TFile("output/output_bck.root");
 //sprintf(hname,"hp_mcKK2_%.4f",Ecm);
 //TH1D *hk = (TH1D*)file->Get(hname);
 //sprintf(hname,"hp_mcmumu_%.4f",Ecm);
 //TH1D *hmu = (TH1D*)file->Get(hname);
// TH1D *hmurebin = (TH1D*)hmu->Rebin(2,hname);
 
   //
   // try to use roofit
   //
   RooRealVar x("p1","momentum",peakvalue,beamlow,beamup,"GeV/c");
   // signal
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.005,peakvalue+0.005);
   RooRealVar sigma("sigma","width of gaussian",0.004,0.0035,0.012);
   //RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.014);
   RooRealVar mean2("mean2","mean of gaussian",peakvalue-0.002,peakvalue-0.005,peakvalue+0.005);
   RooRealVar sigma2("sigma2","width of gaussian",0.007,0.005,0.03);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean2,sigma2);
   RooRealVar frac1("frac1","frac1",0.8,0.75,1.0);

   RooAddPdf sig("sig","signal",RooArgList(gaus,gaus2),RooArgList(frac1));
   
   // background, mainly di-mu
   cout<<"aaaaaaaaa"<<endl;
   //RooDataHist mcmuhist("mcmuhist","mc mu spec",x,hmu);
   RooDataHist mcmuhist("mcmuhist","mc mu spec",x,hbck);
   RooHistPdf mupdf("mupdf","background",x,mcmuhist,4);
   //RooHistPdf bck("bck","background",x,mcmuhist,5);
   RooRealVar mean_sp("mean_sp","mean of gaussian",0.002,-0.005, 0.005);
   RooRealVar sigma_sp("sigma_sp","width of gaussian",0.0024,0.0005,0.003);
   //RooRealVar sigma_sp("sigma_sp","width of gaussian",0.0024,0.0005,0.005);
   RooGaussian gaus_sp("gaus_sp","gauss(x,m,s)",x,mean_sp,sigma_sp);
   if (par_pool!=0){
     mean_sp.setRange(par_pool[0],par_pool[0]);
     mean_sp.setVal(par_pool[0]);
     sigma_sp.setRange(par_pool[2],par_pool[2]);
     sigma_sp.setVal(par_pool[2]);
   }


   // fix some parameters
 //double sigmasp_f = 0.003;
 //sigma_sp.setRange(sigmasp_f, sigmasp_f);
 //sigma_sp.setVal(sigmasp_f);
 //double sigma_f = -0.00497225 + 0.0046221*Ecm;
 //sigma.setRange(sigma_f, sigma_f);
 //sigma.setVal(sigma_f);
 //double sigma2_f = 0.0107904;
 //sigma2.setRange(sigma2_f, sigma2_f);
 //sigma2.setVal(sigma2_f);

   // background, mainly di-mu
   RooFFTConvPdf bck("bck","background",x,mupdf,gaus_sp); 
   //RooHistPdf bck("bck","background",x,mcmuhist,2);
 
   
   RooRealVar signal("signal"," ",1000,0,1000000000);//event number
   RooRealVar background("background"," ",20000,0,100000000);
 
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title(""));
   c1->SetMargin(0.13,0.1,0.13,0.1);
   TGaxis::SetMaxDigits(3);
   xframe->SetTitle("");
   xframe->GetXaxis()->SetTitle("p (GeV/c)");
   xframe->GetXaxis()->SetNdivisions(505);
   xframe->GetXaxis()->SetTitleSize(0.05);
   xframe->GetXaxis()->SetLabelSize(0.05);
   xframe->GetXaxis()->SetTitleOffset(1.1);
   xframe->GetYaxis()->SetNdivisions(505);
   xframe->GetYaxis()->SetTitleSize(0.05);
   xframe->GetYaxis()->SetTitleOffset(1.25);
   xframe->GetYaxis()->SetLabelSize(0.05);
   xframe->GetYaxis()->SetNoExponent(0);
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(sig,bck),RooArgList(signal,background));
   Npar = 10;
   //x.setRange("sigragi",peakvalue-0.1,beame+0.012);
   //sum->fitTo(*dataset);
   sum->fitTo(*dataset);
      
   
   
   
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bck),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt;
   if (Ecm<2.65) pt = new TPaveText(0.5,0.55,0.85,0.90,"BRNDC");
   else pt = new TPaveText(0.15,0.55,0.50,0.90,"BRNDC");
   pt->SetBorderSize(0);
   //pt->SetFillStyle(4000);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.045);
   sprintf(tmpchr,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"signal = %.1f #pm %.1f",signal.getVal(),signal.getError());
   //sprintf(tmpchr,"signal = %.1f #pm %.1f",Nsig_mean,Esig_mean);
   pt->AddText(tmpchr);
   //sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,xframe->chiSquare(Npar));
   sprintf(tmpchr,"#chi^{2} = %5.3f",xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
// sprintf(tmpchr,"p_spectrum_%s",namesfx);
// c1->SetName(tmpchr);
// c1->Write();
   sprintf(tmpchr,"p_spectrum_%s_bmcs.pdf",namesfx);
   c1->Print(tmpchr);

   ofstream ofsig("signal.txt",ios::app);
   ofsig << Ecm <<"\t"<< signal.getVal() <<"\t"<<signal.getError()<<endl;
 
 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   ofstream fitpar("fitpar",std::ios::app);
 //fitpar<<" ene = "<< beame*2 <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
 //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
 //fitpar<<"\t sigNo = " << signal.getVal() << "\t sigNoE = " << signal.getError();
 //fitpar<<"\t bckNo = " << background.getVal() << "\t bckNoE = " << background.getError();
 //fitpar<< "\t chi2 = " << xframe->chiSquare(Npar);
 //fitpar<< std::endl;
   fitpar<<setprecision(4);
   fitpar<< beame*2 <<"\t bacground with mcshape and signal with gaus\n";
   
   fitpar<<setprecision(6);
   fitpar<< mean.getVal() <<"\t ";
   fitpar<< sigma.getVal() <<"\t ";
 //fitpar<< sigma2.getVal() <<"\t ";
 //fitpar<< alpha1.getVal()<<"\t ";
 //fitpar<< nnn1.getVal()<<"\t ";
 //fitpar<< frac1.getVal()<<"\t ";
   
  
   fitpar<< endl;
   
// ofstream fitpar_wid("fitpar_sigwid_sigma.txt", std::ios::app);
// fitpar_wid<<setprecision(4);
// fitpar_wid<< beame*2;
// fitpar_wid<< "\t"<<sigma2.getVal()<<"\t"<<sigma2.getError();
// fitpar_wid<< "\t"<<sigma.getVal()<<"\t"<<sigma.getError()<<endl;
   
// ofstream fitpar_sigma("fitpar_bcksmr_sigma.txt", std::ios::app);
// fitpar_sigma<<setprecision(4);
// fitpar_sigma<< beame*2;
// fitpar_sigma<< "\t"<<sigma_sp.getVal()<<"\t"<<sigma_sp.getError()<<endl;



   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete c1;
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
   //return Nsig_mean;
}





double m_par[12];

bool getPars(double energy)
{
  ifstream par("mcpar.txt");
  double ene[22];
  double parset[22][12];
  for (int iene=0; iene<22; iene++){
    char line[1000];
    par.getline(line,1000);
    istringstream iss;
    iss.str(line);
    iss >> ene[iene];
    for (int i= 0; i<8; i++){
      iss >> parset[iene][i];
    }
  }
  par.close();
  
  for (int iene=0;  iene<22; iene++){
     if (fabs(energy-ene[iene])<1e-4)  { for (int i=0; i<8; i++) m_par[i]= parset[iene][i]; return 1;}
  }
  return false;
}
double FitSpectrum_mcshape_fitpar(TTree *&dataraw, TH1D* hsig, TH1D* hbck, double beame, double &err, const char* namesfx, TFile* fileout)
{
   TCanvas *c1 = new TCanvas("","",800,600);
   
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double mmu = 0.1057;
   double Ecm = beame;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   double pmu = sqrt(pow(beame,2)-pow(mmu,2));
   std::cout<<"Ecm is "<<Ecm<<", Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;
   double beamlow=peakvalue-0.1;
   double beamup=peakvalue+0.2;
   if (Ecm<2.6) {
     beamlow=peakvalue-0.05;
     beamup=peakvalue+0.15;
   }
   
   if (fileout!=0) fileout->cd();
   //
   // try to use roofit
   //
   cout<<"aaaaaaaaa"<<endl;
   RooRealVar x("p1","p",beamlow,beamup,"GeV/c");
    // signal
 
   cout<<"aaaaaaaaa"<<endl;
   RooDataHist mckhist("mcKhist","mc K spec",x,hsig);
   RooDataHist mcmuhist("mcmuhist","mc mu spec",x,hbck);
   cout<<"aaaaaaaaa"<<endl;
   RooHistPdf kpdf("kpdf","signal",x,mckhist,4); 
   RooHistPdf mupdf("mupdf","background",x,mcmuhist,4);
   cout<<"aaaaaaaaa"<<endl;

   getPars(Ecm);
   RooRealVar mean_k("mean_k","mean of gaussian",m_par[0]);
   RooRealVar sigma_k("sigma_k","width of gaussian",m_par[2]);
   RooGaussian gaus_k("gaus_k","gauss(x,m,s)",x,mean_k,sigma_k);

   RooRealVar mean_u("mean","mean of gaussian",m_par[4]);
   RooRealVar sigma_u("sigma","width of gaussian",m_par[6]+m_par[7]);
   RooGaussian gaus_u("gaus_u","gauss(x,m,s)",x,mean_u,sigma_u);
   
   // signal
   RooFFTConvPdf sig("sig","signal",x,kpdf,gaus_k); 
   //RooHistPdf sig("sig","signal",x,mckhist,2); 
   // background, mainly di-mu
   RooFFTConvPdf bck("bck","background",x,mupdf,gaus_u); 
   //RooHistPdf bck("bck","background",x,mcmuhist,2);
   
   RooRealVar signal("signal"," ",1000,0,100000);//event number
   RooRealVar background("background"," ",200,0,100000);
 
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title("fit p"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(sig,bck),RooArgList(signal,background));
   //sum = new RooAddPdf("sum","sum",RooArgList(sig),RooArgList(signal));
   Npar = 4;
   sum->fitTo(*dataset);
   err = signal.getError();
   dataset->plotOn(xframe,Binning(nBins));
   sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bck),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt = new TPaveText(0.15,0.55,0.45,0.90,"BRNDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(4000);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.035);
   sprintf(tmpchr,"mean_{K} = %1.6f #pm %1.6f",mean_k.getVal(),mean_k.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma_{K} = %1.6f #pm %1.6f",sigma_k.getVal(),sigma_k.getError());
   pt->AddText(tmpchr);
 //sprintf(tmpchr,"mean_{u} = %1.6f #pm %1.6f",mean_u.getVal(),mean_u.getError());
 //pt->AddText(tmpchr);
 //sprintf(tmpchr,"#sigma_{u} = %1.6f #pm %1.6f",sigma_u.getVal(),sigma_u.getError());
 //pt->AddText(tmpchr);
   sprintf(tmpchr,"N_{sig} = %.2f #pm %.2f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"N_{bck} = %.2f #pm %.2f",background.getVal(),background.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
   sprintf(tmpchr,"p_spectrum_fixpar_%s.pdf",namesfx);
   c1->Print(tmpchr);
   
   ofstream ofsig("signal.txt",ios::app);
   ofsig << Ecm <<"\t"<< signal.getVal() <<"\t"<<signal.getError()<<endl;
   
   delete xframe;
   delete dataset;
   delete sum;
   delete c1; 
   return signal.getVal();
}

double GetEnergy(int run)
{
  run = abs(run);
  if (run>=39335 && run<=39618) return 3.08;
  if (run>=39711 && run<=39738) return 3.02;
  if (run>=39680 && run<=39710) return 3.00;
  if (run>=39651 && run<=39679) return 2.981;
  if (run>=39619 && run<=39650) return 2.95;
  
  if (run>=39775 && run<=40069) return 2.90;
  if (run>=40128 && run<=40296) return 2.6444;
  if (run>=40300 && run<=40435) return 2.6464;
  if (run>=40436 && run<=40439) return 2.70;
  if (run>=40440 && run<=40443) return 2.80;
  if (run>=40459 && run<=40769) return 2.396;
  if (run>=40771 && run<=40776) return 2.5;
  if (run>=40777 && run<=40804) return 2.6444;//separated beam
  if (run>=40806 && run<=40951) return 2.3864;
  if (run>=40989 && run<=41121) return 2.2;
  if (run>=41122 && run<=41239) return 2.2324;
  if (run>=41240 && run<=41411) return 2.3094;
  if (run>=41416 && run<=41532) return 2.175;
  if (run>=41533 && run<=41570) return 2.15;
  if (run>=41588 && run<=41727) return 2.1;

  if (run>=41729 && run<=41909) return 2.0;
  if (run>=41911 && run<=41958) return 2.05;
  if (run>=41959 && run<=41999) return 2.2324; // separated beam
  
  if (run>=42004) return 2.125;
  if (run>=27147 && run<=27288) return 3.08;
  if (run>=28241 && run<=28266) return 3.08;
  if (run>=28624 && run<=28648) return 2.2324;
  if (run>=28553 && run<=28575) return 2.8;
  return -1;
}


const char* getPureName(const char* name)
{
  int pos=0;
  for (int i=0;;i++){
    if (name[i]=='\0') return &name[pos];
    if (name[i]=='/') pos = i+1;
    if (i>10000) return NULL;
  } 
}

const char* getPureName2(const char* name)
{
  char *name1 = new char[1000];
  int pos=0;
  int ppos=0;

  for (int i=0;;i++){
    if (name[i]=='\0') break; //return &name[pos];
    if (name[i]=='/') pos = i+1;
    if (name[i]=='.') ppos = i;
    if (i>1000) return NULL;
  }
  //std::cout<<"/ pos is "<< pos <<", . pos is "<< ppos <<std::endl;
  for (int i=0;i<ppos-pos;i++){
    name1[i] = name[pos+i];
  }
  name1[ppos-pos]='\0';

  return name1;
}


double getEne(const char* name)
{
  name = getPureName(name);
  std::string str = name;
  int _pos=0;
  int ppos=0;
  for (int i=0;;i++){
    if (name[i]=='\0') break;
    if (name[i]=='_') _pos=i;
    if (name[i]=='.') ppos=i;
  }
  int len = ppos-_pos-1;
  string vstr=str.substr(_pos+1,len);
  if (atof(vstr.c_str())/10000>1) return atof(vstr.c_str())/10000;

  else {// find index number in _*_
    int _pos2[2]={0,0};
    int j=0;
    for (int i=0;;i++){
      if (name[i]=='_') {_pos2[j]=i; j++;}
      if (j>=2) break;
    }
    if (j!=2) return -1;
    string vstr = str.substr(_pos2[0]+1,3);
    int idx = atoi(vstr.c_str());
    std::cout<<vstr<<'\t'<<idx<<std::endl;
    switch (idx) {
      case 1:    return  2.0;   
      case 2:    return  2.05;  
      case 3:    return  2.1;   
      case 4:    return  2.15;  
      case 5:    return  2.175; 
                                
      case 6:    return  2.2;   
      case 7:    return  2.2324;
      case 8:    return  2.3094;
      case 9:    return  2.3864;
      case 10:   return  2.396; 
      case 11:   return  2.5;   
      case 12:   return  2.6444;
      case 13:   return  2.6464;
      case 14:   return  2.7;   
      case 15:   return  2.8;   
      case 16:   return  2.9;   
      case 17:   return  2.95;  
      case 18:   return  2.981; 
      case 19:   return  3.0;   
      case 20:   return  3.02;  
                                
      case 21:   return  3.08;  
    }
  }
  return -2;
}

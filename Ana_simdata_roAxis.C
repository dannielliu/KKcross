//##############################
//### this code is modified ####
//### dliu13@mail.ustc.edu.cn###
//##############################
//##############################
//##############################


#include "TFile.h"
//#include "TFolder.h"
#include <TCanvas.h>
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
  double FitSpectrum3(TTree *dataraw, double beame, double &err, const char* namesfx=0);
//double FitSpectrum3(TEventList *dataraw, double beame, const char* namesfx=0);
double FitSpectrum4(TTree *&dataraw, double beame, const char* namesfx=0, TFile* fileout=0); // fit with MC shape
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
  while (!cutin.eof()){
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
  double sump1p2;
  double difp1p2;
  TTree *elist = new TTree("elist","elist");
  elist->Branch("p1",&p1,"p1/D");
  elist->Branch("p2",&p2,"p2/D");
  elist->Branch("sump1p2",&sump1p2,"sump1p2/D");
  elist->Branch("difp1p2",&difp1p2,"difp1p2/D");
 
  fileout->cd();

  int count0=0,count1=0,count2=0,count3=0;
  int count[10]={0};

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
    double p_exp = sqrt(pow(energys[iene]/2, 2)-pow(mka,2));
    double entries = var_data->GetEntries();
    for (int ien=0;ien<entries;ien++){
      var_data->GetEntry(ien);
      sump1p2 = (p1+p2)/2;
      difp1p2 = (p1-p2)/2;
      if (theta>179 && fabs(difp1p2)<0.04) elist->Fill();
    }
  //sprintf(name,"h2pini_%d",(int)(energys[iene]*10000));
  //TH2D *h2pini = new TH2D(name,name,100,p_exp-0.1,p_exp+0.2,100,p_exp-0.1,p_exp+0.2);
  //sprintf(name,"h2prot_%d",(int)(energys[iene]*10000));
  //TH2D *h2prot = new TH2D(name,name,100,p_exp-0.1,p_exp+0.2,100, -0.1, 0.1);
  //TCanvas *c1 = new TCanvas();
  //sprintf(name,"p1:p2>>h2pini_%d",(int)(energys[iene]*10000));
  //elist->Draw(name);
  //h2pini->Draw("colz");
  //sprintf(name,"h2pini_%d.pdf",(int)(energys[iene]*10000));
  //c1->Print(name);
  //sprintf(name,"difp1p2:sump1p2>>h2prot_%d",(int)(energys[iene]*10000));
  //elist->Draw(name);
  //h2prot->Draw("colz");
  //sprintf(name,"h2prot_%d.pdf",(int)(energys[iene]*10000));
  //c1->Print(name);
  //continue;
    
    
    ofstream cutflow("cutflow2",std::ios::app);
    //ofstream cutflow("cutflow_cmp665andp01",std::ios::app);
    cutflow<<energys[iene]<<" change fitting shape"<<std::endl;
    cutflow<<"Initial size      :"<<var_data->GetEntries()<<std::endl;
    cutflow<<"After theta cut   :"<<var_data->Draw("","theta<179")<<std::endl;
    cutflow<<"p2-exp<3 sigma    :"<<elist->GetEntries()<<std::endl;
    
    double err;
    sprintf(name,"KK_%d",(int)(energys[iene]*10000));
    double sigp1=FitSpectrum3(elist, energys[iene],err, name);
    cutflow<<"p1 fit signal pre :"<<sigp1<<" +/-"<<err <<std::endl;
  }
  // finish program here
  return 0;
 
}

#include <TStyle.h>
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
 //double beamlow=peakvalue-0.1;
 //double beamup=peakvalue+0.2;
 //if (Ecm<2.6) {
 //  beamlow=peakvalue-0.05;
 //  beamup=peakvalue+0.15;
 //}
   double beamlow=peakvalue-0.05;
   double beamup=peakvalue+0.05;
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
   RooRealVar x("sump1p2","momentum",peakvalue,beamlow,beamup,"GeV");
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


#include "RooHistPdf.h"
#include "RooFFTConvPdf.h"
double FitSpectrum4(TTree *&dataraw, double beame, const char* namesfx, TFile* fileout)
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
   sprintf(hname,"hkrebin_%.4f",Ecm);
   TH1D *hkrebin = (TH1D*)hk->Rebin(2,hname);
   sprintf(hname,"hmurebin_%.4f",Ecm);
   TH1D *hmurebin = (TH1D*)hmu->Rebin(2,hname);
   cout<<"hist dim "<< hkrebin->GetDimension()<<endl;
   cout<<"aaaaaaaaa hkrebin pointer "<<hkrebin<<", hmurebin "<<hmurebin<<endl;
   //return 0;
     if (hk==0 || hmu==0) {cout<<"Error : spare pointer!!!!"<< endl; exit(-1);}
     hkrebin->Draw();
   sprintf(hname,"hktmp_%.4f.pdf",Ecm);
     c1->Print(hname);
     hmurebin->Draw();
   sprintf(hname,"hmutmp_%.4f.pdf",Ecm);
     c1->Print(hname);
   //
   if (fileout!=0) fileout->cd();
   //
   // try to use roofit
   //
   cout<<"aaaaaaaaa"<<endl;
   RooRealVar x("p1","momentum",peakvalue,beamlow,beamup,"GeV");
    // signal
 
   cout<<"aaaaaaaaa"<<endl;
   RooDataHist mckhist("mcKhist","mc K spec",x,hkrebin);
   RooDataHist mcmuhist("mcmuhist","mc mu spec",x,hmurebin);
   cout<<"aaaaaaaaa"<<endl;
   RooHistPdf kpdf("kpdf","signal",x,mckhist,4); 
   RooHistPdf mupdf("mupdf","background",x,mcmuhist,4);
   cout<<"aaaaaaaaa"<<endl;

   RooRealVar mean("mean","mean of gaussian",0.0007,-0.005, 0.005);
   RooRealVar sigma("sigma","width of gaussian",0.0024,0.0005,0.003);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);

   // signal
   RooFFTConvPdf sig("sig","signal",x,kpdf,gaus); 
   //RooHistPdf sig("sig","signal",x,mckhist,2); 
   // background, mainly di-mu
   RooFFTConvPdf bck("bck","background",x,mupdf,gaus); 
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
   //x.setRange("sigragi",peakvalue-0.1,beame+0.012);
   //sum->fitTo(*dataset,Range("sigragi"));
   //sum->fitTo(*dataset,Minimizer("Minuit","simplex"),Hesse(false));
   //sum->fitTo(*dataset,Minimizer("GSLSimAn"));
   //sum->fitTo(*dataset,Minimizer("GSLMultiMin"),Hesse(false) );
   //sum->chi2FitTo(*dataset);
   sum->fitTo(*dataset);
   dataset->plotOn(xframe,Binning(nBins));
   //sum->fitTo(mckhist);
   //mckhist.plotOn(xframe);
   sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bck),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt = new TPaveText(0.15,0.65,0.45,0.90,"BRNDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(4000);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.035);
   sprintf(tmpchr,"#mu = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
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
   
   
   file->Close();
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
}




double m_par[12];

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

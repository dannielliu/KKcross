#include "TFile.h"
#include "TTree.h"
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
using namespace RooFit;

#include "RooHistPdf.h"
#include "RooFFTConvPdf.h"
//double FitSpectrum4(TTree *&dataraw, TH1D* hk, TH1D* hmu, double beame, const char* namesfx, TFile* fileout)
double FitSpectrum(TH1D* dataraw, TH1D* hk, TH1D* hmu, double beame, const char* namesfx, TFile* fileout=0)
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
   std::cout<<"Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;
   double beamlow=peakvalue-0.1;
   double beamup=peakvalue+0.2;
   if (Ecm<2.6) {
     beamlow=peakvalue-0.05;
     beamup=peakvalue+0.15;
   }

   TH1D *hkrebin = (TH1D*)hk->Rebin(2,"hkrebin");
   TH1D *hmurebin = (TH1D*)hmu->Rebin(2,"hmurebin");
   
  //
   //fileout->cd();
   //
   // try to use roofit
   //
   RooRealVar x("x1","momentum",peakvalue,beamlow,beamup,"GeV");
    // signal
 //RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.003,peakvalue+0.003);
 //RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.014);
 //RooRealVar alpha1("alpha1","#alpha",1.5,1.0,5.0);
 //RooRealVar nnn1("n1","n",2,1,10);
 //RooCBShape cbshape("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);
 //RooRealVar sigma2("sigma2","width of gaussian",0.007,0.005,0.02);
 //RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma2);
 //RooRealVar frac1("frac1","frac1",0.9,0.8,1.0);

// RooAddPdf sig("sig","signal",RooArgList(cbshape,gaus),RooArgList(frac1));
 
 
 
   RooDataHist mckhist("mcKhist","mc K spec",x,hkrebin);
   RooDataHist mcmuhist("mcmuhist","mc mu spec",x,hmurebin);
   RooHistPdf kpdf("kpdf","signal",x,mckhist,4); 
   RooHistPdf mupdf("mupdf","background",x,mcmuhist,4);

   RooRealVar mean("mean","mean of gaussian",0.0007,-0.005, 0.005);
   RooRealVar sigma("sigma","width of gaussian",0.001,0.00,0.010);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);

   // signal
   RooFFTConvPdf sig("sig","signal",x,kpdf,gaus); 
 //RooHistPdf sig("sig","signal",x,mckhist,4); 
   // background, mainly di-mu
   RooFFTConvPdf bck("bck","background",x,mupdf,gaus); 
 //RooHistPdf bck("bck","background",x,mcmuhist,4);
   
   RooRealVar signal("signal"," ",2000,0,100000);//event number
   RooRealVar background("background"," ",200,0,100000);
 
   RooAddPdf *sum;
   RooPlot *xframe;
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title("fit p"));
   //RooDataSet *dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   RooDataHist *dataset = new RooDataHist(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(sig,bck),RooArgList(signal,background));
   //sum = new RooAddPdf("sum","sum",RooArgList(sig),RooArgList(signal));
   Npar = 2;
   //x.setRange("sigragi",peakvalue-0.1,beame+0.012);
   //sum->fitTo(*dataset,Range("sigragi"));
   //sum->fitTo(*dataset,Minimizer("Minuit","simplex"),Hesse(false));
   //sum->fitTo(*dataset,Minimizer("GSLSimAn"));
   //sum->fitTo(*dataset,Minimizer("GSLMultiMin"),Hesse(false) );
   sum->chi2FitTo(*dataset);
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
   sprintf(tmpchr,"p_spectrum_%s",namesfx);
   c1->SetName(tmpchr);
   c1->Write();
   sprintf(tmpchr,"p_spectrum_%s.pdf",namesfx);
   c1->Print(tmpchr);
   
   
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
}


int main(int argc, char** argv)
{
  TFile *file = new TFile("output/output_bck.root");
  TFile *fileout = new TFile("output/output_mcshape.root","recreate");

  TH1D *hdata;
  TH1D *hmcmu;
  TH1D *hmcKK;
  gStyle->SetOptStat(0);
  double energys[22] = {
  2.0000, 2.0500, 2.1000, 2.1500, 2.1750,
  2.2000, 2.2324, 2.3094, 2.3864, 2.3960,
  2.5000, 2.6444, 2.6464, 2.7000, 2.8000,
  2.9000, 2.9500, 2.9810, 3.0000, 3.0200,
  3.0800, 2.125
  };
  for (int i=0 ;i<22;i++) {
    char hname[100];
    // get data spectrum
    //sprintf(hname,"hp_KK_%d",(int)(energys[i]*10000+0.5));
    sprintf(hname,"hp_mcKK_%.4f",energys[i]);
    hdata = (TH1D*)file->Get(hname);
  //hdata->SetFillColor(0);
  //hdata->SetTitle("");
  //hdata->GetXaxis()->SetTitle("p (GeV/c)");
  //hdata->GetYaxis()->SetTitle("Counts");
    // get dimu spectrum
    sprintf(hname,"hp_mcmumu_%.4f",energys[i]);
    hmcmu = (TH1D*)file->Get(hname);

    // get mc KK spectrum
    sprintf(hname,"hp_mcKK2_%.4f",energys[i]);
    hmcKK = (TH1D*)file->Get(hname);
    sprintf(hname,"mcshape_%.4f",energys[i]);
    FitSpectrum(hdata, hmcKK,hmcmu,energys[i],hname);
    
  //TCanvas *chp = new TCanvas(hname,"Background analysis");
  //chp->SetMargin(0.15,0.1,0.15,0.1);
  //hdata->SetLineColor(1);
  //hdata->SetTitle("");
  //hdata->GetXaxis()->SetNdivisions(505);
  //hdata->GetXaxis()->SetTitleSize(0.05);
  //hdata->GetXaxis()->SetLabelSize(0.05);
  //hdata->GetYaxis()->SetNdivisions(505);
  //hdata->GetYaxis()->SetTitleSize(0.05);
  //hdata->GetYaxis()->SetLabelSize(0.05);
  //hdata->GetYaxis()->SetTitleOffset(1.1);

  //hdata->Draw("E");
  //// scale dimu to data according luminosity
  //hmcmu->Scale(luminosity[i]/(nevt[i]/mccross[i]));
  //hmcmu->SetLineColor(2);
  //hmcmu->SetMarkerColor(2);
  //hmcmu->SetFillColor(0);
  //hmcmu->Draw("same");
  //
  //// scale KK MC
  //hmcKK->Scale((maxdata-hmcmu->GetBinContent(binid))/maxmcKK);
  //hmcKK->SetLineColor(3);
  //hmcKK->SetMarkerColor(3);
  //hmcKK->SetFillColor(0);
  //hmcKK->Draw("same");

  //sprintf(hname,"%.4f GeV",energys[i]);
  //TLegend* leg = new TLegend(0.4,0.55,0.6,0.85,hname);
  //leg->AddEntry(hdata,"data");
  //leg->AddEntry(hmcmu,"#mu^{+} #mu^{-} MC");
  //leg->AddEntry(hmcKK,"K^{+} K^{-} MC");
  //leg->SetBorderSize(0);
  //leg->SetFillStyle(0);
  //leg->Draw();
  //sprintf(hname,"Canvas_%02d.pdf",i+1);
  //chp->Print(hname);
  //delete chp;
    
  }
  return 0;
}

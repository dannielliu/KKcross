#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
using namespace std;
using namespace RooFit;
double FitSpectrum(TH1D* dataraw, TH1D* hk, TH1D* hmu,  const char* namesfx, TFile* fileout=0);

//int main()
int ComparePolarAng3()
{
  double energys[22] = {
  2.0000, 2.0500, 2.1000, 2.1500, 2.1750,
  2.2000, 2.2324, 2.3094, 2.3864, 2.3960,
  2.5000, 2.6444, 2.6464, 2.7000, 2.8000,
  2.9000, 2.9500, 2.9810, 3.0000, 3.0200,
  3.0800, 2.125
  };
  double MCremain[22] = {
  18847  ,    18100.8,    15652  ,    13047.6,    14736  ,    
  17715.2,    19707.4,    16681.7,    11763.8,    11367.1,    
  9847.6 ,    8913.87,    8839.79,    8554.19,    7974.79,    
  7381.3 ,    6969.49,    6752.07,    6666.27,    6563.12,    
  5702.12,    14136.3 
  };
  double Nevt[22] = {
  1829.22,    523.971,    1418.71,    263.035,    1029.09,    
  1697.6 ,    1613.49,    2080.28,    1269.39,    3841.74,    
  53.3638,    1107.33,    1145.1 ,    22.1   ,    27.9761,    
  1997.55,    265.141,    286.96 ,    223.733,    264.62 ,    
  1417.08,    11088.7    
  };
  double luminosity[22] = {
  10074,   3343,   12167,  2841,   10625,
  13699,  11856,   21089, 22549,   66869,
   1098,  33722,   34003,  1034,    1008,
 105253,  15942,   16071, 15881,   17290,
  126188, 108490
  };
  double mccross[22] = {
  22.3763923, 21.2939032, 20.4018357, 19.4157221, 18.9842334,
  18.5074456, 18.0584889, 16.9121088, 15.8285505, 15.6795155,
  14.3957219, 12.9117968, 12.9257851, 12.3962790, 11.5661025,
  10.7698318, 10.4156161, 10.2122681, 10.0576468, 9.9391192,
  9.5904453,  21.1564523
  };
  double nevt[22] = {
  5.0e5, 5.0e5, 5.0e5, 5.0e5, 5.0e5,
  5.0e5, 5.0e5, 5.0e5, 5.0e5, 1.07e6,
  5.0e5, 5.0e5, 5.0e5, 5.0e5, 5.0e5,
  1.13e6, 5.0e5, 5.0e5, 5.0e5, 5.0e5,
  1.16e6, 1.93e6
  };
  double pexp[22]={
  
  };
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","");
  c1->SetMargin(0.15,0.1,0.15,0.1);
  TFile *file = new TFile("output/output_openang.root");
  char name[1000];
  double mka = 0.493677;
  ofstream scale_info("scale_muon.info");
    scale_info<<"Ecm \tNtot \tNsa \tdNsa \tdNsa/Ntot "<<endl;
  for (int i=0; i<21;i++){
    pexp[i] = sqrt(pow(energys[i]/2,2)-pow(mka,2));
    sprintf(name,"KK_%d",(int)(energys[i]*10000));
    TDirectory* dir_data = (TDirectory*)file->Get(name);
    TTree* var_data = (TTree*) dir_data->Get("vars");
    sprintf(name,"mcKK_ecal_%03d_dst",i+1);
    if (energys[i]>2.124) sprintf(name,"mcKK_ecal_%03d_dst",i+2);
    TDirectory* dir_mckk = (TDirectory*)file->Get(name);
    TTree* var_mckk = (TTree*) dir_mckk->Get("vars");
    sprintf(name,"mcdimu_%03d_%.4f",i+1,energys[i]);
    TDirectory* dir_mcmu = (TDirectory*)file->Get(name);
    TTree* var_mcmu = (TTree*) dir_mcmu->Get("vars");
    TH1D* hang_data = new TH1D("hang_data","",100,-1,1);
    TH1D* hang_mckk = new TH1D("hang_mckk","",100,-1,1);
    TH1D* hang_mcmu = new TH1D("hang_mcmu","",100,-1,1);
    char pcut[100];
  //sprintf(pcut,"theta>179 && p1>%f && p1<%f && p1-%f>0.03",pexp[i]-0.1,pexp[i]+0.2, pexp[i]);
  //if (energys[i]<2.6) 
  //  sprintf(pcut,"theta>179 && p1>%f && p1<%f && p1-%f>0.03",pexp[i]-0.05,pexp[i]+0.15,pexp[i]);
    sprintf(pcut,"theta>179 && p1>%f && p1<%f",pexp[i]-0.1,pexp[i]+0.2);
    if (energys[i]<2.6) 
      sprintf(pcut,"theta>179 && p1>%f && p1<%f",pexp[i]-0.05,pexp[i]+0.15);
    var_data->Draw("costheta1>>hang_data",pcut);
    var_data->Draw("costheta2>>+hang_data",pcut);
    double mcentry1 = var_mckk->Draw("costheta1>>hang_mckk",pcut);
    double mcentry2 = var_mckk->Draw("costheta2>>+hang_mckk",pcut);
    double mcentry = mcentry1+mcentry2;
    var_mcmu->Draw("costheta1>>hang_mcmu",pcut);
    var_mcmu->Draw("costheta2>>+hang_mcmu",pcut);
    
    // statistic the 
    //
    sprintf(pcut,"theta>179 && p1>%f && p1<%f && abs(costheta1)>0.8 && abs(costheta2)>0.8",pexp[i]-0.1,pexp[i]+0.2);
    if (energys[i]<2.6) 
      sprintf(pcut,"theta>179 && p1>%f && p1<%f && abs(costheta1)>0.8 && abs(costheta2)>0.8",pexp[i]-0.05,pexp[i]+0.15);
    double SMu_L = luminosity[i]/(nevt[i]/mccross[i]);
    double Ndata_sa = var_data->Draw("",pcut);
    double Nmcmu_sa = var_mcmu->Draw("",pcut);
    double Nmcmu_sa_s = Nmcmu_sa*SMu_L;
    double fit_scale = Ndata_sa/Nmcmu_sa_s;
    double Nmcmu_s = hang_mcmu->GetEntries()*SMu_L/2;
    scale_info<<energys[i] << "\t "<< Nmcmu_s;
    scale_info<<"\t"<< Nmcmu_sa_s;
    scale_info<<", Fitting will give a scale: "<<fit_scale;
    scale_info<<",\t"<<Nmcmu_sa_s - Ndata_sa;
    scale_info<<",\t"<< (Nmcmu_sa_s - Ndata_sa)/Nmcmu_s<< endl;
    //scale_info<<",\tdN_tot: "<<Nmcmu_s*(1-fit_scale)<<endl;
  //scale_info<<"At "<<energys[i] << "\t Total number of muon is "<< Nmcmu_s;
  //scale_info<<", the number of muon at |costheta|>0.8 : "<< Nmcmu_sa_s;
  ////scale_info<<", Fitting will give a scale: "<<fit_scale;
  //scale_info<<",\tdN_sa: "<<Nmcmu_sa_s - Ndata_sa;
  //scale_info<<",\tdN/Ntot is "<< (Nmcmu_sa_s - Ndata_sa)/Nmcmu_s<< endl;
  ////scale_info<<",\tdN_tot: "<<Nmcmu_s*(1-fit_scale)<<endl;

    hang_data->SetLineColor(kBlue);
    hang_mckk->SetLineColor(kRed);
    hang_mcmu->SetLineColor(kGreen);
    hang_mcmu->Scale(luminosity[i]/(nevt[i]/mccross[i]));
    hang_data->GetXaxis()->SetTitle("cos#theta_{K-}");
    hang_data->GetXaxis()->SetNdivisions(505);
    hang_data->GetXaxis()->SetTitleSize(0.05);
    hang_data->GetXaxis()->SetLabelSize(0.05);
    hang_data->GetYaxis()->SetTitle("count");
    hang_data->GetYaxis()->SetNdivisions(505);
    hang_data->GetYaxis()->SetTitleSize(0.05);
    hang_data->GetYaxis()->SetLabelSize(0.05);
    hang_data->GetYaxis()->SetRangeUser(0, 1.1*hang_data->GetMaximum());
    hang_mckk->Scale(Nevt[i]/mcentry);
    TLegend *leg = new TLegend(0.2,0.7,0.5,0.85);
    sprintf(name,"data at %1.4f GeV",energys[i]);
    leg->AddEntry(hang_data,name);
    leg->AddEntry(hang_mckk,"MC KK");
    leg->AddEntry(hang_mcmu,"MC di-#mu");
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    hang_data->Draw();
    hang_mckk->Draw("same");
    hang_mcmu->Draw("same");
    leg->Draw();
    sprintf(name,"cmppolarangKm_%d.pdf",(int)(energys[i]*10000));
    c1->Print(name);
    
    // fit
  //sprintf(name,"%d",(int)(energys[i]*10000));
  //FitSpectrum(hang_data, hang_mckk, hang_mcmu, name);
  }

  return 0;
}

double FitSpectrum(TH1D* dataraw, TH1D* hk, TH1D* hmu,  const char* namesfx, TFile* fileout)
{
   TCanvas *c1 = new TCanvas("c1_1","",800,600);
   c1->SetMargin(0.15,0.05,0.15,0.05);

   int nBins=100;
   int Npar;

   
   TH1D *hkrebin = (TH1D*)hk->Rebin(1,"hkrebin");
   TH1D *hmurebin = (TH1D*)hmu->Rebin(1,"hmurebin");
  //
   //fileout->cd();
   //
   // try to use roofit
   //
   RooRealVar x("x","polar angle",-1.0,1.0);
    // signal
 
   RooDataHist mckhist("mcKhist","mc K spec",x,hkrebin);
   RooDataHist mcmuhist("mcmuhist","mc mu spec",x,hmurebin);
 //RooHistPdf sig("sig","signal",x,mckhist,4); 
   RooHistPdf kpdf("kpdf","signal",x,mckhist,4); 
   RooHistPdf mupdf("mupdf","background",x,mcmuhist,4);

   RooRealVar mean("mean","mean of gaussian",0.0,-0.01, 0.01);
   RooRealVar sigma("sigma","width of gaussian",0.001,0.00,0.1);
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
    
   xframe->SetTitle("");
   xframe->GetXaxis()->SetTitle("polar angle");
   xframe->GetXaxis()->SetNdivisions(505);
   xframe->GetXaxis()->SetTitleSize(0.05);
   xframe->GetXaxis()->SetLabelSize(0.05);
   xframe->GetXaxis()->SetTitleOffset(1.1);
   xframe->GetYaxis()->SetNdivisions(505);
   xframe->GetYaxis()->SetTitleSize(0.05);
   xframe->GetYaxis()->SetTitleOffset(1.25);
   xframe->GetYaxis()->SetLabelSize(0.05);
   xframe->GetYaxis()->SetNoExponent(0);
   
   xframe->Draw();
   TPaveText *pt = new TPaveText(0.17,0.70,0.45,0.93,"BRNDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(4000);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.05);
   sprintf(tmpchr,"@ %1.4f GeV", atof(namesfx)/(double)10000);
   pt->AddText(tmpchr);
 //sprintf(tmpchr,"#mu = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
 //pt->AddText(tmpchr);
 //sprintf(tmpchr,"#sigma = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
 //pt->AddText(tmpchr);
   sprintf(tmpchr,"N_{sig} = %.2f #pm %.2f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"N_{bck} = %.2f #pm %.2f",background.getVal(),background.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
   sprintf(tmpchr,"polarAng_%s",namesfx);
   //c1->SetName(tmpchr);
   //c1->Write();
   sprintf(tmpchr,"polarAng_%s.pdf",namesfx);
   c1->Print(tmpchr);
   
   
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
}



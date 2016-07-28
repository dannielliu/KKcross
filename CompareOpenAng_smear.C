#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TRandom.h"
#include <iostream>
#include <fstream>
using namespace std;

void help(){
  cout<<"help"<<endl;
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
using namespace RooFit;


#include "RooHistPdf.h"
#include "RooFFTConvPdf.h"
double FitSpectrum(TTree *dataraw, TH1D* hsig, TH1D* hbck, double beame, char* namesfx, TFile* fileout)
{
   TCanvas *c2 = new TCanvas("","",800,600);
   
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double mmu = 0.1057;
   double Ecm = beame;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   double pmu = sqrt(pow(beame,2)-pow(mmu,2));
   std::cout<<"Ecm is "<<Ecm<<", Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;
   double beamlow=177;
   double beamup=180;
   
   //
   fileout->cd();
   //
   // try to use roofit
   //
   RooRealVar x("theta","theta",beamlow,beamup);
    // signal
   cout<<"aaaaaaaaa"<<endl;
   RooDataHist mckhist("mcKhist","mc K spec",x,hsig);
   RooDataHist mcmuhist("mcmuhist","mc mu spec",x,hbck);
   cout<<"aaaaaaaaa"<<endl;
   RooHistPdf kpdf("kpdf","signal",x,mckhist,4); 
   RooHistPdf mupdf("mupdf","background",x,mcmuhist,4);
   cout<<"aaaaaaaaa"<<endl;

   RooRealVar mean("mean","mean of gaussian",-0.001 ,-0.2, 0.000);
   RooRealVar sigma("sigma","width of gaussian",0.01,0.0005,0.2);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);

   // signal
   RooFFTConvPdf sig("sig","signal",x,kpdf,gaus); 
   //RooHistPdf sig("sig","signal",x,mckhist,2); 
   // background, mainly di-mu
   RooFFTConvPdf bck("bck","background",x,mupdf,gaus); 
   //RooHistPdf bck("bck","background",x,mcmuhist,2);
   
   RooRealVar signal("signal"," ",1000,0,1000000);//event number
   RooRealVar background("background"," ",200,0,10000000);
 
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
 //sprintf(tmpchr,"ang_%s",namesfx);
 //c1->SetName(tmpchr);
 //c1->Write();
   sprintf(tmpchr,"ang_%s.pdf",namesfx);
   c2->Print(tmpchr);
   
   ofstream otxt("ang_smear.txt",ios::app);
   otxt<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<"\t"<<sigma.getVal()<<"\t"<<sigma.getError()<<endl;



   delete xframe;
   //delete dataset;
   delete sum;
   return signal.getVal();
}

//int CompareOpenAng_smear()
int main()
{
  help();

  double energys[22] = {
  2.0000, 2.0500, 2.1000, 2.125, 2.1500, 2.1750,
  2.2000, 2.2324, 2.3094, 2.3864, 2.3960,
  2.5000, 2.6444, 2.6464, 2.7000, 2.8000,
  2.9000, 2.9500, 2.9810, 3.0000, 3.0200,
  3.0800  };
  double MCremain[22] = {
  18847  ,    18100.8,    15652  ,    14136.3,    13047.6,    14736  ,    
  17715.2,    19707.4,    16681.7,    11763.8,    11367.1,    
  9847.6 ,    8913.87,    8839.79,    8554.19,    7974.79,    
  7381.3 ,    6969.49,    6752.07,    6666.27,    6563.12,    
  5702.12 
  };
  double Nkkend[22]={
  18847, 18100.8, 15652, 14136.3, 13047.6, 14736,
  17715.2, 19707.4, 16681.7, 11763.8, 11367.1,
  9847.6, 8913.87, 8839.79, 8554.19, 7974.79,
  7381.3, 6969.49, 6752.07, 6666.27, 6563.12,
  5702.12
  };
  double Nevt[22] = {
  1829.22,    523.971,    1418.71,    11088.7,    263.035,    1029.09,    
  1697.6 ,    1613.49,    2080.28,    1269.39,    3841.74,    
  53.3638,    1107.33,    1145.1 ,    22.1   ,    27.9761,    
  1997.55,    265.141,    286.96 ,    223.733,    264.62 ,    
  1417.08    
  };
  double Nmuinievt[22] = {
  5.0e5, 5.0e5, 5.0e5, 1.93e6, 5.0e5, 5.0e5,
  5.0e5, 5.0e5, 5.0e5, 5.0e5, 1.07e6,
  5.0e5, 5.0e5, 5.0e5, 5.0e5, 5.0e5,
  1.13e6, 5.0e5, 5.0e5, 5.0e5, 5.0e5,
  1.16e6
  };

  double luminosity[22] = {
  10074,   3343,   12167, 108490,  2841,   10625,
  13699,  11856,   21089, 22549,   66869,
   1098,  33722,   34003,  1034,    1008,
 105253,  15942,   16071, 15881,   17290,
  126188
  };
  double mcmucross[22] = {
  22.3763923, 21.2939032, 20.4018357, 19.90877 , 19.4157221, 18.9842334,
  18.5074456, 18.0584889, 16.9121088, 15.8285505, 15.6795155,
  14.3957219, 12.9117968, 12.9257851, 12.3962790, 11.5661025,
  10.7698318, 10.4156161, 10.2122681, 10.0576468, 9.9391192,
  9.5904453
  };
  double MCobsxs[22]={
  0.953909, 0.871322, 0.764996, 0.715793, 0.675351, 0.67103,  
  0.685954, 0.690814, 0.595784, 0.493849, 0.48468,
  0.414689, 0.351078, 0.35032,  0.33038,  0.297999,
  0.271363, 0.258785, 0.251148, 0.246597, 0.241902,
  0.227902
  };
  double pexp[22]={
  
  };
  double pcut[22]={
  0.883571,   0.912768,   0.94195 ,   0.956409,   0.970868,   0.985444,   
  0.999846,   1.0183  ,   1.06241 ,   1.1062  ,   1.11176 ,   
  1.16991 ,   1.25103 ,   1.25245 ,   1.28212 ,   1.33714 ,   
  1.39212 ,   1.4198  ,   1.43688 ,   1.44717 ,   1.45769 ,   
  1.49096  
  };
  double psigma[22];

  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","");
  TPad *pad1 = new TPad("pad1","pad1",0.05,0.3,1.0,0.98);
  TPad *pad2 = new TPad("pad2","pad2",0.05,0.02,1.0,0.3);
  pad1->SetTopMargin(0.05);
  pad1->SetBottomMargin(0.02);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.3);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  char name[1000];
  double mka = 0.493677;
  for (int i=0;i<22;i++){
    pexp[i] = sqrt(pow(energys[i]/2,2)-pow(mka,2));
    psigma[i] = (pcut[i]-pexp[i])/3;
    //cout << "Ecm = "<< energys[i]<<", pexp = "<< pexp[i]<<", sigma = "<< (pcut[i]-pexp[i])/3<<endl;
  }
 
  TFile *fileout = new TFile("output/output_noF_1.root","recreate");
  double totp, costheta1, costheta2;
  double theta;
  double p1,p2;
  TTree *elist = new TTree("elist","elist");
  elist->Branch("theta",&theta,"theta/D");
  
  TFile *file = new TFile("output/output_openang.root");
  ofstream of("openanglediff.txt");
  of<< "Energy\t Ndata_179 \t Kscale_1795 \t diff "<<endl;
  for (int i=0; i<22;i++){
    sprintf(name,"KK_%d",(int)(energys[i]*10000));
    TDirectory* dir_data = (TDirectory*)file->Get(name);
    TTree* var_data = (TTree*) dir_data->Get("vars");
    var_data->SetBranchAddress("p1",&p1);
    var_data->SetBranchAddress("p2",&p2);
    var_data->SetBranchAddress("theta",&theta);

    sprintf(name,"mcKK_ecal_%03d_dst",i+1);
    TDirectory* dir_mckk = (TDirectory*)file->Get(name);
    TTree* var_mckk = (TTree*) dir_mckk->Get("vars");
    var_mckk->SetBranchAddress("p1",&p1);
    var_mckk->SetBranchAddress("p2",&p2);
    var_mckk->SetBranchAddress("theta",&theta);
    sprintf(name,"mcdimu_%03d_%.4f",i+1,energys[i]);
    if (i==3) sprintf(name,"mcdimu_%.4f",energys[i]);
    if (i>3) sprintf(name,"mcdimu_%03d_%.4f",i,energys[i]);
    TDirectory* dir_mcmu = (TDirectory*)file->Get(name);
    TTree* var_mcmu = (TTree*) dir_mcmu->Get("vars");
  //sprintf(name,"bckhadron_%03d_%.4f",i+1,energys[i]);
  //TDirectory* dir_mchad = (TDirectory*)file->Get(name);
  //TTree* var_mchad = (TTree*) dir_mchad->Get("vars");
    
    //cout<<"pointers: data "<<var_data<<", mckk "<<var_mckk<<", mcmu "<<var_mcmu<<endl;
    pad1->cd();

    //double pexp = sqrt(pow(energys[i]/2,2)-pow(mka,2));
    char evt_flt[1000];
    //sprintf(evt_flt,"theta>177 && theta<180 && p1>%f && p1<%f",pexp-0.03,pexp+0.03);
    sprintf(evt_flt,"theta>177 && theta<180 && p1>%f && p1<%f && p2<%f ",pexp[i]-0.1,pexp[i]+0.2,pcut[i]);
    if (energys[i]<2.6) 
      sprintf(evt_flt,"theta>177 && theta<180 && p1>%f && p1<%f && p2<%f ",pexp[i]-0.05,pexp[i]+0.15,pcut[i]);
    TH1D* hang_data = new TH1D("hang_data","",100,177,180);
    TH1D* hang_mckk = new TH1D("hang_mckk","",100,177,180);
    TH1D* hang_mcmu = new TH1D("hang_mcmu","",100,177,180);
    TH1D* hang_mchad = new TH1D("hang_mchad","",100,177,180);
    TH1D* hang_mcsum = new TH1D("hang_mcsum","",100,177,180);
    var_data->Draw("theta>>hang_data",evt_flt);
    double mcentry = var_mckk->Draw("theta>>hang_mckk",evt_flt);
    var_mcmu->Draw("theta>>hang_mcmu",evt_flt);
    //var_mchad->Draw("theta>>hang_mchad",evt_flt);
    double entries = var_mckk->GetEntries();

    // cout before smear and after smear
  //double cnt_bf=0;
  //double cnt_af=0;
  //for (int ien=0; ien<entries; ien++){
  //  var_mckk->GetEntry(ien);
  //  if (theta>179 && theta<180 && p2<pcut[i]) {
  //    if (energys[i]>2.6 && p1>pexp[i]-0.1 && p1<pexp[i]+0.2 ) cnt_bf++;
  //    if (energys[i]<2.6 && p1>pexp[i]-0.05 && p1<pexp[i]+0.15 ) cnt_bf++;
  //  }

  //  theta = theta+gRandom->Gaus(0.0085972-0.00652*energys[i],0.04);
  //  if (theta>179 && theta<180 && p2<pcut[i]) {
  //    if (energys[i]>2.6 && p1>pexp[i]-0.1 && p1<pexp[i]+0.2 ) cnt_af++;
  //    if (energys[i]<2.6 && p1>pexp[i]-0.05 && p1<pexp[i]+0.15 ) cnt_af++;
  //  }
  //}
  //cout<<energys[i]<<"\t"<<cnt_bf<<"\t"<<cnt_af<<endl;
  //continue;


    fileout->cd();
    elist->Reset();
    for (int ien=0; ien<entries; ien++){
      var_data->GetEntry(ien);
      if (theta>177 && theta<180 && p2<pcut[i]) {
      if (energys[i]>2.6 && p1>pexp[i]-0.1 && p1<pexp[i]+0.2 ) elist->Fill();
      if (energys[i]<2.6 && p1>pexp[i]-0.05 && p1<pexp[i]+0.15 ) elist->Fill();
      }
    }
    cout<<"size of elist "<<elist->GetEntries()<<endl;



    // count number
    double muscale_l = luminosity[i]/(Nmuinievt[i]/mcmucross[i]);
    hang_data->SetLineColor(2);
    hang_mckk->SetLineColor(3);
    hang_mcmu->SetLineColor(4);
    //hang_mchad->SetLineColor(6);
    hang_mcsum->SetLineColor(7);
    hang_data->SetLineWidth(2);
    hang_mckk->SetLineWidth(2);
    hang_mcmu->SetLineWidth(2);
    //hang_mchad->SetLineWidth(2);
    hang_mcsum->SetLineWidth(2);
    hang_data->GetXaxis()->SetTitle("angle");
    hang_data->GetXaxis()->SetNdivisions(505);
    hang_data->GetXaxis()->SetTitleSize(0.05);
    hang_data->GetXaxis()->SetLabelSize(0.05);
    hang_data->GetXaxis()->SetLabelOffset(1.5);
    hang_data->GetYaxis()->SetTitle("count");
    hang_data->GetYaxis()->SetNdivisions(505);
    hang_data->GetYaxis()->SetTitleSize(0.05);
    hang_data->GetYaxis()->SetLabelSize(0.05);
    hang_data->GetYaxis()->SetRangeUser(0, 1.1*hang_data->GetMaximum());
    hang_mckk->Scale(Nevt[i]/Nkkend[i]);
    hang_mcmu->Scale(luminosity[i]/(Nmuinievt[i]/mcmucross[i]));
    //hang_mchad->Scale(hang_mckk->GetMaximum()/hang_mchad->GetMaximum());
    hang_mcsum->Add(hang_mckk);
    //hang_mcsum->Add(hang_mchad);
    hang_mcsum->Add(hang_mcmu);
    TLegend *leg = new TLegend(0.3,0.6,0.6,0.88);
    leg->AddEntry(hang_data,"data");
    leg->AddEntry(hang_mckk,"MC KK");
    leg->AddEntry(hang_mcmu,"MC #mu#mu");
    //leg->AddEntry(hang_mchad,"MC hadron");
    leg->AddEntry(hang_mcsum,"MC sum");
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    hang_data->Draw("E");
    hang_mckk->Draw("same");
    hang_mcmu->Draw("same");
    //hang_mchad->Draw("same");
    hang_mcsum->Draw("same");
    leg->Draw();

    pad2->cd();

    TH1D* hang_diff = new TH1D("hang_diff","",100,177,180);
    for (int ibin=1; ibin<=100; ibin++){
      double ndata = hang_data->GetBinContent(ibin);
      if (fabs(ndata)<1) ndata=1; 
      double bindif = (hang_mcsum->GetBinContent(ibin) - hang_data->GetBinContent(ibin))/ndata;
      double binerr = sqrt(pow(hang_mcsum->GetBinError(ibin)/ndata,2)
                          +pow(hang_data->GetBinError(ibin)/ndata,2))*(1-bindif);
      hang_diff->SetBinContent(ibin,bindif);
      hang_diff->SetBinError(ibin,binerr);
    }
    hang_diff->GetXaxis()->SetTitle("angle");
    hang_diff->GetYaxis()->SetTitle("#Delta_{relative}");
    //hang_diff->GetYaxis()->SetTitle("#frac{#epsilon_{data} - #epsilon_{MC}}{#epsilon_{data}}");
    hang_diff->GetXaxis()->SetNdivisions(505);
    hang_diff->GetYaxis()->SetNdivisions(502);
    hang_diff->GetXaxis()->SetTitleSize(0.12);
    hang_diff->GetXaxis()->SetLabelSize(0.12);
    hang_diff->GetYaxis()->SetTitleSize(0.15);
    hang_diff->GetYaxis()->SetLabelSize(0.12);
    hang_diff->GetYaxis()->SetRangeUser(-1,1);
    hang_diff->GetYaxis()->SetTitleOffset(0.3);
    hang_diff->Draw();
    
    TF1 f1("f1","0",177,180);
    f1.SetLineStyle(kDashed);
    f1.Draw("same");
    c1->cd();
    //pad1->Draw();
    //pad2->Draw();
    
    sprintf(name,"cmpang_%d.pdf",(int)(energys[i]*10000));
    c1->Print(name);

  //cout<<"entries: "<<elist->GetEntries()<<" "<<hang_mckk->GetEntries()<<" "<< hang_mcmu->GetEntries()<<endl;
  //sprintf(name,"KKang_%d",(int)(energys[i]*10000));
  //FitSpectrum(elist, hang_mckk, hang_mcmu, energys[i],name,fileout);
  }

}



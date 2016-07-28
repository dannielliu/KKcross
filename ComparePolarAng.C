int ComparePolarAng()
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
  double Nkkend[22]={
  18847, 18100.8, 15652, 13047.6, 14736,
  17715.2, 19707.4, 16681.7, 11763.8, 11367.1,
  9847.6, 8913.87, 8839.79, 8554.19, 7974.79,
  7381.3, 6969.49, 6752.07, 6666.27, 6563.12,
  5702.12, 14136.3
  };
  double Nevt[22] = {
  1829.22,    523.971,    1418.71,    263.035,    1029.09,    
  1697.6 ,    1613.49,    2080.28,    1269.39,    3841.74,    
  53.3638,    1107.33,    1145.1 ,    22.1   ,    27.9761,    
  1997.55,    265.141,    286.96 ,    223.733,    264.62 ,    
  1417.08,    11088.7    
  };
  double Nmuinievt[22] = {
  5.0e5, 5.0e5, 5.0e5, 5.0e5, 5.0e5,
  5.0e5, 5.0e5, 5.0e5, 5.0e5, 1.07e6,
  5.0e5, 5.0e5, 5.0e5, 5.0e5, 5.0e5,
  1.13e6, 5.0e5, 5.0e5, 5.0e5, 5.0e5,
  1.16e6, 1.93e6
  };
  double luminosity[22] = {
  10074,   3343,   12167,  2841,   10625,
  13699,  11856,   21089, 22549,   66869,
   1098,  33722,   34003,  1034,    1008,
 105253,  15942,   16071, 15881,   17290,
  126188, 108490
  };
  double mcmucross[22] = {
  22.3763923, 21.2939032, 20.4018357, 19.4157221, 18.9842334,
  18.5074456, 18.0584889, 16.9121088, 15.8285505, 15.6795155,
  14.3957219, 12.9117968, 12.9257851, 12.3962790, 11.5661025,
  10.7698318, 10.4156161, 10.2122681, 10.0576468, 9.9391192,
  9.5904453,  21.1564523
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
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","");
  c1->SetMargin(0.15,0.1,0.15,0.1);
  TFile *file = new TFile("output/output_openang.root");
  double mka = 0.493677;
  char name[1000];
  for (int i=0; i<21;i++){
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
    TH1D* hang_data = new TH1D("hang_data","",50,-1,1);
    TH1D* hang_mckk = new TH1D("hang_mckk","",50,-1,1);
    TH1D* hang_mcmu = new TH1D("hang_mcmu","",50,-1,1);
    double pexp = sqrt(pow(energys[i]/2,2)-pow(mka,2));
    char evtflt[100];
    sprintf(evtflt,"theta>179 && p1>%f",pexp+0.03);
    var_data->Draw("costheta2>>hang_data",evtflt);
    double mcentry = var_mckk->Draw("costheta2>>hang_mckk",evtflt);
    var_mcmu->Draw("costheta2>>hang_mcmu",evtflt);
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
    //hang_mckk->Scale(Nevt[i]/mcentry);
    hang_mckk->Scale(Nevt[i]/Nkkend[i]);
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
    sprintf(name,"cmppolarang_%d_Km.pdf",(int)(energys[i]*10000));
    c1->Print(name);
  }

}

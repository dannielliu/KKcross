int CompareOpenCos()
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

  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","");
  c1->SetMargin(0.15,0.1,0.15,0.1);
  TFile *file = new TFile("output/output_openang.root");
  char name[1000];
  for (int i=0; i<21;i++){
    sprintf(name,"KK_%d",(int)(energys[i]*10000));
    TDirectory* dir_data = (TDirectory*)file->Get(name);
    TTree* var_data = (TTree*) dir_data->Get("vars");
    sprintf(name,"mcKK_ecal_%03d_dst",i+1);
    if (energys[i]>2.124) sprintf(name,"mcKK_ecal_%03d_dst",i+2);
    TDirectory* dir_mckk = (TDirectory*)file->Get(name);
    TTree* var_mckk = (TTree*) dir_mckk->Get("vars");
    TH1D* hang_data = new TH1D("hang_data","",100,-1,-0.999);
    TH1D* hang_mckk = new TH1D("hang_mckk","",100,-1,-0.999);
    var_data->Draw("costheta>>hang_data","costheta<-0.999");
    double mcentry = var_mckk->Draw("costheta>>hang_mckk","costheta<-0.999");
    hang_data->SetLineColor(kBlue);
    hang_mckk->SetLineColor(kRed);
    hang_data->GetXaxis()->SetTitle("cos#theta");
    hang_data->GetXaxis()->SetNdivisions(505);
    hang_data->GetXaxis()->SetTitleSize(0.05);
    hang_data->GetXaxis()->SetLabelSize(0.05);
    hang_data->GetYaxis()->SetTitle("count");
    hang_data->GetYaxis()->SetNdivisions(505);
    hang_data->GetYaxis()->SetTitleSize(0.05);
    hang_data->GetYaxis()->SetLabelSize(0.05);
    hang_data->GetYaxis()->SetRangeUser(0, 1.1*hang_data->GetMaximum());
    hang_mckk->Scale(Nevt[i]/mcentry);
    TLegend *leg = new TLegend(0.3,0.6,0.6,0.75);
    leg->AddEntry(hang_data,"data");
    leg->AddEntry(hang_mckk,"MC");
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    hang_data->Draw();
    hang_mckk->Draw("same");
    leg->Draw();
    sprintf(name,"cmpcos_%d.pdf",(int)(energys[i]*10000));
    c1->Print(name);
  }

}

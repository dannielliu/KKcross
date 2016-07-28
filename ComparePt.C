int ComparePt()
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
  double pcut[22] = {
  0.883571	 ,    0.912768  ,    0.94195 ,    0.970868 ,    0.985444 ,
  0.999846	 ,    1.0183    ,    1.06241 ,    1.1062   ,    1.11176  ,
  1.16991	 ,    1.25103   ,    1.25245 ,    1.28212  ,    1.33714  ,
  1.39212	 ,    1.4198    ,    1.43688 ,    1.44717  ,    1.45769  ,
  1.49096	 ,    0.956409 
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
    
    TH1D* hang_datap = new TH1D("hang_datap","",100,0,1.6);
    TH1D* hang_datam = new TH1D("hang_datam","",100,0,1.6);
    TH1D* hang_mckkp = new TH1D("hang_mckkp","",100,0,1.6);
    TH1D* hang_mckkm = new TH1D("hang_mckkm","",100,0,1.6);
    sprintf(name,"theta>179 && p1 < %f && p1> %f && p2 > %f", pcut[i], pcut[i]-0.2, pcut[i]-0.2);
    //sprintf(name,"theta>179 && p1 < %f", pcut[i]);
    var_data->Draw("pt1>>hang_datap", name);
    var_data->Draw("pt2>>hang_datam", name);
    double mcentry = var_mckk->Draw("pt1>>hang_mckkp", name);
    var_mckk->Draw("pt2>>hang_mckkm", name);
    hang_datap->SetLineColor(2);
    hang_datam->SetLineColor(3);
    hang_mckkp->SetLineColor(4);
    hang_mckkm->SetLineColor(6);
    hang_mckkp->Scale(Nevt[i]/mcentry);
    hang_mckkm->Scale(Nevt[i]/mcentry);
    hang_datap->GetXaxis()->SetTitle("p_{t} (GeV/c)");
    hang_datap->GetXaxis()->SetNdivisions(505);
    hang_datap->GetXaxis()->SetTitleSize(0.05);
    hang_datap->GetXaxis()->SetLabelSize(0.05);
    hang_datap->GetYaxis()->SetTitle("count");
    hang_datap->GetYaxis()->SetNdivisions(505);
    hang_datap->GetYaxis()->SetTitleSize(0.05);
    hang_datap->GetYaxis()->SetLabelSize(0.05);
    hang_datap->GetYaxis()->SetRangeUser(0, 1.2*hang_datap->GetMaximum());
    //hang_mckk->Scale(Nevt[i]/mcentry);
    TPaveText *pt = new TPaveText(0.15,0.9*hang_datap->GetMaximum(),0.6,hang_datap->GetMaximum());
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    sprintf(name, "@ %.4f GeV", energys[i]);
    pt->AddText(name);
    TLegend *leg = new TLegend(0.2,0.55,0.5,0.8);
    sprintf(name,"data Pt_{K^{+}}",energys[i]);
    leg->AddEntry(hang_datap,name);
    sprintf(name,"data Pt_{K^{-}}",energys[i]);
    leg->AddEntry(hang_datam,name);
    sprintf(name,"MC Pt_{K^{+}}",energys[i]);
    leg->AddEntry(hang_mckkp,name);
    sprintf(name,"MC Pt_{K^{-}}",energys[i]);
    leg->AddEntry(hang_mckkm,name);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    hang_datap->Draw();
    hang_datam->Draw("same");
    hang_mckkp->Draw("same");
    hang_mckkm->Draw("same");
    pt->Draw();
    leg->Draw();
    sprintf(name,"pt_%d.pdf",(int)(energys[i]*10000));
    c1->Print(name);
  }

}

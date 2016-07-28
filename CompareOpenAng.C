
int CompareOpenAng()
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
  double MCobsxs[22]={
  0.953909, 0.871322, 0.764996, 0.675351, 0.67103,  
  0.685954, 0.690814, 0.595784, 0.493849, 0.48468,
  0.414689, 0.351078, 0.35032,  0.33038,  0.297999,
  0.271363, 0.258785, 0.251148, 0.246597, 0.241902,
  0.227902, 0.715793
  };
  double pexp[22]={
  
  };
  double pcut[22]={
  0.883571,   0.912768,   0.94195 ,   0.970868,   0.985444,   
  0.999846,   1.0183  ,   1.06241 ,   1.1062  ,   1.11176 ,   
  1.16991 ,   1.25103 ,   1.25245 ,   1.28212 ,   1.33714 ,   
  1.39212 ,   1.4198  ,   1.43688 ,   1.44717 ,   1.45769 ,   
  1.49096 ,   0.956409 
  };
  double psigma[21];

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
  TFile *file = new TFile("output/output_openang.root");
  char name[1000];
  double mka = 0.493677;
  for (int i=0;i<21;i++){
    pexp[i] = sqrt(pow(energys[i]/2,2)-pow(mka,2));
    psigma[i] = (pcut[i]-pexp[i])/3;
    //cout << "Ecm = "<< energys[i]<<", pexp = "<< pexp[i]<<", sigma = "<< (pcut[i]-pexp[i])/3<<endl;
  }
 
  ofstream of("openanglediff.txt");
  of<< "Energy\t Ndata_179 \t Kscale_1795 \t diff "<<endl;
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
    sprintf(name,"bckhadron_%03d_%.4f",i+1,energys[i]);
    TDirectory* dir_mchad = (TDirectory*)file->Get(name);
    TTree* var_mchad = (TTree*) dir_mchad->Get("vars");
    
    cout<<"pointers: data "<<var_data<<", mckk "<<var_mckk<<", mcmu "<<var_mcmu<<endl;
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
    var_mchad->Draw("theta>>hang_mchad",evt_flt);
    
    // count number
    double muscale_l = luminosity[i]/(Nmuinievt[i]/mcmucross[i]);
    sprintf(evt_flt,"theta>179 && theta<180 && p1>%f && p1<%f && p2<%f ",pexp[i]-0.1,pexp[i]+0.2,pcut[i]);
    if (energys[i]<2.6)
      sprintf(evt_flt,"theta>179 && theta<180 && p1>%f && p1<%f && p2<%f ",pexp[i]-0.05,pexp[i]+0.15,pcut[i]);
    double Ndata_179 = var_data->Draw("",evt_flt);
    double Nmcmu_179 = var_mcmu->Draw("",evt_flt);
    double Nmcka_179 = var_mckk->Draw("",evt_flt);
    double Ndata_k179 = Ndata_179 - Nmcmu_179*muscale_l;
    sprintf(evt_flt,"theta>179.5 && theta<180 && p1>%f && p1<%f && p2<%f",pexp[i]-0.1,pexp[i]+0.2,pcut[i]);
    if (energys[i]<2.6)
      sprintf(evt_flt,"theta>179.5 && theta<180 && p1>%f && p1<%f && p2<%f",pexp[i]-0.05,pexp[i]+0.15,pcut[i]);
    double Ndata_1795 = var_data->Draw("",evt_flt);
    double Nmcmu_1795 = var_mcmu->Draw("",evt_flt);
    double Nmcka_1795 = var_mckk->Draw("",evt_flt);
    double Kscale_1795 = (Ndata_1795-Nmcmu_1795*muscale_l)/Nmcka_1795;
    double Nk_scalewith1795 = Nmcka_179*Kscale_1795;
    cout<<"test"<<endl;
    //cout<< "At "<< energys[i]<<",\t Ndata_179 "<<Ndata_k179 <<",\t Kscale_1795 "<< Kscale_1795 <<",\t diff "<< Ndata_k179-Nmcka_179*Kscale_1795<<endl;
    cout<< energys[i]<< "\t "<<Nevt[i] << ",\t "<< Ndata_k179-Nmcka_179*Kscale_1795 << "\t"<<(Ndata_k179-Nmcka_179*Kscale_1795)/Nevt[i]<<endl;
    //of<< energys[i]<< "\t "<<Nevt[i] << ",\t "<< Ndata_k179-Nmcka_179*Kscale_1795 << "\t"<<(Ndata_k179-Nmcka_179*Kscale_1795)/Nevt[i]<<endl;
    cout<<"test2"<<endl;
    //return 0;
    hang_data->SetLineColor(2);
    hang_mckk->SetLineColor(3);
    hang_mcmu->SetLineColor(4);
    hang_mchad->SetLineColor(6);
    hang_mcsum->SetLineColor(7);
    hang_data->SetLineWidth(2);
    hang_mckk->SetLineWidth(2);
    hang_mcmu->SetLineWidth(2);
    hang_mchad->SetLineWidth(2);
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
    hang_mchad->Scale(hang_mckk->GetMaximum()/hang_mchad->GetMaximum());
    hang_mcsum->Add(hang_mckk);
    //hang_mcsum->Add(hang_mchad);
    hang_mcsum->Add(hang_mcmu);
    TLegend *leg = new TLegend(0.3,0.6,0.6,0.88);
    leg->AddEntry(hang_data,"data");
    leg->AddEntry(hang_mckk,"MC KK");
    leg->AddEntry(hang_mcmu,"MC #mu#mu");
    leg->AddEntry(hang_mchad,"MC hadron");
    leg->AddEntry(hang_mcsum,"MC sum");
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    hang_data->Draw("E");
    hang_mckk->Draw("same");
    hang_mcmu->Draw("same");
    hang_mchad->Draw("same");
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
  }

}

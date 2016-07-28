double mumuxs(double *_Ecm, double *par)
{
  double Ecm = _Ecm[0];
  double E = Ecm/2.;
  double mmu = 0.105658;
  double alpha = 1/137.;
  double GeV2Tombar = 0.38938;
  double GeV2Tonbar = GeV2Tombar*1e6;
  double xs = 4*TMath::Pi()*pow(alpha,2)/(3*pow(Ecm,2))*(1-3./8*pow(mmu/E,4)) *GeV2Tonbar;
  return xs;
}

double mumuxsang(double *ang, double *_Ecm)
{
  double Ecm = _Ecm[0];
  double E = Ecm/2.;
  double mmu = 0.105658;
  double alpha = 1/137.;
  double pi = TMath::Pi();
  double GeV2Tombar = 0.38938;
  double GeV2Tonbar = GeV2Tombar*1e6;
  double xs = pow(alpha,2)/(4*pow(Ecm,2))*(1+pow(cos(ang[0]),2))*2*pi*sin(ang[0]) *GeV2Tonbar;
  return xs;
}

double mumuxsdang(double *_Ecm, double *_ang)
{
  double Ecm = _Ecm[0];
  TF1 f1("dimuxs",mumuxsang,0,TMath::Pi(),1);
  f1.SetParameter(0,Ecm);
  double xs = f1.Integral(20./180.*TMath::Pi(),160./180.*TMath::Pi());
  return xs;
}

double bhabhaxsang(double *ang, double *_Ecm)
{
  double Ecm = _Ecm[0];
  double alpha = 1/137.;
  double pi = TMath::Pi();
  double GeV2Tombar = 0.38938;
  double GeV2Tonbar = GeV2Tombar*1e6;
  double xs = pow(alpha,2)/(16*pow(Ecm,2))*pow(3+pow(cos(ang[0]),2),2)/pow(sin(ang[0]/2),4) *2*pi*sin(ang[0]) *GeV2Tonbar;
  return xs;
}

double bhabhaxsdang(double *_Ecm, double *_ang)
{
  double Ecm = _Ecm[0];
  TF1 f1("BBxs",bhabhaxsang,0,TMath::Pi(),1);
  f1.SetParameter(0,Ecm);
  double xs = f1.Integral(20./180.*TMath::Pi(),160./180.*TMath::Pi());
  return xs;
}

double getMax(TH1D* h, int i_low, int i_up, int &binid)
{
  double max=0;
  for (int i=i_low; i<=i_up;i++){ 
    double v = h->GetBinContent(i);
    if (v>max) {max=v;binid=i;}
  }
  return max;
}



int mumucross()
{
  TCanvas *c1 = new TCanvas("c1","di-mu cross section");
  TF1 *f1 = new TF1("f1",mumuxsdang,2,3.1,0);
  f1->GetXaxis()->SetTitle("Ecm (GeV)");
  f1->GetYaxis()->SetTitle("#sigma (nbar)");
  f1->Draw();
  
  TCanvas *c2 = new TCanvas("c2","BhaBha cross section");
  TF1 *f2 = new TF1("f2",bhabhaxsdang,2,3.1,0);
  f2->GetXaxis()->SetTitle("Ecm (GeV)");
  f2->GetYaxis()->SetTitle("#sigma (nbar)");
  f2->Draw();
  //return 0;

  double energys[22] = {
  2.0000, 2.0500, 2.1000, 2.1500, 2.1750,
  2.2000, 2.2324, 2.3094, 2.3864, 2.3960,
  2.5000, 2.6444, 2.6464, 2.7000, 2.8000,
  2.9000, 2.9500, 2.9810, 3.0000, 3.0200,
  3.0800, 2.125
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

  TFile *file = new TFile("output/output_bck.root");

  TH1D *hdata;
  TH1D *hmcmu;
  TH1D *hmcKK;
  gStyle->SetOptStat(0);
  for (int i=0 ;i<22;i++) {
    //double nmumu = luminosity[i]*f1->Eval(energys[i]);
    //cout<<nmumu<<"\t events at "<< energys[i] <<" GeV."<<endl;
    //double nee = luminosity[i]*f2->Eval(energys[i]);
    //cout<<nee<<"\t events at "<< energys[i] <<" GeV."<<endl;
    cout<<"ene "<<energys[i]<<", scale is "<< luminosity[i]/(nevt[i]/mccross[i])<<endl;
    //continue;
    char hname[100];
    // get data spectrum
    //sprintf(hname,"hp_KK_%d",(int)(energys[i]*10000+0.5));
    sprintf(hname,"hp_mcKK_%.4f",energys[i]);
    hdata = (TH1D*)file->Get(hname);
    hdata->SetFillColor(0);
    //hdata->SetTitle("momentum spectrum of data and di-#mu MC");
    hdata->SetTitle("");
    hdata->GetXaxis()->SetTitle("p (GeV/c)");
    hdata->GetYaxis()->SetTitle("Counts");
    // get dimu spectrum
    sprintf(hname,"hp_mcmumu_%.4f",energys[i]);
    hmcmu = (TH1D*)file->Get(hname);

    // get mc KK spectrum
    sprintf(hname,"hp_mcKK2_%.4f",energys[i]);
    hmcKK = (TH1D*)file->Get(hname);
    int binid;
    double maxdata = getMax(hdata,1,50,binid);
    double maxmcKK = getMax(hmcKK,1,50,binid);

    cout<< hdata <<"\t"<<hmcmu<<" ndata,nmcKK "<< maxdata<<" "<<maxmcKK<<endl;
    sprintf(hname,"Canvas_%02d",i+1);
    
    TCanvas *chp = new TCanvas(hname,"Background analysis");
    chp->SetMargin(0.15,0.1,0.15,0.1);
    hdata->SetLineColor(1);
    hdata->SetTitle("");
    hdata->GetXaxis()->SetNdivisions(505);
    hdata->GetXaxis()->SetTitleSize(0.05);
    hdata->GetXaxis()->SetLabelSize(0.05);
    hdata->GetYaxis()->SetNdivisions(505);
    hdata->GetYaxis()->SetTitleSize(0.05);
    hdata->GetYaxis()->SetLabelSize(0.05);
    hdata->GetYaxis()->SetTitleOffset(1.1);

    hdata->Draw("E");
    // scale dimu to data according luminosity
    hmcmu->Scale(luminosity[i]/(nevt[i]/mccross[i]));
    hmcmu->SetLineColor(2);
    hmcmu->SetMarkerColor(2);
    hmcmu->SetFillColor(0);
    hmcmu->Draw("same");
    
    // scale KK MC
    hmcKK->Scale((maxdata-hmcmu->GetBinContent(binid))/maxmcKK);
    hmcKK->SetLineColor(3);
    hmcKK->SetMarkerColor(3);
    hmcKK->SetFillColor(0);
    hmcKK->Draw("same");

    sprintf(hname,"%.4f GeV",energys[i]);
    TLegend* leg = new TLegend(0.4,0.55,0.6,0.85,hname);
    leg->AddEntry(hdata,"data");
    leg->AddEntry(hmcmu,"#mu^{+} #mu^{-} MC");
    leg->AddEntry(hmcKK,"K^{+} K^{-} MC");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();
    sprintf(hname,"Canvas_%02d.pdf",i+1);
    chp->Print(hname);
    delete chp;
    
  }
  return 0;
}

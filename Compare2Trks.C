
double getEne(int idx);

int Compare2Trks()
{
  TFile *file = new TFile("output/output3.root","update");
  TDirectory *dir = (TDirectory*)file->Get("comparasion");
  if (dir==0) dir = file->mkdir("comparasion");
  dir->cd();
  TCanvas *c1 = new TCanvas();
  gStyle->SetOptStat(0);
  for (int enei=1;enei<=21;enei++){
    double ene = getEne(enei);
    char dirname[100];
    sprintf(dirname,"KK_%03d_665",enei);
    std::cout<<"dirname: "<<dirname<<std::endl;
    TDirectory *dirkk = (TDirectory*) file->Get(dirname);
    sprintf(dirname,"pipi_%03d_%5d",enei,(int)(ene*10000));
    std::cout<<"dirname: "<<dirname<<std::endl;
    TDirectory *dirpipi = (TDirectory*) file->Get(dirname);
    sprintf(dirname,"mumu_%03d_%5d",enei,(int)(ene*10000));
    std::cout<<"dirname: "<<dirname<<std::endl;
    TDirectory *dirmumu = (TDirectory*) file->Get(dirname);
    if (dirkk==0 || dirpipi==0 || dirmumu==0) {
      std::cout<<"At ene: "<<ene <<"\tdirkk: "<<dirkk<<"\tdirpipi: "<<dirpipi<<"\tdirmumu: "<<dirmumu<<std::endl;
      return -1;
    }

    TH1D *h1pk  = (TH1D*)dirkk->Get("hpNocut");
    TH1D *h1ppi = (TH1D*)dirpipi->Get("hpNocut");
    TH1D *h1pmu = (TH1D*)dirmumu->Get("hpNocut");
    std::cout<<"1p Nocut: "<<h1pk<<' '<<h1ppi<<' '<<h1pmu<<std::endl;
    h1pk->SetLineColor(1);
    h1pk->SetMarkerColor(1);
    h1pk->SetFillColor(0);
    h1ppi->SetLineColor(2);
    h1ppi->SetMarkerColor(2);
    h1ppi->SetFillColor(0);
    h1pmu->SetLineColor(3);
    h1pmu->SetMarkerColor(3);
    h1pmu->SetFillColor(0);
    h1pk->Draw();
    h1ppi->Draw("same");
    h1pmu->Draw("same");
    TLegend * legend = new TLegend(0.7,0.7,0.9,0.9);
    legend->AddEntry(h1pk,"p_{K}");
    legend->AddEntry(h1ppi,"p_{#pi}");
    legend->AddEntry(h1pmu,"p_{#mu}");
    legend->Draw();
    sprintf(dirname,"pDisNocut_%f",ene);
    dir->WriteTObject(c1,dirname);
    delete legend;
    
    TH2D *h2pk  = (TH2D*)dirkk->Get("h2pNocut");
    TH2D *h2ppi = (TH2D*)dirpipi->Get("h2pNocut");
    TH2D *h2pmu = (TH2D*)dirmumu->Get("h2pNocut");
    std::cout<<"2p Nocut: "<<h2pk<<' '<<h2ppi<<' '<<h2pmu<<std::endl;
    h2pk->SetLineColor(1);
    h2pk->SetMarkerColor(1);
    //h2pk->SetMarkerStyle(2);
    h2pk->SetFillColor(0);
    h2ppi->SetLineColor(2);
    h2ppi->SetMarkerColor(2);
    //h2ppi->SetMarkerStyle(4);
    h2ppi->SetFillColor(0);
    h2pmu->SetLineColor(3);
    h2pmu->SetMarkerColor(3);
    //h2pmu->SetMarkerStyle(5);
    h2pmu->SetFillColor(0);
    h2pk->Draw();
    h2ppi->Draw("same");
    h2pmu->Draw("same");
    legend = new TLegend(0.7,0.7,0.9,0.9);
    legend->AddEntry(h1pk,"p_{K}");
    legend->AddEntry(h1ppi,"p_{#pi}");
    legend->AddEntry(h1pmu,"p_{#mu}");
    legend->Draw();
    sprintf(dirname,"pDisNocut2D_%f",ene);
    dir->WriteTObject(c1,dirname);
    delete legend;

    h1pk  = (TH1D*)dirkk->Get("hpcutTOF");
    h1ppi = (TH1D*)dirpipi->Get("hpcutTOF");
    h1pmu = (TH1D*)dirmumu->Get("hpcutTOF");
    std::cout<<"1p cut: "<<h1pk<<' '<<h1ppi<<' '<<h1pmu<<std::endl;
    h1pk->SetLineColor(1);
    h1pk->SetMarkerColor(1);
    h1pk->SetFillColor(0);
    h1ppi->SetLineColor(2);
    h1ppi->SetMarkerColor(2);
    h1ppi->SetFillColor(0);
    h1pmu->SetLineColor(3);
    h1pmu->SetMarkerColor(3);
    h1pmu->SetFillColor(0);
    h1pk->Draw();
    h1ppi->Draw("same");
    h1pmu->Draw("same");
    legend = new TLegend(0.7,0.7,0.9,0.9);
    legend->AddEntry(h1pk,"p_{K}");
    legend->AddEntry(h1ppi,"p_{#pi}");
    legend->AddEntry(h1pmu,"p_{#mu}");
    legend->Draw();
    sprintf(dirname,"pDis_%f",ene);
    dir->WriteTObject(c1,dirname);
    delete legend;
    
    h2pk  = (TH2D*)dirkk->Get("h2pcutTOF");
    h2ppi = (TH2D*)dirpipi->Get("h2pcutTOF");
    h2pmu = (TH2D*)dirmumu->Get("h2pcutTOF");
    std::cout<<"2p cut: "<<h2pk<<' '<<h2ppi<<' '<<h2pmu<<std::endl;
    h2pk->SetLineColor(1);
    h2pk->SetMarkerColor(1);
    h2pk->SetFillColor(0);
    h2ppi->SetLineColor(2);
    h2ppi->SetMarkerColor(2);
    h2ppi->SetFillColor(0);
    h2pmu->SetLineColor(3);
    h2pmu->SetMarkerColor(3);
    h2pmu->SetFillColor(0);
    h2pk->Draw();
    h2ppi->Draw("same");
    h2pmu->Draw("same");
    legend = new TLegend(0.7,0.7,0.9,0.9);
    legend->AddEntry(h1pk,"p_{K}");
    legend->AddEntry(h1ppi,"p_{#pi}");
    legend->AddEntry(h1pmu,"p_{#mu}");
    legend->Draw();
    sprintf(dirname,"pDis2D_%f",ene);
    dir->WriteTObject(c1,dirname);
    delete legend;
  }

  file->Write();
  delete c1;
  delete file;
  exit(0);
  return 0;
}

double getEne(int idx)
{
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
  return -2;
}

// this code is used to show 3 σ separation of 2 prong events.  bhabha, dimu, pi pi, K K, p p

int CompareP()
{
  char wkdir[1000]={"/Volumes/IR1_CPRA_X6/output"};
  TFile *outfile = new TFile("/Volumes/IR1_CPRA_X6/output/pcom.root","recreate");
  char file[21][100];
  double ene[21];
  sprintf(file[0], "001_2.0000"); ene[0] = 2.0000;
  sprintf(file[1], "002_2.0500"); ene[1] = 2.0500;
  sprintf(file[2], "003_2.1000"); ene[2] = 2.1000;
  sprintf(file[3], "004_2.1500"); ene[3] = 2.1500;
  sprintf(file[4], "005_2.1750"); ene[4] = 2.1750;
  sprintf(file[5], "006_2.2000"); ene[5] = 2.2000;
  sprintf(file[6], "007_2.2324"); ene[6] = 2.2324;
  sprintf(file[7], "008_2.3094"); ene[7] = 2.3094;
  sprintf(file[8], "009_2.3864"); ene[8] = 2.3864;
  sprintf(file[9], "010_2.3960"); ene[9] = 2.3960;
  sprintf(file[10],"011_2.5000"); ene[10]= 2.5000;
  sprintf(file[11],"012_2.6444"); ene[11]= 2.6444;
  sprintf(file[12],"013_2.6464"); ene[12]= 2.6464;
  sprintf(file[13],"014_2.7000"); ene[13]= 2.7000;
  sprintf(file[14],"015_2.8000"); ene[14]= 2.8000;
  sprintf(file[15],"016_2.9000"); ene[15]= 2.9000;
  sprintf(file[16],"017_2.9500"); ene[16]= 2.9500;
  sprintf(file[17],"018_2.9810"); ene[17]= 2.9810;
  sprintf(file[18],"019_3.0000"); ene[18]= 3.0000;
  sprintf(file[19],"020_3.0200"); ene[19]= 3.0200;
  sprintf(file[20],"021_3.0800"); ene[20]= 3.0800;
  double mka = 0.493677;
  double mmu = 0.1057;

  double enee[21]={0};
  double pka[21],pkae[21],pka3e[21],pka5e[21];
  double pmu[21],pmue[21],pmu3e[21],pmu5e[21];

  gStyle->SetOptStat(0);
  for (int enei=0; enei<21; enei++){
    char filename[100];
    sprintf(filename,"%s/mcBB_%s.root",wkdir,file[enei]);
    TFile *filei1 = new TFile(filename);
    TH1D *hpee   = (TH1D*)filei1->Get("hpNocut");
    TH2D *h2pee  = (TH2D*)filei1->Get("h2pNocut");
    TH1D *hpeec  = (TH1D*)filei1->Get("hpcutTOF");
    TH2D *h2peec = (TH2D*)filei1->Get("h2pcutTOF");

    sprintf(filename,"%s/mcdimu_%s.root",wkdir,file[enei]);
    TFile *filei2 = new TFile(filename);
    TH1D *hpmu   = (TH1D*)filei2->Get("hpNocut");
    TH2D *h2pmu  = (TH2D*)filei2->Get("h2pNocut");
    TH1D *hpmuc  = (TH1D*)filei2->Get("hpcutTOF");
    TH2D *h2pmuc = (TH2D*)filei2->Get("h2pcutTOF");

    sprintf(filename,"%s/mchad_%s.root",wkdir,file[enei]);
    TFile *filei3 = new TFile(filename);
    TH1D *hpha   = (TH1D*)filei3->Get("hpNocut");
    TH2D *h2pha  = (TH2D*)filei3->Get("h2pNocut");
    TH1D *hphac  = (TH1D*)filei3->Get("hpcutTOF");
    TH2D *h2phac = (TH2D*)filei3->Get("h2pcutTOF");

    sprintf(filename,"%s/KKsel_%03d_myxs.root",wkdir,enei+1);
    TFile *filei4 = new TFile(filename);
    TH1D *hpka   = (TH1D*)filei4->Get("hpNocut");
    TH2D *h2pka  = (TH2D*)filei4->Get("h2pNocut");
    TH1D *hpkac  = (TH1D*)filei4->Get("hpcutTOF");
    TH2D *h2pkac = (TH2D*)filei4->Get("h2pcutTOF");

    TCanvas *c1 = new TCanvas();
    TLegend *legend = new TLegend(0.15, 0.6,0.4,0.85);
    if (hpka!=0) {hpka->SetFillColor(0);hpka->SetLineColor(1);hpka->SetMarkerColor(1);hpka->Draw(); legend->AddEntry(hpka,"p_{K}");} else continue;
    if (hpee!=0) {hpee->SetFillColor(0);hpee->SetLineColor(2);hpee->SetMarkerColor(2);hpee->Draw("same");legend->AddEntry(hpee,"p_{e}");}
    if (hpmu!=0) {hpmu->SetFillColor(0);hpmu->SetLineColor(3);hpmu->SetMarkerColor(3);hpmu->Draw("same");legend->AddEntry(hpmu,"p_{#mu}");}
    if (hpha!=0) {hpha->SetFillColor(0);hpha->SetLineColor(4);hpha->SetMarkerColor(4);hpha->Draw("same");legend->AddEntry(hpha,"p_{had}");}
    legend->Draw();
    
    TCanvas *c3 = new TCanvas();
    TLegend *lgd3 = new TLegend(0.15, 0.6,0.4,0.85);
    if (hpkac!=0) {hpkac->SetFillColor(0);hpkac->SetLineColor(1);hpkac->SetMarkerColor(1);hpkac->Draw(); lgd3->AddEntry(hpkac,"p_{K}");} else continue;
    if (hpeec!=0) {hpeec->SetFillColor(0);hpeec->SetLineColor(2);hpeec->SetMarkerColor(2);hpeec->Draw("same");lgd3->AddEntry(hpeec,"p_{e}");}
    if (hpmuc!=0) {hpmuc->SetFillColor(0);hpmuc->SetLineColor(3);hpmuc->SetMarkerColor(3);hpmuc->Draw("same");lgd3->AddEntry(hpmuc,"p_{#mu}");}
    if (hphac!=0) {hphac->SetFillColor(0);hphac->SetLineColor(4);hphac->SetMarkerColor(4);hphac->Draw("same");lgd3->AddEntry(hphac,"p_{had}");}
    lgd3->Draw();
    
    TCanvas *c2 = new TCanvas();
    TLegend *lgd2 = new TLegend(0.15, 0.6,0.4,0.85);
    if(h2pka!=0){h2pka->SetFillColor(0);h2pka->SetLineColor(1);h2pka->SetMarkerColor(1);h2pka->Draw();lgd2->AddEntry(h2pka,"p_{K}");} else continue;
    if(h2pee!=0){h2pee->SetFillColor(0);h2pee->SetLineColor(2);h2pee->SetMarkerColor(2);h2pee->Draw("same");lgd2->AddEntry(h2pee,"p_{e}");}
    if(h2pmu!=0){h2pmu->SetFillColor(0);h2pmu->SetLineColor(3);h2pmu->SetMarkerColor(3);h2pmu->Draw("same");lgd2->AddEntry(h2pmu,"p_{#mu}");}
    if(h2pha!=0){h2pha->SetFillColor(0);h2pha->SetLineColor(4);h2pha->SetMarkerColor(4);h2pha->Draw("same");lgd2->AddEntry(h2pha,"p_{had}");}
    lgd2->Draw();
    
    TCanvas *c4 = new TCanvas();
    TLegend *lgd4 = new TLegend(0.15, 0.6,0.4,0.85);
    if (h2pkac!=0) {h2pkac->SetFillColor(0);h2pkac->SetLineColor(1);h2pkac->SetMarkerColor(1);h2pkac->Draw(); lgd4->AddEntry(h2pkac,"p_{K}");} else continue;
    if (h2peec!=0) {h2peec->SetFillColor(0);h2peec->SetLineColor(2);h2peec->SetMarkerColor(2);h2peec->Draw("same");lgd4->AddEntry(h2peec,"p_{e}");}
    if (h2pmuc!=0) {h2pmuc->SetFillColor(0);h2pmuc->SetLineColor(3);h2pmuc->SetMarkerColor(3);h2pmuc->Draw("same");lgd4->AddEntry(h2pmuc,"p_{#mu}");}
    if (h2phac!=0) {h2phac->SetFillColor(0);h2phac->SetLineColor(4);h2phac->SetMarkerColor(4);h2phac->Draw("same");lgd4->AddEntry(h2phac,"p_{had}");}
    lgd4->Draw();
    
    sprintf(filename,"%s_1p",file[enei]);
    outfile->WriteTObject(c1,filename);
    sprintf(filename,"%s_2p",file[enei]);
    outfile->WriteTObject(c2,filename);
    sprintf(filename,"%s_1p_cut",file[enei]);
    outfile->WriteTObject(c3,filename);
    sprintf(filename,"%s_2p_cut",file[enei]);
    outfile->WriteTObject(c4,filename);

    double p_ka = sqrt(pow(ene[enei]/2,2)-mka*mka);
    double p_mu = sqrt(pow(ene[enei]/2,2)-mmu*mmu);
    //std::cout<<p_ka<<' '<<p_mu<<std::endl;
    hpka->Fit("gaus","R","",p_ka-0.1,p_ka+0.1);
    hpmu->Fit("gaus","R","",p_mu-0.1,p_mu+0.1);
    pka[enei] = hpka->GetFunction("gaus")->GetParameter(1);
    pkae[enei]= hpka->GetFunction("gaus")->GetParameter(2);
    pka3e[enei]= 3*pkae[enei]/2.6*3;
    pka5e[enei]= 5*pkae[enei]/2.6*3;
    pmu[enei] = hpmu->GetFunction("gaus")->GetParameter(1);
    pmue[enei]= hpmu->GetFunction("gaus")->GetParameter(2);
    pmu3e[enei]= 3*pmue[enei]/2.6*3;
    pmu5e[enei]= 5*pmue[enei]/2.6*3;
    std::cout<<p_ka<<' ' << pka[enei]<<" "<< pkae[enei]<<std::endl;
    std::cout<<p_mu<<' ' << pmu[enei]<<" "<< pmue[enei]<<std::endl;
    
    delete filei1;
    delete filei2;
    delete filei3;
    delete filei4;
    delete c1;
    delete c2;
    delete c3;
    delete c4;
    delete legend;
    delete lgd2;
    delete lgd3;
    delete lgd4;
  }

  new TCanvas();
  TGraphErrors *graph1 = new TGraphErrors(21,ene,pka,enee,pkae);
  graph1->SetFillColor(0);
  graph1->SetMarkerStyle(24);
  TGraphErrors *graph1_1 = new TGraphErrors(21,ene,pka,enee,pkae);
  graph1_1->SetFillColor(3);
  TGraphErrors *graph1_3 = new TGraphErrors(21,ene,pka,enee,pka3e);
  graph1_3->SetFillColor(4);
  TGraphErrors *graph1_5 = new TGraphErrors(21,ene,pka,enee,pka5e);
  graph1_5->SetFillColor(2);
  //graph1->Draw("A3LP");

  TGraphErrors *graph2 = new TGraphErrors(21,ene,pmu,enee,pmue);
  graph2->SetFillColor(0);
  graph2->SetMarkerStyle(25);
  TGraphErrors *graph2_1 = new TGraphErrors(21,ene,pmu,enee,pmue);
  graph2_1->SetFillColor(3);
  TGraphErrors *graph2_3 = new TGraphErrors(21,ene,pmu,enee,pmu3e);
  graph2_3->SetFillColor(4);
  TGraphErrors *graph2_5 = new TGraphErrors(21,ene,pmu,enee,pmu5e);
  graph2_5->SetFillColor(2);
  //graph2->Draw("L3P");
  
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("p of 2 prong process");
  mg->Add(graph1_5);
  mg->Add(graph2_5);
  mg->Add(graph1_3);
  mg->Add(graph2_3);
  mg->Add(graph1_1);
  mg->Add(graph2_1);
  mg->Draw("ALP3");
  mg->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  mg->GetYaxis()->SetTitle("p (GeV/c)");
  graph1->Draw("PL");
  graph2->Draw("PL");

  TLegend *legend = new TLegend(0.15,0.65,0.45,0.85);
  legend->AddEntry(graph1,"p_{K}");
  legend->AddEntry(graph2,"p_{#mu}");
  legend->AddEntry(graph1_1,"1 #sigma");
  legend->AddEntry(graph1_3,"3 #sigma");
  legend->AddEntry(graph1_5,"5 #sigma");
  legend->SetMargin(0.6);
  legend->Draw();

}

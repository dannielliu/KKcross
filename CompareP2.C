// this code is used to show 3 σ separation of 2 prong events.  bhabha, dimu, pi pi, K K, p p

int CompareP2()
{
  //char wkdir[1000]={"/Volumes/IR1_CPRA_X6/output"};
  //TFile *outfile = new TFile("/Volumes/IR1_CPRA_X6/output/pcom.root","recreate");
  //char file[21][100];
  double ene[21];
  ene[0] = 2.0000;
  ene[1] = 2.0500;
  ene[2] = 2.1000;
  ene[3] = 2.1500;
  ene[4] = 2.1750;
  ene[5] = 2.2000;
  ene[6] = 2.2324;
  ene[7] = 2.3094;
  ene[8] = 2.3864;
  ene[9] = 2.3960;
  ene[10]= 2.5000;
  ene[11]= 2.6444;
  ene[12]= 2.6464;
  ene[13]= 2.7000;
  ene[14]= 2.8000;
  ene[15]= 2.9000;
  ene[16]= 2.9500;
  ene[17]= 2.9810;
  ene[18]= 3.0000;
  ene[19]= 3.0200;
  ene[20]= 3.0800;
  double mka = 0.493677;
  double mmu = 0.1057;

  double enee[21]={0};
  double pka[21],pkae[21],pka3e[21],pka5e[21];
  double pmu[21],pmue[21],pmu3e[21],pmu5e[21];

  ifstream infile("fitpar");
  gStyle->SetOptStat(0);
  for (int enei=0; enei<21; enei++){
    char filename[100];
    char line[1000];
    infile.getline(line,1000);
    istringstream iss;
    iss.str(line);
    iss >> ene[enei] >> pka[enei] >> pkae[enei] >> pmu[enei] >> pmue[enei];

    double p_ka = sqrt(pow(ene[enei]/2,2)-mka*mka);
    double p_mu = sqrt(pow(ene[enei]/2,2)-mmu*mmu);
    //std::cout<<p_ka<<' '<<p_mu<<std::endl;
    pka3e[enei]= 3*pkae[enei];
    pka5e[enei]= 5*pkae[enei];
    pmu3e[enei]= 3*pmue[enei];
    pmu5e[enei]= 5*pmue[enei];
    //std::cout<<p_ka<<'\t' << pka[enei]<<"\t"<< pkae[enei] <<
    std::cout<<ene[enei]<<"\t pka+3sigma "<<pka[enei]+3*pkae[enei]<<std::endl;
    //std::cout<<p_mu<<'\t' << pmu[enei]<<"\t"<< pmue[enei]<<std::endl;
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

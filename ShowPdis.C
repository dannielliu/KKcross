int ShowPdis()
{
  char line[100];
  double ene[22];
  double pe[22];
  double pmu[22];
  double ppi[22];
  double pk[22];
  double pp[22];
  double deltakpi[22];
  double deltakp[22];
  ifstream infile("pdis");
  infile.getline(line,1000);
  for (int i=0;i<21;i++){
    infile >> ene[i] >> pe[i] >> pmu[i] >> ppi[i] >> pk[i] >> pp[i];
    deltakpi[i] = ppi[i]-pk[i];
    deltakp[i] = pk[i]-pp[i];
  }
  int n=21;
  TGraph *graph1 = new TGraph(n,ene,pe);
  TGraph *graph2 = new TGraph(n,ene,pmu);
  TGraph *graph3 = new TGraph(n,ene,ppi);
  TGraph *graph4 = new TGraph(n,ene,pk);
  TGraph *graph5 = new TGraph(n,ene,pp);
  
  graph1->SetFillColor(0);
  graph2->SetFillColor(0);
  graph3->SetFillColor(0);
  graph4->SetFillColor(0);
  graph5->SetFillColor(0);
  graph1->SetLineColor(1);
  graph1->SetMarkerColor(1);
  graph1->SetMarkerStyle(25);
  graph2->SetLineColor(2);
  graph2->SetMarkerColor(2);
  graph2->SetMarkerStyle(2);
  graph3->SetLineColor(3);
  graph3->SetMarkerColor(3);
  graph3->SetMarkerStyle(27);
  graph4->SetLineColor(4);
  graph4->SetMarkerColor(4);
  graph4->SetMarkerStyle(4);
  graph5->SetLineColor(6);
  graph5->SetMarkerColor(6);
  graph5->SetMarkerStyle(5);

  graph1->GetYaxis()->SetRangeUser(0,2);
  graph1->GetYaxis()->SetTitle("p (GeV/c)");
  graph1->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  graph1->Draw("AP");
  graph2->Draw("P");
  graph3->Draw("P");
  graph4->Draw("P");
  graph5->Draw("P");
  
  TLegend *legend = new TLegend(0.13,0.65,0.23,0.85);
  legend->AddEntry(graph1,"p_{e}");
  legend->AddEntry(graph2,"p_{#mu}");
  legend->AddEntry(graph3,"p_{#pi}");
  legend->AddEntry(graph4,"p_{K}");
  legend->AddEntry(graph5,"p_{p}");
  legend->Draw();

  
  new TCanvas();
  TGraph *graph6 = new TGraph(n,ene,deltakpi);
  TGraph *graph7 = new TGraph(n,ene,deltakp);
  
  graph6->SetFillColor(0);
  graph7->SetFillColor(0);
  graph6->SetLineColor(1);
  graph6->SetMarkerColor(1);
  graph6->SetMarkerStyle(25);
  graph7->SetLineColor(2);
  graph7->SetMarkerColor(2);
  graph7->SetMarkerStyle(2);
  
  graph6->GetYaxis()->SetRangeUser(0,1);
  graph6->GetYaxis()->SetTitle("#delta p (GeV/c)");
  graph6->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  graph6->Draw("AP");
  graph7->Draw("P");
  
  TLegend *legend2 = new TLegend(0.13,0.65,0.23,0.85);
  legend2->AddEntry(graph6,"|p_{#pi} - p_{K}|");
  legend2->AddEntry(graph7,"|p_{K} - p_{p}|");
  legend2->Draw();




}

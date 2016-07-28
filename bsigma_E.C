
int bsigma_E()
{
  //ifstream ist("es.txt");
  ifstream ist1("fitpar_sigwid_sigma.txt");
  ifstream ist2("fitpar_bcksmr_sigma.txt");
  ifstream ist3("fitpar_bckwid_sigma.txt");
  ifstream ist4("fitpar_sigsmr_sigma.txt");
  
  TFile* file = new TFile("sigma.root","recreate");
  double ene[22];
  double eer[22];
  double sig[22];
  double err[22];

  istringstream iss;
  char line[1000];
  int id=0;
  while (!ist1.eof() && id<22){
    double energy, sigma, error;
    ist1.getline(line,1000);
    if (line[0]=='\0') continue;
    iss.clear();
    iss.str(line);
    iss>>energy>>sigma>>error;
    ene[id] = energy;
    eer[id] = 0;
    sig[id] = sigma;
    err[id] = error;
    id++;
  }
  TGraphErrors* graph1 = new TGraphErrors(22,ene,sig, eer, err);
  graph1->SetMarkerStyle(25);
  //graph1->Draw("AP");
  file->WriteTObject(graph1,"signal_width_gaus");
  
  id = 0;
  while (!ist2.eof() && id<22){
    double energy, sigma, error;
    ist2.getline(line,1000);
    if (line[0]=='\0') continue;
    iss.clear();
    iss.str(line);
    iss>>energy>>sigma>>error;
    ene[id] = energy;
    eer[id] = 0;
    sig[id] = sigma;
    err[id] = error;
    id++;
  }
  TGraphErrors* graph2 = new TGraphErrors(22,ene,sig, eer, err);
  graph2->SetMarkerStyle(25);
  //graph1->Draw("AP");
  file->WriteTObject(graph2,"back_sigma_smear");

  id = 0;
  while (!ist3.eof() && id<22){
    double energy, sigma, error;
    ist3.getline(line,1000);
    if (line[0]=='\0') continue;
    iss.clear();
    iss.str(line);
    iss>>energy>>sigma>>error;
    ene[id] = energy;
    eer[id] = 0;
    sig[id] = sigma;
    err[id] = error;
    id++;
  }
  TGraphErrors* graph3 = new TGraphErrors(22,ene,sig, eer, err);
  graph1->SetMarkerStyle(25);
  //graph1->Draw("AP");
  file->WriteTObject(graph3,"back_width_gaus");

  id = 0;
  while (!ist4.eof() && id<22){
    double energy, sigma, error;
    ist4.getline(line,1000);
    if (line[0]=='\0') continue;
    iss.clear();
    iss.str(line);
    iss>>energy>>sigma>>error;
    ene[id] = energy;
    eer[id] = 0;
    sig[id] = sigma;
    err[id] = error;
    id++;
  }
  TGraphErrors* graph4 = new TGraphErrors(22,ene,sig, eer, err);
  graph4->SetMarkerStyle(25);
  //graph1->Draw("AP");
  file->WriteTObject(graph4,"signal_sigma_smear");

  id = 0;
  ist1.clear();
  ist1.seekg(0,ios::beg);
  while (!ist1.eof() && id<22){
    double energy, sigma, error, cbsigma, cbsigmaerr;
    ist1.getline(line,1000);
    if (line[0]=='\0') continue;
    iss.clear();
    iss.str(line);
    iss>>energy>>sigma>>error>>cbsigma>>cbsigmaerr;
    ene[id] = energy;
    eer[id] = 0;
    sig[id] = cbsigma;
    err[id] = cbsigmaerr;
    id++;
  }
  TGraphErrors* graph5 = new TGraphErrors(22,ene,sig, eer, err);
  graph5->SetMarkerStyle(25);
  //graph1->Draw("AP");
  file->WriteTObject(graph5,"signal_width_cbsh");

  id = 0;
  ist3.clear();
  ist3.seekg(0,ios::beg);
  while (!ist3.eof() && id<22){
    double energy, sigma, error, cbsigma, cbsigmaerr;
    ist3.getline(line,1000);
    if (line[0]=='\0') continue;
    iss.clear();
    iss.str(line);
    iss>>energy>>sigma>>error>>cbsigma>>cbsigmaerr;
    ene[id] = energy;
    eer[id] = 0;
    sig[id] = cbsigma;
    err[id] = cbsigmaerr;
    id++;
  }
  TGraphErrors* graph6 = new TGraphErrors(22,ene,sig, eer, err);
  graph6->SetMarkerStyle(25);
  //graph1->Draw("AP");
  file->WriteTObject(graph6,"back_width_cbsh");

}

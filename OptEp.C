int OptEp()
{
  char wkdir[1000] = {"/Volumes/IR1_CPRA_X6/output"};
  TFile *outfile = new TFile("/Volumes/IR1_CPRA_X6/output_cut/epdis.root","recreate");
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

  TFile* filein = new TFile("/Volumes/IR1_CPRA_X6/output_cut/output_mc.root");
//ofstream cutpar("cutpar");
//cutpar<<"ene \tp \tE/p"<<std::endl;

  double ratio1[1000];
  double ratio2[1000];
  double eta[1000];

  gStyle->SetOptStat(0);
  for (int enei=0; enei<21; enei++){
  //cutpar<<ene[enei];
    
    char filename[100];
    sprintf(filename,"mcBB_%s",file[enei]);
    TDirectory *filei1 = (TDirectory*)filein->Get(filename);
    TH1D *hepee    = (TH1D*)filei1->Get("hep");
    
    sprintf(filename,"mcdimu_%s",file[enei]);
    TDirectory *filei2 = (TDirectory*)filein->Get(filename);
    TH1D *hepmu    = (TH1D*)filei2->Get("hep");

    sprintf(filename,"mchad_%s",file[enei]);
    TDirectory *filei3 = (TDirectory*)filein->Get(filename);
    TH1D *hepha    = (TH1D*)filei3->Get("hep");

    sprintf(filename,"KKsel_%03d_myxs",enei+1);
    TDirectory *filei4 = (TDirectory*)filein->Get(filename);
    TH1D *hepka    = (TH1D*)filei4->Get("hep");
    
    sprintf(filename,"mcKK_%03d_665",enei+1);
    TDirectory *filei5 = (TDirectory*)filein->Get(filename);
    TH1D *hepk2    = (TH1D*)filei5->Get("hep");

    TCanvas *c3 = new TCanvas();
    TLegend *lgd3 = new TLegend(0.25, 0.6,0.45,0.85);
    double hepkamax = hepka->GetMaximum();
    double hepeemax = hepee->GetMaximum();
    double hepmumax = hepmu->GetMaximum();
    double hephamax = hepha->GetMaximum();
    double hepk2max = hepk2->GetMaximum();
    hepee->Scale(hepkamax/hepeemax);
    hepmu->Scale(hepkamax/hepmumax);
    hepk2->Scale(hepkamax/hepk2max);
    if(hepka!=0){hepka->SetFillColor(0);hepka->SetLineColor(1);hepka->SetMarkerColor(1);hepka->Draw();lgd3->AddEntry(hepka,"(E/p)_{K}");} else continue;
    if(hepee!=0){hepee->SetFillColor(0);hepee->SetLineColor(2);hepee->SetMarkerColor(2);hepee->Draw("same");lgd3->AddEntry(hepee,"(E/p)_{e}");}
    if(hepmu!=0){hepmu->SetFillColor(0);hepmu->SetLineColor(3);hepmu->SetMarkerColor(3);hepmu->Draw("same");lgd3->AddEntry(hepmu,"(E/p)_{#mu}");}
    if(hepha!=0){hepha->SetFillColor(0);hepha->SetLineColor(4);hepha->SetMarkerColor(4);hepha->Draw("same");lgd3->AddEntry(hepha,"(E/p)_{had}");}
    if(hepk2!=0){hepk2->SetFillColor(0);hepk2->SetLineColor(7);hepk2->SetMarkerColor(7);hepk2->Draw("same");lgd3->AddEntry(hepk2,"(E/p)_{K} (no ISR)");}
    lgd3->Draw();
    double Noka3 = hepk2->Integral(1,1000);
    double Nokat = hepee->Integral(1,1000);
    if (Noka3<10) return;
    double Ninety;
    for (int i=1;i<=1000;i++){
      ratio1[i-1] = hepk2->Integral(1,i)/Noka3;
      ratio2[i-1] = hepee->Integral(1,i)/Nokat;
      //Nokat *=10;
      eta[i-1]    = sqrt(Noka3*ratio1[i-1]+Nokat*100*ratio2[i-1])/(Noka3*ratio1[i-1]);
    }
    TGraph *graph = new TGraph(1000,ratio1,eta);
    graph->SetMarkerStyle(5);
    graph->Draw("AP");
    graph->GetXaxis()->SetTitle("#epsilon_{sig}");
    graph->GetYaxis()->SetTitle("#eta");
    std::cout<<"signal No: "<< Noka3 << "background No:" << Nokat << std::endl;

    sprintf(filename,"%s_ep",file[enei]);
    outfile->WriteTObject(c3,filename);
    sprintf(filename,"output_cut/%s_ep.pdf",file[enei]);
    c3->SetLogy();
    c3->Print(filename);
    sprintf(filename,"%s_2ep",file[enei]);
   
    delete filei1;
    delete filei2;
    delete filei3;
    delete filei4;
    delete c3;
    delete lgd3;
  //cutpar<<std::endl;
  }

  outfile->Close(); 
}

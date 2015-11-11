int SetCuts()
{
  char wkdir[1000] = {"/Volumes/IR1_CPRA_X6/output"};
  TFile *outfile = new TFile("/Volumes/IR1_CPRA_X6/output_cut/vardis.root","recreate");
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
  TFile* filein = new TFile("/Volumes/IR1_CPRA_X6/output_cut/output_mc.root");
  //ofstream cutpar("cutpar");
  //cutpar<<"ene \tp \tE/p"<<std::endl;

  gStyle->SetOptStat(0);
  for (int enei=0; enei<1; enei++){
    //cutpar<<ene[enei];
    
    char filename[100];
    sprintf(filename,"mcBB_%s",file[enei]);
    TDirectory *filei1 = (TDirectory*)filein->Get(filename);
    TH1D *hpee     = (TH1D*)filei1->Get("h1p");
    TH2D *h2pee    = (TH2D*)filei1->Get("h2p");
    TH1D *hepee    = (TH1D*)filei1->Get("hep");
    TH2D *h2epee   = (TH2D*)filei1->Get("h2ep");
    TH1D *hcos1ee  = (TH1D*)filei1->Get("hcostheta1");
    TH2D *hcos2ee  = (TH2D*)filei1->Get("hcostheta2");
    TH1D *hthetaee = (TH1D*)filei1->Get("htheta");
    TH1D *hdtofee  = (TH1D*)filei1->Get("hdtof");

    sprintf(filename,"mcdimu_%s",file[enei]);
    TDirectory *filei2 = (TDirectory*)filein->Get(filename);
    TH1D *hpmu     = (TH1D*)filei2->Get("h1p");
    TH2D *h2pmu    = (TH2D*)filei2->Get("h2p");
    TH1D *hepmu    = (TH1D*)filei2->Get("hep");
    TH2D *h2epmu   = (TH2D*)filei2->Get("h2ep");
    TH1D *hcos1mu  = (TH1D*)filei2->Get("hcostheta1");
    TH2D *hcos2mu  = (TH2D*)filei2->Get("hcostheta2");
    TH1D *hthetamu = (TH1D*)filei2->Get("htheta");
    TH1D *hdtofmu  = (TH1D*)filei2->Get("hdtof");

    sprintf(filename,"mchad_%s",file[enei]);
    TDirectory *filei3 = (TDirectory*)filein->Get(filename);
    TH1D *hpha     = (TH1D*)filei3->Get("h1p");
    TH2D *h2pha    = (TH2D*)filei3->Get("h2p");
    TH1D *hepha    = (TH1D*)filei3->Get("hep");
    TH2D *h2epha   = (TH2D*)filei3->Get("h2ep");
    TH1D *hcos1ha  = (TH1D*)filei3->Get("hcostheta1");
    TH2D *hcos2ha  = (TH2D*)filei3->Get("hcostheta2");
    TH1D *hthetaha = (TH1D*)filei3->Get("htheta");
    TH1D *hdtofha  = (TH1D*)filei3->Get("hdtof");

    sprintf(filename,"KKsel_%03d_myxs",enei+1);
    TDirectory *filei4 = (TDirectory*)filein->Get(filename);
    TH1D *hpka     = (TH1D*)filei4->Get("h1p");
    TH2D *h2pka    = (TH2D*)filei4->Get("h2p");
    TH1D *hepka    = (TH1D*)filei4->Get("hep");
    TH2D *h2epka   = (TH2D*)filei4->Get("h2ep");
    TH1D *hcos1ka  = (TH1D*)filei4->Get("hcostheta1");
    TH2D *hcos2ka  = (TH2D*)filei4->Get("hcostheta2");
    TH1D *hthetaka = (TH1D*)filei4->Get("htheta");
    TH1D *hdtofka  = (TH1D*)filei4->Get("hdtof");
    
    sprintf(filename,"KKsel_%03d_radflag0",enei+1);
    TDirectory *filei5 = (TDirectory*)filein->Get(filename);
    TH1D *hpk2     = (TH1D*)filei5->Get("h1p");
    TH2D *h2pk2    = (TH2D*)filei5->Get("h2p");
    TH1D *hepk2    = (TH1D*)filei5->Get("hep");
    TH2D *h2epk2   = (TH2D*)filei5->Get("h2ep");
    TH1D *hcos1k2  = (TH1D*)filei5->Get("hcostheta1");
    TH2D *hcos2k2  = (TH2D*)filei5->Get("hcostheta2");
    TH1D *hthetak2 = (TH1D*)filei5->Get("htheta");
    TH1D *hdtofk2  = (TH1D*)filei5->Get("hdtof");

    TCanvas *c1 = new TCanvas();
    TLegend *lgd1 = new TLegend(0.15, 0.75,0.4,0.9);
    double hpkamax = hpka->GetMaximum();
    double hpk2max = hpk2->GetMaximum();
    double hpeemax = hpee->GetMaximum();
    double hpmumax = hpmu->GetMaximum();
    double hphamax = hpha->GetMaximum();
    hpee->Scale(hpkamax/hpeemax);
    hpmu->Scale(hpkamax/hpmumax);
    hpk2->Scale(hpkamax/hpk2max);
    if (hpka!=0) {hpka->SetFillColor(0);hpka->SetLineColor(1);hpka->SetMarkerColor(1);hpka->Draw(); lgd1->AddEntry(hpka,"p_{K}");} else continue;
    if (hpee!=0) {hpee->SetFillColor(0);hpee->SetLineColor(2);hpee->SetMarkerColor(2);hpee->Draw("same");lgd1->AddEntry(hpee,"p_{e}");}
    if (hpmu!=0) {hpmu->SetFillColor(0);hpmu->SetLineColor(3);hpmu->SetMarkerColor(3);hpmu->Draw("same");lgd1->AddEntry(hpmu,"p_{#mu}");}
    if (hpha!=0) {hpha->SetFillColor(0);hpha->SetLineColor(4);hpha->SetMarkerColor(4);hpha->Draw("same");lgd1->AddEntry(hpha,"p_{had}");}
    if (hpk2!=0) {hpk2->SetFillColor(0);hpk2->SetLineColor(7);hpk2->SetMarkerColor(7);hpk2->Draw("same");lgd1->AddEntry(hpk2,"p_{K} (no ISR)");}
    lgd1->Draw();
    double Noka1 = hpk2->Integral(1,1000);
    double pexp = sqrt(pow(ene[enei]/2,2)-mka*mka);
    hpka->Fit("gaus","R","",pexp-0.03,pexp+0.03);
    double meanp = hpka->GetFunction("gaus")->GetParameter(1);
    double sigmap = hpka->GetFunction("gaus")->GetParameter(2);
    double limit1 = meanp+3*sigmap;
    TArrow *arr1 = new TArrow(limit1,0.5*hpkamax,limit1,0,0.05);
    arr1->SetLineColor(6);
    arr1->SetLineWidth(2);
    arr1->Draw();
    std::cout<<file[enei]<<" with 3 sigma upper limit at "<<limit1<<", Total Entries is "<< Noka1<<std::endl;
    //cutpar<<"\t"<<limit1;
    double surv1 = hpk2->Integral(1,(int)(limit1/2.0*1000));
    std::cout<<"p cut, Survive ratio is "<<hpka->Integral(1,(int)(limit1/2.0*1000))/hpka->Integral(1,1000)<<", noISR "<<surv1/Noka1<<std::endl;

    TCanvas *c2 = new TCanvas();
    TLegend *lgd2 = new TLegend(0.15, 0.6,0.4,0.85);
    if(h2pka!=0){h2pka->SetFillColor(0);h2pka->SetLineColor(1);h2pka->SetMarkerColor(1);h2pka->Draw();lgd2->AddEntry(h2pka,"p_{K}");} else continue;
    if(h2pee!=0){h2pee->SetFillColor(0);h2pee->SetLineColor(2);h2pee->SetMarkerColor(2);h2pee->Draw("same");lgd2->AddEntry(h2pee,"p_{e}");}
    if(h2pmu!=0){h2pmu->SetFillColor(0);h2pmu->SetLineColor(3);h2pmu->SetMarkerColor(3);h2pmu->Draw("same");lgd2->AddEntry(h2pmu,"p_{#mu}");}
    if(h2pha!=0){h2pha->SetFillColor(0);h2pha->SetLineColor(4);h2pha->SetMarkerColor(4);h2pha->Draw("same");lgd2->AddEntry(h2pha,"p_{had}");}
    if(h2pk2!=0){h2pk2->SetFillColor(0);h2pk2->SetLineColor(7);h2pk2->SetMarkerColor(7);h2pk2->Draw("same");lgd2->AddEntry(h2pk2,"p_{K} (no ISR)");}
    lgd2->Draw();
    
    TCanvas *c3 = new TCanvas();
    TLegend *lgd3 = new TLegend(0.75, 0.6,0.9,0.85);
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
    double Noee3 = hepee->Integral(1,1000);
    if (Noka3<10) return;
    double Ninety;
  //for (int i=1;i<=1000;i++){
  //  if (hepk2->Integral(1,i)<0.9*Noka3) continue;
  //  Ninety = 1.5*i/1000;
  //  TArrow *arr3 = new TArrow(Ninety,0.5*hepkamax,Ninety,0,0.05);
  //  arr3->SetLineColor(6);
  //  arr3->SetLineWidth(2);
  //  arr3->Draw();
  //  std::cout<<file[enei]<<" with 90\% E/p at "<<Ninety<<", Total Entries is "<< Noka3<<std::endl;
  //  break;
  //}
    for (int i=1;i<=1000;i++){
      if (hepee->Integral(1,i)<0.005*Noee3) continue;
      Ninety = 1.5*i/1000;
      TArrow *arr3 = new TArrow(Ninety,0.5*hepkamax,Ninety,0,0.05);
      arr3->SetLineColor(6);
      arr3->SetLineWidth(2);
      arr3->Draw();
      std::cout<<file[enei]<<" with 90\% ee veto, E/p at "<<Ninety<<", Total Entries is "<< Noka3<<std::endl;
      break;
    }
    //cutpar<<"\t"<<Ninety;
    double surv3 = hepk2->Integral(1,(int)(Ninety/1.5*1000));
    std::cout<<"E/p cut, Survive ratio is "<<hepka->Integral(1,(int)(Ninety/1.5*1000))/hepka->Integral(1,1000)<<", noISR " << surv3/Noka3<<std::endl;
    std::cout<<"E/p cut, bhabha survive ratio is "<<hepee->Integral(1,(int)(Ninety/1.5*1000))/hepee->Integral(1,1000)<<std::endl;
    
    TCanvas *c4 = new TCanvas();
    TLegend *lgd4 = new TLegend(0.15, 0.6,0.4,0.85);
    if(h2epka!=0){h2epka->SetFillColor(0);h2epka->SetLineColor(1);h2epka->SetMarkerColor(1);h2epka->Draw();lgd4->AddEntry(h2epka,"(E/p)_{K}");} else continue;
    if(h2epee!=0){h2epee->SetFillColor(0);h2epee->SetLineColor(2);h2epee->SetMarkerColor(2);h2epee->Draw("same");lgd4->AddEntry(h2epee,"(E/p)_{e}");}
    if(h2epmu!=0){h2epmu->SetFillColor(0);h2epmu->SetLineColor(3);h2epmu->SetMarkerColor(3);h2epmu->Draw("same");lgd4->AddEntry(h2epmu,"(E/p)_{#mu}");}
    if(h2epha!=0){h2epha->SetFillColor(0);h2epha->SetLineColor(4);h2epha->SetMarkerColor(4);h2epha->Draw("same");lgd4->AddEntry(h2epha,"(E/p)_{had}");}
    if(h2epk2!=0){h2epk2->SetFillColor(0);h2epk2->SetLineColor(7);h2epk2->SetMarkerColor(7);h2epk2->Draw("same");lgd4->AddEntry(h2epk2,"(E/p)_{K} (no ISR)");}
    lgd4->Draw();
    
    TCanvas *c5 = new TCanvas();
    TLegend *lgd5 = new TLegend(0.15, 0.6,0.4,0.85);
    double hcos1kamax = hcos1ka->GetMaximum();
    double hcos1eemax = hcos1ee->GetMaximum();
    double hcos1mumax = hcos1mu->GetMaximum();
    double hcos1hamax = hcos1ha->GetMaximum();
    double hcos1k2max = hcos1k2->GetMaximum();
    hcos1ee->Scale(hcos1kamax/hcos1eemax);
    hcos1mu->Scale(hcos1kamax/hcos1mumax);
    hcos1ha->Scale(hcos1kamax/hcos1hamax);
    hcos1k2->Scale(hcos1kamax/hcos1k2max);
    if(hcos1ka!=0){hcos1ka->SetFillColor(0);hcos1ka->SetLineColor(1);hcos1ka->SetMarkerColor(1);hcos1ka->Draw();lgd5->AddEntry(hcos1ka,"cos#theta1_{K}");} else continue;
    if(hcos1ee!=0){hcos1ee->SetFillColor(0);hcos1ee->SetLineColor(2);hcos1ee->SetMarkerColor(2);hcos1ee->Draw("same");lgd5->AddEntry(hcos1ee,"cos#theta1_{e}");}
    if(hcos1mu!=0){hcos1mu->SetFillColor(0);hcos1mu->SetLineColor(3);hcos1mu->SetMarkerColor(3);hcos1mu->Draw("same");lgd5->AddEntry(hcos1mu,"cos#theta1_{#mu}");}
    if(hcos1ha!=0){hcos1ha->SetFillColor(0);hcos1ha->SetLineColor(4);hcos1ha->SetMarkerColor(4);hcos1ha->Draw("same");lgd5->AddEntry(hcos1ha,"cos#theta1_{had}");}
    if(hcos1k2!=0){hcos1k2->SetFillColor(0);hcos1k2->SetLineColor(7);hcos1k2->SetMarkerColor(7);hcos1k2->Draw("same");lgd5->AddEntry(hcos1k2,"cos#theta1_{K} (no ISR)");}
    lgd5->Draw();
    
    TCanvas *c6 = new TCanvas();
    TLegend *lgd6 = new TLegend(0.15, 0.6,0.4,0.85);
    double hcos2kamax = hcos2ka->GetMaximum();
    double hcos2eemax = hcos2ee->GetMaximum();
    double hcos2mumax = hcos2mu->GetMaximum();
    double hcos2hamax = hcos2ha->GetMaximum();
    double hcos2k2max = hcos2k2->GetMaximum();
    hcos2ee->Scale(hcos2kamax/hcos2eemax);
    hcos2mu->Scale(hcos2kamax/hcos2mumax);
    hcos2ha->Scale(hcos2kamax/hcos2hamax);
    hcos2k2->Scale(hcos2kamax/hcos2k2max);
    if(hcos2ka!=0){hcos2ka->SetFillColor(0);hcos2ka->SetLineColor(1);hcos2ka->SetMarkerColor(1);hcos2ka->Draw();lgd6->AddEntry(hcos2ka,"cos#theta2_{K}");} else continue;
    if(hcos2ee!=0){hcos2ee->SetFillColor(0);hcos2ee->SetLineColor(2);hcos2ee->SetMarkerColor(2);hcos2ee->Draw("same");lgd6->AddEntry(hcos2ee,"cos#theta2_{e}");}
    if(hcos2mu!=0){hcos2mu->SetFillColor(0);hcos2mu->SetLineColor(3);hcos2mu->SetMarkerColor(3);hcos2mu->Draw("same");lgd6->AddEntry(hcos2mu,"cos#theta2_{#mu}");}
    if(hcos2ha!=0){hcos2ha->SetFillColor(0);hcos2ha->SetLineColor(4);hcos2ha->SetMarkerColor(4);hcos2ha->Draw("same");lgd6->AddEntry(hcos2ha,"cos#theta2_{had}");}
    if(hcos2k2!=0){hcos2k2->SetFillColor(0);hcos2k2->SetLineColor(7);hcos2k2->SetMarkerColor(7);hcos2k2->Draw("same");lgd6->AddEntry(hcos2k2,"cos#theta2_{K} (no ISR)");}
    lgd6->Draw();
    
    TCanvas *c7 = new TCanvas();
    TLegend *lgd7 = new TLegend(0.15, 0.6,0.4,0.85);
    double hthetakamax = hthetaka->GetMaximum();
    double hthetaeemax = hthetaee->GetMaximum();
    double hthetamumax = hthetamu->GetMaximum();
    double hthetahamax = hthetaha->GetMaximum();
    double hthetak2max = hthetak2->GetMaximum();
    hthetaee->Scale(hthetakamax/hthetaeemax);
    hthetamu->Scale(hthetakamax/hthetamumax);
    hthetaha->Scale(hthetakamax/hthetahamax);
    hthetak2->Scale(hthetakamax/hthetak2max);
    hthetaka->SetBins(100,170,180);
    hthetaee->SetBins(100,170,180);
    hthetamu->SetBins(100,170,180);
    hthetaha->SetBins(100,170,180);
    hthetak2->SetBins(100,170,180);
    if(hthetaka!=0){hthetaka->SetFillColor(0);hthetaka->SetLineColor(1);hthetaka->SetMarkerColor(1);hthetaka->Draw();lgd7->AddEntry(hthetaka,"#theta_{K}");} else continue;
    if(hthetaee!=0){hthetaee->SetFillColor(0);hthetaee->SetLineColor(2);hthetaee->SetMarkerColor(2);hthetaee->Draw("same");lgd7->AddEntry(hthetaee,"#theta_{e}");}
    if(hthetamu!=0){hthetamu->SetFillColor(0);hthetamu->SetLineColor(3);hthetamu->SetMarkerColor(3);hthetamu->Draw("same");lgd7->AddEntry(hthetamu,"#theta_{#mu}");}
    if(hthetaha!=0){hthetaha->SetFillColor(0);hthetaha->SetLineColor(4);hthetaha->SetMarkerColor(4);hthetaha->Draw("same");lgd7->AddEntry(hthetaha,"#theta_{had}");}
    if(hthetak2!=0){hthetak2->SetFillColor(0);hthetak2->SetLineColor(7);hthetak2->SetMarkerColor(7);hthetak2->Draw("same");lgd7->AddEntry(hthetak2,"#theta_{K} (no ISR)");}
    lgd7->Draw();
    TArrow *arr7 = new TArrow(179,0.5*hthetakamax,179,0,0.05);
    arr7->SetLineColor(6);
    arr7->SetLineWidth(2);
    arr7->Draw();
    int Nbin = hthetaee->GetNbinsX();
    double Noka7 = hthetak2->Integral(1,Nbin);
    double surv7 = hthetak2->Integral(1,(int)((179.-170.)/(180-170)*Nbin));
    std::cout<<"theta cut, Survive ratio is "<<1-hthetaka->Integral(1,(int)0.9*Nbin)/hthetaka->Integral(1,Nbin)<<", noISR "<<1-surv7/Noka7<<std::endl;
    
    TCanvas *c8 = new TCanvas();
    TLegend *lgd8 = new TLegend(0.15, 0.6,0.4,0.85);
    double hdtofkamax = hdtofka->GetMaximum();
    if(hdtofka!=0){hdtofka->SetFillColor(0);hdtofka->SetLineColor(1);hdtofka->SetMarkerColor(1);hdtofka->Draw();lgd8->AddEntry(hdtofka,"#delta TOF_{K}");} else continue;
    if(hdtofee!=0){hdtofee->SetFillColor(0);hdtofee->SetLineColor(2);hdtofee->SetMarkerColor(2);hdtofee->Draw("same");lgd8->AddEntry(hdtofee,"#delta TOF_{e}");}
    if(hdtofmu!=0){hdtofmu->SetFillColor(0);hdtofmu->SetLineColor(3);hdtofmu->SetMarkerColor(3);hdtofmu->Draw("same");lgd8->AddEntry(hdtofmu,"#delta TOF_{#mu}");}
    if(hdtofha!=0){hdtofha->SetFillColor(0);hdtofha->SetLineColor(4);hdtofha->SetMarkerColor(4);hdtofha->Draw("same");lgd8->AddEntry(hdtofha,"#delta TOF_{had}");}
    if(hdtofk2!=0){hdtofk2->SetFillColor(0);hdtofk2->SetLineColor(7);hdtofk2->SetMarkerColor(7);hdtofk2->Draw("same");lgd8->AddEntry(hdtofk2,"#delta TOF_{K} (no ISR)");}
    lgd8->Draw();  
    TArrow *arr8 = new TArrow(-3,0.5*hdtofkamax,-3,0,0.05);
    arr8->SetLineColor(6);
    arr8->SetLineWidth(2);
    arr8->Draw();
    TArrow *arr82 = new TArrow(3,0.5*hdtofkamax,3,0,0.05);
    arr82->SetLineColor(6);
    arr82->SetLineWidth(2);
    arr82->Draw();
    double Noka8 = hdtofk2->Integral(1,1000);
    double surv8 = hdtofk2->Integral(1,(int)((179.-170.)/(180-170)*1000));
    std::cout<<"tof cut, Survive ratio is "<<hdtofka->Integral(1,900)/hdtofka->Integral(1,1000)<<", noISR "<<surv8/Noka8<<std::endl;

    sprintf(filename,"%s_1p",file[enei]);
    outfile->WriteTObject(c1,filename);
    sprintf(filename,"output_cut/%s_1p.pdf",file[enei]);
    c1->SetLogy();
    c1->Print(filename);
    sprintf(filename,"%s_2p",file[enei]);
    outfile->WriteTObject(c2,filename);
  //sprintf(filename,"output_cut/%s_2p.pdf",file[enei]);
  //c2->Print(filename);
    sprintf(filename,"%s_ep",file[enei]);
    outfile->WriteTObject(c3,filename);
    sprintf(filename,"output_cut/%s_ep.pdf",file[enei]);
    c3->SetLogy();
    c3->Print(filename);
    sprintf(filename,"%s_2ep",file[enei]);
    outfile->WriteTObject(c4,filename);
  //sprintf(filename,"output_cut/%s_2ep.pdf",file[enei]);
  //c4->Print(filename);
    sprintf(filename,"%s_cos1",file[enei]);
    outfile->WriteTObject(c5,filename);
    sprintf(filename,"output_cut/%s_cos1.pdf",file[enei]);
    c5->Print(filename);
    sprintf(filename,"%s_cos2",file[enei]);
    outfile->WriteTObject(c6,filename);
    sprintf(filename,"output_cut/%s_cos2.pdf",file[enei]);
    c6->Print(filename);
    sprintf(filename,"%s_cos",file[enei]);
    outfile->WriteTObject(c7,filename);
    sprintf(filename,"output_cut/%s_cos.pdf",file[enei]);
    c7->SetLogy();
    c7->Print(filename);
    sprintf(filename,"%s_dtof",file[enei]);
    outfile->WriteTObject(c8,filename);
    sprintf(filename,"output_cut/%s_dtof.pdf",file[enei]);
    c8->Print(filename);
   
  //delete filei1;
  //delete filei2;
  //delete filei3;
  //delete filei4;
  //delete c1;
  //delete c2;
  //delete c3;
  //delete c4;
  //delete c5;
  //delete c6;
  //delete c7;
  //delete c8;
  //delete lgd1;
  //delete lgd2;
  //delete lgd3;
  //delete lgd4;
  //delete lgd5;
  //delete lgd6;
  //delete lgd7;
  //delete lgd8;
    //cutpar<<std::endl;
  }

  outfile->Close(); 
}

#include <iostream>
#include <fstream>
#include <sstream>
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAxis.h"
#include "TLegend.h"

int compareGen()
{
  ifstream oldgen("inicut.txt.oldmceff");
  ifstream newgen("inicut.txt");
  istringstream iss;
  iss.clear();
  
  // old generator
  double ene[100], totNo[100], gTrk[100], etgood[100];
  double cos[100], the[100], ep1[100], ep2[100], tof[100], p2g[100], nsur[100];
  double isr[100], ndata[100], nedata[100], lum[100];
  double eff[100], cor[100], cross[100], croe[100];
  int np=0;
  char line[1000];
  
  oldgen.getline(line,1000);
  while (!oldgen.eof()){
    oldgen.getline(line,1000);
    if (line[0]=='\0'||line[0]==' '|| line[0]=='#') continue;
    iss.clear();
    iss.str(line);
    iss>>ene[np]>>totNo[np]>>gTrk[np]>>etgood[np];
    iss>>cos[np]>>the[np]>>ep1[np]>>ep2[np]>>tof[np]>>p2g[np]>>nsur[np];
    iss>>isr[np]>>ndata[np]>>nedata[np]>>lum[np];
    gTrk[np] = gTrk[np]/totNo[np];
    eff[np] = nsur[np]/totNo[np];
    cor[np] = eff[np]*isr[np];
    cross[np] = ndata[np]/(lum[np]*eff[np]*isr[np]);
    croe[np] = nedata[np]/(lum[np]*eff[np]*isr[np]);
    np++;
  }


  // new generator
  double ene_new[100], totNo_new[100], gTrk_new[100], etgood_new[100];
  double cos_new[100], the_new[100], ep1_new[100], ep2_new[100], tof_new[100], p2g_new[100], nsur_new[100];
  double isr_new[100], ndata_new[100], nedata_new[100], lum_new[100];
  double eff_new[100], cor_new[100], cross_new[100], croe_new[100];
  int np_new=0;
  //char line[1000];
  
  newgen.getline(line,1000);
  while (!newgen.eof()){
    newgen.getline(line,1000);
    if (line[0]=='\0'||line[0]==' '|| line[0]=='#') continue;
    iss.clear();
    iss.str(line);
    iss>>ene_new[np_new]>>totNo_new[np_new]>>gTrk_new[np_new]>>etgood_new[np_new];
    iss>>cos_new[np_new]>>the_new[np_new]>>ep1_new[np_new]>>ep2_new[np_new]>>tof_new[np_new]>>p2g_new[np_new]>>nsur_new[np_new];
    iss>>isr_new[np_new]>>ndata_new[np_new]>>nedata_new[np_new]>>lum_new[np_new];
    gTrk_new[np_new] = gTrk_new[np_new]/totNo_new[np_new];
    eff_new[np_new] = nsur_new[np_new]/totNo_new[np_new];
    cor_new[np_new] = eff_new[np_new]*isr_new[np_new];
    cross_new[np_new] = ndata_new[np_new]/(lum_new[np_new]*eff_new[np_new]*isr_new[np_new]);
    croe_new[np_new] = nedata_new[np_new]/(lum_new[np_new]*eff_new[np_new]*isr_new[np_new]);
    np_new++;
  }

  cout<<np<<" "<<np_new<<endl;
  double eneerr[100] = {0};
  TCanvas *c1 = new TCanvas("c1","gen comparasion");
  c1->Divide(2,2);
  c1->cd(1);
  TGraphErrors *goldxs = new TGraphErrors(np,ene,cross,eneerr,croe);
  goldxs->SetMarkerStyle(25);
  goldxs->SetMarkerColor(2);
  goldxs->SetFillColor(0);
  goldxs->SetTitle("cross section");
  goldxs->GetXaxis()->SetTitle("#sqrt{s} GeV");
  goldxs->GetYaxis()->SetTitle("#sigma pb");
  TGraphErrors *gnewxs = new TGraphErrors(np_new,ene_new,cross_new,eneerr,croe_new);
  gnewxs->SetMarkerStyle(26);
  gnewxs->SetMarkerColor(3);
  gnewxs->SetFillColor(0);
  TLegend *legxs = new TLegend(0.55,0.5,0.8,0.8);
  legxs->AddEntry(goldxs,"with ConExc v2-88");
  legxs->AddEntry(gnewxs,"with ConExc v3-18");
  //TCanvas *c1 = new TCanvas("cxs","cross section");
  goldxs->Draw("AP");
  gnewxs->Draw("P");
  legxs->Draw();

  c1->cd(2);
  TGraph *goldeff = new TGraph(np,ene,eff);
  goldeff->SetMarkerStyle(25);
  goldeff->SetMarkerColor(2);
  goldeff->SetFillColor(0);
  goldeff->SetTitle("efficiency");
  goldeff->GetXaxis()->SetTitle("#sqrt{s} GeV");
  goldeff->GetYaxis()->SetTitle("#epsilon");
  TGraph *gneweff = new TGraph(np_new,ene_new,eff_new);
  gneweff->SetMarkerStyle(26);
  gneweff->SetMarkerColor(3);
  gneweff->SetFillColor(0);
  TLegend *legeff = new TLegend(0.55,0.5,0.8,0.8);
  legeff->AddEntry(goldeff,"with ConExc v2-88");
  legeff->AddEntry(gneweff,"with ConExc v3-18");
  //TCanvas *c2 = new TCanvas("ceff","efficiency");
  goldeff->Draw("AP");
  gneweff->Draw("P");
  legeff->Draw();
  
  c1->cd(3);
  TGraph *goldisr = new TGraph(np,ene,isr);
  goldisr->SetMarkerStyle(25);
  goldisr->SetMarkerColor(2);
  goldisr->SetFillColor(0);
  goldisr->SetTitle("ISR factor");
  goldisr->GetXaxis()->SetTitle("#sqrt{s} GeV");
  goldisr->GetYaxis()->SetTitle("1+#delta");
  TGraph *gnewisr = new TGraph(np_new,ene_new,isr_new);
  gnewisr->SetMarkerStyle(26);
  gnewisr->SetMarkerColor(3);
  gnewisr->SetFillColor(0);
  TLegend *legisr = new TLegend(0.55,0.5,0.8,0.8);
  legisr->AddEntry(goldisr,"with ConExc v2-88");
  legisr->AddEntry(gnewisr,"with ConExc v3-18");
  //TCanvas *c3 = new TCanvas("cisr","correction factor");
  goldisr->Draw("AP");
  gnewisr->Draw("P");
  legisr->Draw();
  
  c1->cd(4);
  TGraph *goldcor = new TGraph(np,ene,cor);
  goldcor->SetMarkerStyle(25);
  goldcor->SetMarkerColor(2);
  goldcor->SetFillColor(0);
  goldcor->SetTitle("#epsilon*(1+#delta)");
  goldcor->GetXaxis()->SetTitle("#sqrt{s} GeV");
  goldcor->GetYaxis()->SetTitle("#epsilon*(1+#delta)");
  TGraph *gnewcor = new TGraph(np_new,ene_new,cor_new);
  gnewcor->SetMarkerStyle(26);
  gnewcor->SetMarkerColor(3);
  gnewcor->SetFillColor(0);
  TLegend *legcor = new TLegend(0.55,0.5,0.8,0.8);
  legcor->AddEntry(goldcor,"with ConExc v2-88");
  legcor->AddEntry(gnewcor,"with ConExc v3-18");
  //TCanvas *c4 = new TCanvas("ccor","coriciency");
  goldcor->Draw("AP");
  gnewcor->Draw("P");
  legcor->Draw();
  
  TGraph *goldtrk = new TGraph(np,ene,gTrk);
  goldtrk->SetMarkerStyle(25);
  goldtrk->SetMarkerColor(2);
  goldtrk->SetFillColor(0);
  goldtrk->SetTitle("2 good track selection efficiency");
  goldtrk->GetXaxis()->SetTitle("#sqrt{s} GeV");
  goldtrk->GetYaxis()->SetTitle("#epsilon_{2 track}");
  TGraph *gnewtrk = new TGraph(np_new,ene_new,gTrk_new);
  gnewtrk->SetMarkerStyle(26);
  gnewtrk->SetMarkerColor(3);
  gnewtrk->SetFillColor(0);
  TLegend *legtrk = new TLegend(0.55,0.5,0.8,0.8);
  legtrk->AddEntry(goldtrk,"with ConExc v2-88");
  legtrk->AddEntry(gnewtrk,"with ConExc v3-18");
  TCanvas *c4 = new TCanvas("ctrk","trk");
  goldtrk->Draw("AP");
  gnewtrk->Draw("P");
  legtrk->Draw();


  return 0;
}

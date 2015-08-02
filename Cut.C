{
  TFile *_file0 = new TFile("output/output3.root");
  TTree *mctree = (TTree*)mcKK_001_665->Get("vars");
  TTree *datatree = (TTree*)KK_20000->Get("vars");
  gStyle->SetOptStat(0);

  //TCut cut("abs(costheta1)<0.8 & abs(costheta2)<0.8 & ep1<0.85 & ep2<0.85 & costheta<-0.99");
  TCut cut("abs(costheta1)<0.8 & abs(costheta2)<0.8 & ep1<0.85 & ep2<0.85 & costheta<-0.99 & abs(tof1-tof2)<0.7");
  
//TH1D *htof = new TH1D("htof","htof",1000,-10,10);
//TH1D *htof2 = new TH1D("htof2","htof",1000,-10,10);
//htof2->SetLineColor(2);
//new TCanvas();
//datatree->Draw("tof1-tof2>>htof",cut);
//mctree->Draw("tof1-tof2>>htof2",cut);
//htof->Draw();
//htof2->Draw("same");

//TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
//legend->AddEntry(htof,"data");
//legend->AddEntry(htof2,"MC");
//legend->Draw();

//TArrow *line1 = new TArrow(-0.7,1000,-0.7,0);
//TArrow *line2 = new TArrow( 0.7,1000, 0.7,0);
//line1->Draw();
//line2->Draw();

  TH1D *hp2 = new TH1D("hp2","hp",1000,0,1.5);
  TH1D *hp1 = new TH1D("hp1","hp",1000,0,1.5);
  TH1D *hp0 = new TH1D("hp0","hp",1000,0,1.5);
  hp1->SetLineColor(1);
  hp2->SetLineColor(3);
  hp0->SetLineColor(2);
  new TCanvas();
  datatree->Draw("p1>>hp1",cut);
  TCut cut("abs(costheta1)<0.8 & abs(costheta2)<0.8 & ep1<0.85 & ep2<0.85 & costheta<-0.99 & abs(tof1-tof2)<0.7 & abs(p2-0.86965)<0.03");
  datatree->Draw("p1>>hp2",cut);
  mctree->Draw("p1>>hp0",cut);
  hp1->Draw();
  hp2->Draw("same");
  hp0->Draw("same");
  
  TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
  legend->AddEntry(hp1,"Before cut");
  legend->AddEntry(hp2,"After cut");
  legend->AddEntry(hp0,"MC sample");
  legend->Draw();
}

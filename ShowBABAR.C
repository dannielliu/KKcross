#include <sstream>
#include <fstream>
#include <iostream>
#include <math.h>
#include <iomanip>

#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMath.h"
#include "TF1.h"

using namespace std;
int EffAndCross(char* name=0);
double BreExp(double *x, double *par)
{
  return par[0]*TMath::BreitWigner(x[0],par[1],par[2])+par[3]/TMath::Power(x[0],6);
}

int main(int argc, char** argv)
{
  if (argc>1) EffAndCross(argv[1]);
  else EffAndCross();
  return 0;
}

int EffAndCross(char* name)
{
 
  char line[1000];
  TFile *file = new TFile("cross.root","recreate");
  
  double xx[1000];
  double xxe[1000];
  double yy[1000];
  double yye[1000];
  double emat[1000][1000];

  int n2=0;
  ifstream infile2("KKBabar.txt");
  if (!infile2.is_open()) return -100;
  while (!infile2.eof())
  { 
    double energy1, energy2;
    double obs, obse;

    infile2.getline(line,1000);
    istringstream iss;
    iss.str(line);
  //if (iss.peek()<48 || iss.peek()>57) // not a number
  //{
  //  continue;  
  //}
    //else 
    {
      char a;
      iss >> energy1 >> a >> energy2;
      iss >> obs >> obse;
    }
    //std::cout<<"E1: "<< energy1 << "\tE2: "<< energy2 << "\tsigma: "<< obs <<"\terr: "<<obse<<std::endl;

    xx[n2] = (energy2+energy1)/2.;
    xxe[n2] = (energy2-energy1)/2.;
    yy[n2] = obs;
    yye[n2] = 0;
    n2 ++;
  }

  ifstream infile3("epapsFileKpKm.txt");
  if (!infile3.is_open()) return -200;
  int xid=0,yid=0;
  while(!infile3.eof())
  {
    infile3.getline(line,1000);
    istringstream iss;
    iss.str(line);
    if ((iss.peek()<48 || iss.peek()>57) && iss.peek()!='-') // not a number
    {
      continue;  
    }
    iss >> emat[xid][yid];
    yid++;
    if (yid==159){
      xid++;
      yid=0;
    }
    if (xid==159) break;
  }
  //std::cout<<"xid: "<<xid<<"\tyid: "<<yid<<std::endl;
  for (int i=0; i<159;i++)
  {
    yye[i] = sqrt(emat[i][i]);
  }

  TCanvas *c2 = new TCanvas();
  TGraphErrors *graph1 = new TGraphErrors(n2,xx,yy,xxe,yye);
  graph1->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  graph1->GetYaxis()->SetTitle("#sigma (nb)");
  graph1->SetTitle("Cross section");
  //graph1->SetMarkerStyle(5);
  graph1->SetFillColor(0);
  graph1->Draw("AP");

  ofstream fitxs("fitxsbabar.txt");
  double gap = 0.01;
  TF1 *fit;
  for (int i=0; xx[i]<1.1-0.0000001;i++){
    double xs = yy[i];
    fitxs<<xx[i]<<"\t"<<xs<<"\t"<<0<<endl;
    //fitxs << "xx.push_back(" << xx[i] << ");\tyy.push_back("<< yy[i] << ");\t er.push_back(0);" <<endl;
  }
 
//fit = new TF1("f0","[0]/(pow(x-[1],2)+pow([2],2))",0.99,1.05);
//fit->SetParameter(0,1.0);
//fit->SetParameter(1,1.02);
//fit->SetParameter(2,0.005);
//graph1->Fit(fit,"R","",0.99,1.05);
//for (double x=0.99; x<1.05-0.0000001;x+=gap){
//  double xs = fit->Eval(x);
//  fitxs << "xx.push_back(" << x << ");\tyy.push_back("<< xs << ");\t er.push_back(0);" <<endl;
//}

  fit = new TF1("f1","pol2",1.1,1.6);
  graph1->Fit(fit,"R+","",1.1,1.6);
  for (double x=1.1; x<1.6-0.0000001;x+=gap){
    double xs = fit->Eval(x);
    fitxs<<x<<"\t"<<xs<<"\t"<<0<<endl;
    //fitxs << "xx.push_back(" << x << ");\tyy.push_back("<< xs << ");\t er.push_back(0);" <<endl;
  }
  
  fit = new TF1("f2","pol2",1.6,1.8);
  graph1->Fit(fit,"R+","",1.6,1.8);
  for (double x=1.6; x<1.8-0.0000001;x+=gap){
    double xs = fit->Eval(x);
    fitxs<<x<<"\t"<<xs<<"\t"<<0<<endl;
    //fitxs << "xx.push_back(" << x << ");\tyy.push_back("<< xs << ");\t er.push_back(0);" <<endl;
  }
  
  fit = new TF1("f3","pol2",1.8,2.16);
  graph1->Fit(fit,"R+","",1.8,2.16);
  for (double x=1.8; x<2.16-0.0000001;x+=gap){
    double xs = fit->Eval(x);
    fitxs<<x<<"\t"<<xs<<"\t"<<0<<endl;
    //fitxs << "xx.push_back(" << x << ");\tyy.push_back("<< xs << ");\t er.push_back(0);" <<endl;
  }

  fit = new TF1("f4","gaus",2.16,2.4);
  graph1->Fit(fit,"R+","",2.16,2.4);
  for (double x=2.16; x<2.4-0.0000001;x+=gap){
    double xs = fit->Eval(x);
    fitxs<<x<<"\t"<<xs<<"\t"<<0<<endl;
    //fitxs << "xx.push_back(" << x << ");\tyy.push_back("<< xs << ");\t er.push_back(0);" <<endl;
  }
  
  fit = new TF1("f5","expo",2.4,5);
  graph1->Fit(fit,"R+","",2.4,5);
  for (double x=2.4; x<5-0.0000001;x+=gap){
    double xs = fit->Eval(x);
    fitxs<<x<<"\t"<<xs<<"\t"<<0<<endl;
    //fitxs << "xx.push_back(" << x << ");\tyy.push_back("<< xs << ");\t er.push_back(0);" <<endl;
  }



//graph = new TGraph(poNo,x,y);
////graph->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
////graph->GetYaxis()->SetTitle("Cross section (nb)");
//graphe->SetMarkerStyle(5);
//graphe->SetLineColor(2);
//graphe->SetFillColor(0);
//graphe->SetMarkerColor(2);
//graphe->Draw("P");

//TLegend *legend = new TLegend(0.6,0.6,0.85,0.85);
//legend->AddEntry(graph1,"BaBar 2013");
//legend->AddEntry(graphe, "BESIII");
//legend->Draw();

  file->WriteTObject(c2);

// output BaBar KK cross section
//for (int i=0;i<n2;i++){
//  std::cout<<xx[i]<<"   "<<yy[i]<<"   "<<yye[i]<<std::endl;
//}

  return 0;
}

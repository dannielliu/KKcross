#include <sstream>
#include <fstream>
#include <iostream>
#include <math.h>

#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TAxis.h"

using namespace std;
int EffAndCross();

int main()
{
  EffAndCross();
  return 0;
}

int EffAndCross()
{
  char filename[1000]={"inicut.txt"};
  ifstream infile(filename);

  TGraph *graph;
  double x[100];
  double y[100];
  double xe[100];
  double ye[100];

  double yeff[100];
  double ycor[100];
  double yec[100];
  
  int poNo=0;

  char line[1000];
  std::cout<<"Energy " << "\tefficiency "<< "\tcross1 " << "\tcross2" <<std::endl;
  while (!infile.eof())
  { 
    double energy, totNo, noISR, cut2trk, cntISR, cutcos1, cutcos2, cutcos, cutep1, cutep2, cuttof, cutp;
    double np1, npdata, nerr;
    double isrcor, lum;
    //double obs1, obs2;
    //double nm, nmerr;

    infile.getline(line,1000);
    istringstream iss;
    iss.str(line);
    if (iss.peek()<48 || iss.peek()>57) // not a number
    {
      continue;  
    }
    else {
      iss >> energy >> totNo >> noISR >> cut2trk >> cntISR >> cutcos1 >> cutcos2 >> cutcos >> cutep1 >> cutep2 >> cuttof >> cutp >> np1 >>isrcor >> npdata >> nerr >> lum;
      //iss >> isrcor; 
      //char a; iss >> a >> a >> a; double err; iss >> err;
      //iss >> obs1 >> nm >> nmerr;
      //iss >> lum;
    }
    double eff = np1/totNo;
    double cross = npdata/(lum*eff*isrcor);
    double crosserr = nerr/(lum*eff*isrcor);
    std::cout<<energy<<"\t"<< eff << "\t"<< cross <<"\t"<< crosserr << "\t" << isrcor*eff <<std::endl;

    x[poNo] = energy;
    xe[poNo]=0;
    y[poNo] = cross;
    ye[poNo]= crosserr;
    yeff[poNo] = eff;
    ycor[poNo] = isrcor;
    yec[poNo]  = isrcor*eff;
    poNo ++;
  }
  TCanvas *c1 = new TCanvas();
  TGraphErrors* graphe = new TGraphErrors(poNo,x,y,xe,ye);
  graphe->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  graphe->GetYaxis()->SetTitle("Cross section (nb)");
  graphe->SetMarkerStyle(5);
  graphe->Draw("AP");

  TFile *file = new TFile("cross.root","recreate");
  file->WriteTObject(c1);
  
  graph = new TGraph(poNo,x,yeff);
  graph->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  graph->GetYaxis()->SetTitle("#epsilon");
  graph->SetMarkerStyle(5);
  graph->Draw("AP");

  file->WriteTObject(c1);
  
//graph = new TGraph(poNo,x,ycor);
//graph->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
//graph->GetYaxis()->SetTitle("1+#delta");
//graph->SetMarkerStyle(5);
//graph->Draw("AP");

//file->WriteTObject(c1);
//
  graph = new TGraph(poNo,x,yec);
  graph->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  graph->GetYaxis()->SetTitle("(1+#delta)*#epsilon");
  graph->SetTitle("(1+#delta)*#epsilon");
  graph->SetMarkerStyle(5);
  graph->Draw("AP");

  file->WriteTObject(c1);

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
    std::cout<<"E1: "<< energy1 << "\tE2: "<< energy2 << "\tsigma: "<< obs <<"\terr: "<<obse<<std::endl;

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
  std::cout<<"xid: "<<xid<<"\tyid: "<<yid<<std::endl;
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

//graph = new TGraph(poNo,x,y);
////graph->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
////graph->GetYaxis()->SetTitle("Cross section (nb)");
  graphe->SetMarkerStyle(5);
  graphe->SetLineColor(2);
  graphe->SetFillColor(0);
  graphe->SetMarkerColor(2);
  graphe->Draw("P");

  TLegend *legend = new TLegend(0.6,0.6,0.85,0.85);
  legend->AddEntry(graph1,"BaBar 2013");
  legend->AddEntry(graphe, "BESIII");
  legend->Draw();

  file->WriteTObject(c2);

// output BaBar KK cross section
  for (int i=0;i<n2;i++){
    std::cout<<xx[i]<<"   "<<yy[i]<<"   "<<yye[i]<<std::endl;
  }

  return 0;
}

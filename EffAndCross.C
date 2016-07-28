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
double BreE6(double *x, double *par)
{
  return par[0]*TMath::BreitWigner(x[0],par[1],par[2])+par[3]/TMath::Power(x[0],6);
}
double BreExp(double *x, double *par)
{
  return par[0]*TMath::BreitWigner(x[0],par[1],par[2])+par[3]*TMath::Exp(-pow(x[0]-par[4],1)/par[5]) + par[6] ;
}
double Breinvx(double *x, double *par)
{
  return par[0]*TMath::BreitWigner(x[0],par[1],par[2])+par[3]/pow(x[0]-par[4],1.5) + par[5] ;
}



int main(int argc, char** argv)
{
  if (argc>1) EffAndCross(argv[1]);
  else EffAndCross();
  return 0;
}

int EffAndCross(char* name)
{
  char filename[1000]={"inicut_new.txt"};
  if (name!=0) sprintf(filename,"%s",name);
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
  ofstream ofbes("cross.txt_665p01");
  ofbes<<"Energy \t TotNo \t fitp \t eff \t isr \t e*(1+δ) \t cross \t err \t MCstat \t p_cut \t Stat.err" <<std::endl;
  while (!infile.eof())
  { 
    double energy, totNo, noISR, cut2trk, cntISR, cutcos1, cutcos2, cutcos, cutep1, cutep2, cuttof, cutp;
    double np1, npdata, nerr;
    double isrcor, lum;
    //double obs1, obs2;
    //double nm, nmerr;
    double np1_cut2;

    infile.getline(line,1000);
    istringstream iss;
    iss.str(line);
    if (iss.peek()<48 || iss.peek()>57) // not a number
    {
      continue;  
    }
    else {
      iss >> energy >> totNo ;
      //iss >>/* cntISR >>*/ cutcos1 >> cutcos2 >> cutcos >> cutep1 >> cutep2 >> cuttof >> cutp ;
      iss >> np1 >>isrcor >> npdata >> nerr >> lum;
      //iss >> np1_cut2;; 
      //char a; iss >> a >> a >> a; double err; iss >> err;
      //iss >> obs1 >> nm >> nmerr;
      //iss >> lum;
    }
    double eff = np1/totNo;
    double cross = npdata/(lum*eff*isrcor)/1000;
    double crosserr = nerr/(lum*eff*isrcor)/1000;
    //std::cout<<setiosflags(ios::fixed);
    ofbes<<setprecision(6);
    ofbes<<energy<< "\t" <<totNo  << "\t" << np1;
    ofbes<<fixed;
    ofbes<<setprecision(4);
    ofbes<<"\t"<< eff << "\t" << isrcor << "\t" << isrcor*eff  << "\t"<< cross<<"\t"<< crosserr ;
    ofbes<<'\t'<< sqrt((1-eff)/(totNo*eff))*100<<'\t'<< fabs(np1_cut2-npdata)/npdata*100 <<'\t'<< crosserr/cross*100;// <<std::endl;
    ofbes<<defaultfloat;

    x[poNo] = energy;
    xe[poNo]=0;
    y[poNo] = cross;
    ye[poNo]= crosserr;
    yeff[poNo] = eff;
    ycor[poNo] = isrcor;
    yec[poNo]  = isrcor*eff;

    // form factor
    double mk = 0.493677;
    double beta = sqrt(1-4*pow(mk/energy,2));
    double alpha = 1./137;
    double fsq = 3*pow(energy,2)/(TMath::Pi()*pow(alpha,2)*pow(beta,3))*cross; // unit: GeV^2 * nb
    double fsqerr = 3*pow(energy,2)/(TMath::Pi()*pow(alpha,2)*pow(beta,3))*crosserr; // unit: GeV^2 * nb
    // convert GeV^2 * nb to natural unit, (hc)^2 = 0.389 379 338 GeV^2*mbar = 389 379.338 GeV^2*nb
    fsq = fsq/389379.338;
    fsqerr = fsqerr/389379.338;
    ofbes<< "\t"<< energy<<"\t"<< fsq << "\t"<< fsqerr<<endl;

    poNo ++;
    if (poNo>=22) break;
  }
  TCanvas *c1 = new TCanvas();
  TGraphErrors* graphe = new TGraphErrors(poNo,x,y,xe,ye);
  graphe->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  graphe->GetYaxis()->SetTitle("Cross section (nb)");
  graphe->SetMarkerStyle(5);
  graphe->Draw("AP");
  
  TF1 *fit = new TF1("fit",Breinvx,2.0,3.1,6);
  fit->SetParName(0,"Const1");
  fit->SetParName(1,"mean");
  fit->SetParName(2,"sigma");
  fit->SetParName(3,"Const2");
  fit->SetParName(4,"shiftx");
  fit->SetParName(5,"shifty");
  fit->SetParLimits(0,0.01,0.05);
  fit->SetParLimits(1,2.2,2.3);
  fit->SetParLimits(2,0.08,0.15);
  fit->SetParLimits(3,0.01,3);
  fit->SetParLimits(4,1.85,2.1);
  fit->SetParLimits(5,-0.05,0.05);
  fit->SetParameters(0.02, 2.25, 0.085, 0.05, 1.95, 0.0);
  graphe->Fit(fit,"","",2.02,3.1);
  c1->Print("cross2.pdf");
  return 0;

  TFile *file = new TFile("cross.root","recreate");
  file->WriteTObject(c1,"CrossSection");
  
  graph = new TGraph(poNo,x,yeff);
  graph->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  graph->GetYaxis()->SetTitle("#epsilon");
  graph->SetMarkerStyle(5);
  graph->Draw("AP");

  file->WriteTObject(c1,"Efficiency");
  
  graph = new TGraph(poNo,x,ycor);
  graph->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  graph->GetYaxis()->SetTitle("1+#delta");
  graph->SetMarkerStyle(5);
  graph->Draw("AP");

  file->WriteTObject(c1,"Cor");
//
  graph = new TGraph(poNo,x,yec);
  graph->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  graph->GetYaxis()->SetTitle("(1+#delta)*#epsilon");
  graph->SetTitle("(1+#delta)*#epsilon");
  graph->SetMarkerStyle(5);
  graph->Draw("AP");

  file->WriteTObject(c1,"EffCor");

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

  for (int i=0; i<159;i++)
  {
    cout<<xx[i]<<"\t"<<yy[i]<<"\t"<<yye[i]<<endl;
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
//for (int i=0;i<n2;i++){
//  std::cout<<xx[i]<<"   "<<yy[i]<<"   "<<yye[i]<<std::endl;
//}

  return 0;
}

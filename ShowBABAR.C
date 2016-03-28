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
  return par[0]*TMath::BreitWigner(x[0],par[1],par[2])+par[3]/TMath::Power(x[0]-par[4],6);
}

double LandauExp(double *x, double *par)
{
  return par[0]*TMath::Landau(x[0],par[1],par[2])+par[3]/TMath::Exp(par[4]*x[0]);
}

double GausInvx(double *x, double *par)
{
  return par[0]*TMath::Gaus(x[0],par[1],par[2])+par[3]/pow(x[0]-par[4],par[5]);
}
double DGausInvx(double *x, double *par)
{
  return par[0]*TMath::Gaus(x[0],par[1],par[2])+par[3]/pow(x[0]-par[4],par[5])
        +par[6]*TMath::Gaus(x[0],par[7],par[8]);
}

double DGausExp(double *x, double *par)
{
  return par[0]*TMath::Gaus(x[0],par[1],par[2])
        +par[3]*TMath::Gaus(x[0],par[4],par[5])
	+par[6]/TMath::Exp(par[7]*x[0]);
}

double TGaus(double *x, double *par)
{
  return par[0]*TMath::Gaus(x[0],par[1],par[2])
        +par[3]*TMath::Gaus(x[0],par[4],par[5])
        +par[6]*TMath::Gaus(x[0],par[7],par[8]);
}

double maxwell(double *x, double *par)
{
  double xx = x[0]-par[1];
  if (xx<0) return 0;
  else return par[0]*TMath::Exp(-pow(xx,1.5)*par[2])*pow(xx,2.);
}

double DGausAExp(double *x, double *par)
{
  return        maxwell(x,&par[0])
        +par[3]*TMath::Gaus(x[0],par[4],par[5])
        +par[6]*TMath::Gaus(x[0],par[7],par[8]);
}

double AGausDExp(double *x, double *par)
{
  return        maxwell(x,&par[0])
        +par[3]*TMath::Gaus(x[0],par[4],par[5])
        +par[6]*TMath::Exp(-(x[0]-par[7])*par[8]);
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


  ifstream inbescross("cross.txt_665p01");
  double xbes[22];
  double xebes[22];
  double ybes[22];
  double yebes[22];
  inbescross.getline(line,1000);
  int idbes=0;
  while (!inbescross.eof() && idbes<21){
    double tmpe, tmpeff, tmpcross, tmpcrosse, tmp;
    inbescross >> tmpe >> tmp >> tmp>> tmpeff >> tmp >> tmp >> tmpcross >> tmpcrosse >> tmp >> tmp >> tmp ;
    xbes[idbes] = tmpe;
    xebes[idbes] = 0;
    ybes[idbes] = tmpcross;
    yebes[idbes] = tmpcrosse;
    idbes++;
  }

  TGraphErrors *gbes = new TGraphErrors(idbes,xbes,ybes,xebes,yebes);
  gbes->SetLineColor(2);
  gbes->SetMarkerColor(2);
  gbes->SetFillColor(0);
  //gbes->Draw("P");








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

  //---------------------------
  //add bes result
  for (int i=0;i<21;i++){
    xx[n2] = xbes[i];
    xxe[n2] = 0.001;
    yy[n2] = ybes[i];
    yye[n2] = yebes[i];
    n2 ++; 
  }
  //----------------------------

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
    //fitxs<<xx[i]<<"\t"<<xs<<"\t"<<0<<endl;
    fitxs << "xx.push_back(" << xx[i] << ");\tyy.push_back("<< yy[i] << ");\t er.push_back(0);" <<endl;
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
    //fitxs<<x<<"\t"<<xs<<"\t"<<0<<endl;
    fitxs << "xx.push_back(" << x << ");\tyy.push_back("<< xs << ");\t er.push_back(0);" <<endl;
  }
  
  fit = new TF1("f2","pol2",1.6,1.76);
  graph1->Fit(fit,"R+","",1.6,1.779);
  for (double x=1.6; x<1.779-0.0000001;x+=gap){
    double xs = fit->Eval(x);
    //fitxs<<x<<"\t"<<xs<<"\t"<<0<<endl;
    fitxs << "xx.push_back(" << x << ");\tyy.push_back("<< xs << ");\t er.push_back(0);" <<endl;
  }
  
//fit = new TF1("f3","pol2",1.8,2.16);
//graph1->Fit(fit,"R+","",1.8,2.16);
//for (double x=1.8; x<2.16-0.0000001;x+=gap){
//  double xs = fit->Eval(x);
//  fitxs<<x<<"\t"<<xs<<"\t"<<0<<endl;
//  //fitxs << "xx.push_back(" << x << ");\tyy.push_back("<< xs << ");\t er.push_back(0);" <<endl;
//}
  
  
//
//fit = new TF1("f4","gaus",2.16,2.35);
//graph1->Fit(fit,"R+","",2.16,2.4);
//for (double x=2.16; x<2.4-0.0000001;x+=gap){
//  double xs = fit->Eval(x);
//  fitxs<<x<<"\t"<<xs<<"\t"<<0<<endl;
//  //fitxs << "xx.push_back(" << x << ");\tyy.push_back("<< xs << ");\t er.push_back(0);" <<endl;
//}
//
//fit = new TF1("f5","expo",2.35,5);
//graph1->Fit(fit,"R+","",2.4,5);
//for (double x=2.4; x<5-0.0000001;x+=gap){
//  double xs = fit->Eval(x);
//  fitxs<<x<<"\t"<<xs<<"\t"<<0<<endl;
//  //fitxs << "xx.push_back(" << x << ");\tyy.push_back("<< xs << ");\t er.push_back(0);" <<endl;
//}
  
//fit = new TF1("f4",GausInvx,2.,5,6);
//fit->SetNpx(1000);
//fit->SetParameters(1,2.25,0.05,0.03,1.8,1.5);
//fit->SetParLimits(0, 0.1,  10);
//fit->SetParLimits(1, 2.2,  2.4);
//fit->SetParLimits(2, 0.03, 0.05);
//fit->SetParLimits(3, 0.02,  10);
//fit->SetParLimits(4, 1.5,    2);
//fit->SetParLimits(5, 1,    2);
  
//fit = new TF1("f4",DGausInvx,2.,5,9);
//fit->SetNpx(1000);
//fit->SetParameters(0.2,   2.25, 0.05,
//                   0.03,  1.8,  5,
//                   0.1,   2.0,  0.05);
//
//fit->SetParLimits(0, 0.1,  10);
//fit->SetParLimits(1, 2.2,  2.3);
//fit->SetParLimits(2, 0.03, 0.05);
//
//fit->SetParLimits(3, 0.02, 10);
//fit->SetParLimits(4, -1,    2);
//fit->SetParLimits(5, 5,    20);
//
//fit->SetParLimits(6, 0.1,  0.2);
//fit->SetParLimits(7, 1.9,  2.05);
//fit->SetParLimits(8, 0.04, 0.08);
//graph1->Fit(fit,"R+","",1.8,5);

//fit = new TF1("f4",DGausExp,2.,5,8);
//fit->SetNpx(1000);
//fit->SetParameters(0.2,   2.25, 0.05,
//                   0.1,   2.0,  0.06,
//                   0.03,  1.8);
//
//fit->SetParLimits(0, 0.1,  10);
//fit->SetParLimits(1, 2.2,  2.4);
//fit->SetParLimits(2, 0.03, 0.052);
//
//fit->SetParLimits(3, 0.1,  1.0);
//fit->SetParLimits(4, 1.9,    2.0);
//fit->SetParLimits(5, 0.06,   0.1);
//
//fit->SetParLimits(6, 0.02,  9);
//fit->SetParLimits(7, 0,    1.8);
//
//graph1->Fit(fit,"R+","",1.75,5);
 
//fit = new TF1("f4",TGaus,2.,5,9);
//fit->SetNpx(1000);
//fit->SetParameters(0.2,   2.25, 0.05,
//                   0.1,   2.0,  0.06,
//                   0.03,  1.8,  1.0);
//
//fit->SetParLimits(0, 0.1,  10);
//fit->SetParLimits(1, 2.2,  2.4);
//fit->SetParLimits(2, 0.03, 0.052);
//
//fit->SetParLimits(3, 0.1,  1.0);
//fit->SetParLimits(4, 1.9,    2.0);
//fit->SetParLimits(5, 0.06,   0.1);
//
//fit->SetParLimits(6, 0.02,  90);
////fit->SetParLimits(7, -10,    2);
//fit->SetParLimits(7, -10,    2.4);
//fit->SetParLimits(8, 0.2,    2);
  
//fit = new TF1("f4",DGausAExp,2.,5,9);
//fit->SetNpx(1000);
//fit->SetParameters(0.5,   2.2,  30,
//                   0.1,   2.0,  0.06,
//                   0.03,  1.8,  1.0);
//
//fit->SetParLimits(0, 0.1,  100);
//fit->SetParLimits(1, 2.1,  2.3);
//fit->SetParLimits(2, 10,   60);
//
//fit->SetParLimits(3, 0.1,  1.0);
//fit->SetParLimits(4, 1.9,    2.0);
//fit->SetParLimits(5, 0.06,   0.1);
//
//fit->SetParLimits(6, 0.02,  90);
//fit->SetParLimits(7, -10,    2.4);
//fit->SetParLimits(8, 0.2,    2);
 
  fit = new TF1("f4",AGausDExp,2.,5,9);
  fit->SetNpx(1000);
  fit->SetParameters(50,   2.2,  30,
                     0.1,  2.0,  0.06,
                     0.1,  2.2,  1.8);
  
  fit->SetParLimits(0, 0.1,  100);
  fit->SetParLimits(1, 2.1,  2.3);
  fit->SetParLimits(2, 10,   60);
  
  fit->SetParLimits(3, 0.1,  1.0);
  fit->SetParLimits(4, 1.9,    2.1);
  fit->SetParLimits(5, 0.06,   0.1);
  
  fit->SetParLimits(6, 0.02,  90);
  fit->SetParLimits(7, -10,    2.4);
  fit->SetParLimits(8, 0.2,    2);
  
  graph1->Fit(fit,"R+","",1.77,5);



  for (double x=1.77; x<5-0.0000001;x+=gap){
    double xs = fit->Eval(x);
    //fitxs<<x<<"\t"<<xs<<"\t"<<0<<endl;
    fitxs << "xx.push_back(" << x << ");\tyy.push_back("<< xs << ");\t er.push_back(0);" <<endl;
  }

  gbes->Draw("P");


  file->WriteTObject(c2);
  graph1->GetXaxis()->SetRangeUser(1.6,3.5);
  graph1->GetYaxis()->SetRangeUser(0,0.5);
  c2->Print("cross.pdf");

  return 0;
}

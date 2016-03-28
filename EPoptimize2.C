#include <iostream>
#include <fstream>
#include <sstream>
#include "./TminFitGra.h"
#include "TF1.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
using namespace std;

double func1(double *x, double *par)
{
  double chp=par[0];
  double xx = x[0];
  double par5 = par[3]*par[0]+par[4] - par[1]*pow(par[0]-par[2],2);
  if (xx> chp) return par[1]*pow(xx-par[2],2)+par5;
  else return par[3]*xx+par[4];
}

int EPoptimize2()
{
  ifstream infile("epSNR2.dat");
  ofstream outfile("OptEP.dat",ios::app);

  istringstream iss;
  char line[1000];

  cout<<"aaaaaaaa"<<endl;
  const int nene=21;
  const int nep = 100;
  double epn[nene][nep];
  double epnE[nene][nep];
  double sigN[nene][nep];
  double sigE[nene][nep];
  double bckN[nene][nep];
  double bckE[nene][nep];
  double sumN[nene][nep];
  double sumE[nene][nep];
  double snr[nene][nep];
  double snrE[nene][nep];
  int iene=-1;
  int iep=-1;
  double enes[nene];
  double eneerrs_low[nene];
  double eneerrs_up[nene];
  double opteps[nene];
  double opteperrs_low[nene];
  double opteperrs_up[nene];
  char title[nene][1000];

  bool start=0;
  if (!infile.is_open()) return -1;
  while (!infile.eof()){
    //cout<<"b"<<endl;
    infile.getline(line,1000);
	if (line[0]!='K') 
	{
	  if (!start) continue;
	  if (line[0]=='#' || line[0]==' ' || line[0]=='\0') continue;
	}
	else {
	  iep=-1; start=1;
	  iene++;
	  if (iene>=nene) break;
	  string tit;
	  iss.clear();
	  iss.str(line);
              iss >> tit;
	  std::cout<<tit<<endl;
	  sprintf(title[iene],"%s",tit.c_str());
	  continue;
	}
	iep++;
	//cout<<line<<endl;
	iss.clear();
	iss.str(line);
	double ratio,signal,signE,background,backE,tot,totE;
	iss >> ratio >> signal >> signE >> background >> backE; // 3 sigma
	iss >> tot >> totE;
	//iss >> signal >> background; // 5 sigma
        //cout<<"c "<< ratio<<endl;

	epn[iene][iep] = ratio;
	epnE[iene][iep] = 0;
	sigN[iene][iep] = signal;
	sigE[iene][iep] = signE;
	bckN[iene][iep] = background;
	bckE[iene][iep] = backE;
	sumN[iene][iep] = tot;
	sumE[iene][iep] = totE;

	
	if (sigN[iene][iep] < 0.1) {snr[iene][iep] = 0; snrE[iene][iep]=0;}
	else {
	  snr[iene][iep] = sigN[iene][iep]/sqrt(sigN[iene][iep]+bckN[iene][iep]);
	  snrE[iene][iep] = sqrt(pow(signE,2)/tot+pow(signal*totE,2)/pow(tot,3));
	}
    //cout<<"d"<<endl;
	//cout<< epn[iene][i]<<sigN[iene][i]<<bckN[iene][i]<<endl;
  }

  cout<<"e"<<endl;
  //for (iene=0;iene<1;iene++){
  for (iene=0;iene<nene;iene++){
	TGraphErrors *graph = new TGraphErrors(nep,epn[iene],snr[iene],epnE[iene],snrE[iene]);
	graph->SetTitle(title[iene]);
	graph->SetMarkerStyle(25);
	graph->GetYaxis()->SetTitle("S/#sqrt{S+N}");
	graph->GetXaxis()->SetTitle("E/p");
	TCanvas *ca1 = new TCanvas(title[iene], title[iene]);
	graph->Draw("AP");
	cout<<"a  np "<<graph->GetN()<<endl;
	ClearGraph(graph);
	cout<<"b  np "<<graph->GetN()<<endl;
	double max_y = GetMaximum(graph);
	cout<<"b  maximum "<<max_y<<endl;
	graph->GetYaxis()->SetRangeUser(0,1.2*max_y);
	//graph->Draw("AP");
	
	  graph->Fit("pol5","","",0.4,0.9);
	
////////////TminFitGra &mingra = gTmin;
////////////mingra.SetGraph(graph);
            TF1 *f1 = new TF1("f1","pol5",0.4,0.95);
            f1->SetParameters(graph->GetFunction("pol5")->GetParameters());
////////////TF1 *f1 = new TF1("f1",func1,0.4,0.99,5);
////////////f1->SetParameters(0.8, -1, 0.8, 50, 10);
////////////mingra.setFunc(f1);
//////////////mingra.showdata();
////////////mingra.fit();
	
	f1->SetLineColor(kGreen);
	f1->Draw("same");
	//graph->Fit("pol5","","",0.4,0.9);
	double max = f1->GetMaximum();
	double aveerr = GetAveErr(graph);
	double optep = f1->GetMaximumX(0.4,0.95);
	double optep_min = f1->GetX(max-aveerr,0.3, optep);
	double optep_max = f1->GetX(max-aveerr, optep, 0.95);
	outfile<<"Optmized E/p at "<< title[iene] << " is "<<optep << " : " << optep_min<<" - "<< optep_max<<std::endl;
	enes[iene] = atof(&title[iene][3])/10000.;
	opteps[iene] = optep;
	opteperrs_low[iene] = 0;
	opteperrs_up[iene] = 0;
	opteperrs_low[iene] = optep-optep_min;
	opteperrs_up[iene] = optep_max-optep;
	
	char ofname[1000];
	sprintf(ofname,"EPopt_%s.pdf",title[iene]);
	ca1->Print(ofname);
        
  }
  //return 1;
    cout<<"f"<<endl;
  TGraphAsymmErrors* gep = new TGraphAsymmErrors(iene,enes,opteps,eneerrs_low, eneerrs_up, opteperrs_low, opteperrs_up);
  //gep->SetTitle("Optimized E/p");
  gep->SetTitle("");
  gep->SetMarkerStyle(25);
  gep->GetYaxis()->SetTitle("E/p");
  gep->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
  TCanvas *ca1 = new TCanvas("gep","Optimized E/p");
  ca1->SetMargin(0.15,0.1,0.15,0.1);
  gep->GetXaxis()->SetTitleSize(0.05);
  gep->GetXaxis()->SetLabelSize(0.05);
  gep->GetXaxis()->SetTitleOffset(1.2);
  gep->GetYaxis()->SetTitleSize(0.05);
  gep->GetYaxis()->SetLabelSize(0.05);
  gep->GetYaxis()->SetTitleOffset(1.2);
  gep->Draw("AP");
  gep->Fit("pol2");
  ca1->Print("output/EPopt.png");
  double c0 = gep->GetFunction("pol2")->GetParameter(0);
  double c1 = gep->GetFunction("pol2")->GetParameter(1);
  double c2 = gep->GetFunction("pol2")->GetParameter(2);
  outfile << "Ep = "<< c0 << " + " << c1 << "*x + " << c2 << " *x^2"<<std::endl;
  std::cout << "Ep = "<< c0 << " + " << c1 << "*x + " << c2 << " *x^2"<<std::endl;
  
  return 0;
}


int main()
{
  EPoptimize2();
  return 0;
}

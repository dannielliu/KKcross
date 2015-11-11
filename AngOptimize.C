#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

int AngOptimize()
{
  ifstream infile("angSNR.dat");
  ofstream outfile("OptAng.dat",ios::app);

  istringstream iss;
  char line[1000];

  cout<<"aaaaaaaa"<<endl;
  const int nene=21;
  const int nep = 28;
  double epn[nene][nep];
  double sigN[nene][nep];
  double bckN[nene][nep];
  double snr[nene][nep];
  int iene=-1;
  int iep=-1;
  double enes[nene];
  double opteps[nene];
  char title[nene][1000];

  TF1 *fit = new TF1("fit","[0]*(x-[2])+[1]/(x-[2])+[3]",178,179.6);
  fit->SetParameter(0,1.5);
  fit->SetParameter(1,1);
  fit->SetParameter(2,180);
  fit->SetParLimits(2,180,180)
  fit->SetParameter(3,40);

  bool start=0;
  if (!infile.is_open()) return -1;
  while (!infile.eof()){
    //cout<<"b"<<endl;
    infile.getline(line,1000);
	if (line[0]!='K') 
	{
	  if (!start) continue;
	  if (line[0]=='#' || line[0]==' ') continue;
	}
	else {
	  iep=-1; start=1;
	  iene++;
	  if (iene>=21) break;
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
	double ratio,signal,background;
	iss >> ratio >> signal >> background; // 3 sigma
	//iss >> signal >> background; // 5 sigma

	epn[iene][iep] = ratio;
	sigN[iene][iep] = signal;
	bckN[iene][iep] = background;
	
	if (sigN[iene][iep]==0) snr[iene][iep] = 0;
	else snr[iene][iep] = sigN[iene][iep]/sqrt(sigN[iene][iep]+bckN[iene][iep]);
	//cout<< epn[iene][i]<<sigN[iene][i]<<bckN[iene][i]<<endl;
  }

  for (iene=0;iene<nene;iene++){
	TGraph *graph = new TGraph(nep,epn[iene],snr[iene]);
	graph->SetTitle(title[iene]);
	graph->SetMarkerStyle(25);
	graph->GetYaxis()->SetTitle("S/#sqrt{S+N}");
	graph->GetXaxis()->SetTitle("#theta");
	TCanvas *ca1 = new TCanvas(title[iene], title[iene]);
	graph->Draw("AP");
	//graph->Fit("pol4","","",178,179.9);
////////TF1 *fpol2 = graph->GetFunction("pol2");
////////fit->SetParameter(0,fpol2->GetParameter(0));
////////fit->SetParameter(1,fpol2->GetParameter(1));
////////fit->SetParameter(2,fpol2->GetParameter(2));
	//graph->Fit(fit,"","",178,179.6);
	graph->Fit(fit,"","",172.0+iene*6.0/21.0,179.9);
	double optep = fit->GetMaximumX(176,179.6);
	outfile<<"Optmized theta at "<< title[iene] << " is "<<optep<<std::endl;
	enes[iene] = atof(&title[iene][3])/10000.;
	opteps[iene] = optep;
	char ofname[1000];
	sprintf(ofname,"output/Angopt_%s.pdf",title[iene]);
	ca1->Print(ofname);
  } 
  TGraph* gep = new TGraph(iene,enes,opteps);
  gep->SetTitle("Optimized #theta");
  gep->SetMarkerStyle(25);
  gep->GetYaxis()->SetTitle("S/#sqrt{S+N}");
  gep->GetXaxis()->SetTitle("#theta");
  TCanvas *ca1 = new TCanvas("gep","Optimized #theta");
  gep->Draw("AP");
  gep->Fit("pol1");
  ca1->Print("output/Angopt.pdf");
  double c0 = gep->GetFunction("pol1")->GetParameter(0);
  double c1 = gep->GetFunction("pol1")->GetParameter(1);
  //double c2 = gep->GetFunction("pol2")->GetParameter(2);
  outfile << "Theta = "<< c0 << " + " << c1 << "*x + "  <<std::endl;
  std::cout << "Theta = "<< c0 << " + " << c1 << "*x + "<<std::endl;
  return 0;
}

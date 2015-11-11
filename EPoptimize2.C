#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

int EPoptimize2()
{
  ifstream infile("epSNR_eperr.dat");
  ofstream outfile("OptEP.dat",ios::app);

  istringstream iss;
  char line[1000];

  cout<<"aaaaaaaa"<<endl;
  const int nene=1;
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
  double opteps[nene];
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
  for (iene=0;iene<nene;iene++){
	TGraphErrors *graph = new TGraphErrors(nep,epn[iene],snr[iene],epnE[iene],snrE[iene]);
	graph->SetTitle(title[iene]);
	graph->SetMarkerStyle(25);
	graph->GetYaxis()->SetTitle("S/#sqrt{S+N}");
	graph->GetXaxis()->SetTitle("E/p");
	TCanvas *ca1 = new TCanvas(title[iene], title[iene]);
	graph->Draw("AP");
	graph->Fit("pol5","","",0.4,0.9);
	double optep = graph->GetFunction("pol5")->GetMaximumX(0.4,0.95);
	outfile<<"Optmized E/p at "<< title[iene] << " is "<<optep<<std::endl;
	enes[iene] = atof(&title[iene][3])/10000.;
	opteps[iene] = optep;
	char ofname[1000];
	sprintf(ofname,"output/EPopt_%s.pdf",title[iene]);
	ca1->Print(ofname);
        
  }
    cout<<"f"<<endl;
  TGraph* gep = new TGraph(iene,enes,opteps);
  gep->SetTitle("Optimized E/p");
  gep->SetMarkerStyle(25);
  gep->GetYaxis()->SetTitle("S/#sqrt{S+N}");
  gep->GetXaxis()->SetTitle("E/p");
  TCanvas *ca1 = new TCanvas("gep","Optimized E/p");
  gep->Draw("AP");
  gep->Fit("pol2");
  ca1->Print("output/EPopt.pdf");
  double c0 = gep->GetFunction("pol2")->GetParameter(0);
  double c1 = gep->GetFunction("pol2")->GetParameter(1);
  double c2 = gep->GetFunction("pol2")->GetParameter(2);
  outfile << "Ep = "<< c0 << " + " << c1 << "*x + " << c2 << " *x^2"<<std::endl;
  std::cout << "Ep = "<< c0 << " + " << c1 << "*x + " << c2 << " *x^2"<<std::endl;
  
  return 0;
}

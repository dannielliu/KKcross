//##############################
//### this code is modified ####
//### dliu13@mail.ustc.edu.cn###
//##############################
//##############################
//##############################


#include "TFile.h"
//#include "TFolder.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TLorentzVector.h"
#include "TH2D.h"
#include "TH1.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

using namespace std;
int Ana(const char *filename=0, const char* outdir=0, TFile* fileout=0);
double GetEnergy(int run);
double FitSpectrum(TTree *&dataraw, double beame, const char* namesfx=0); // cb signal, gaus back
double FitSpectrum1(TTree *&dataraw, double beame, const char* namesfx=0); // cb signal, poly bg
double FitSpectrum2(TTree *&dataraw, double beame, const char* namesfx=0);
double FitSpectrum3(TTree *&dataraw, double beame, const char* namesfx=0);
double FitSpectrum4(TTree *&dataraw, double beame, const char* namesfx=0, TFile* fileout=0); // fit with MC shape
double FitSpectrum5(TTree *&dataraw, double beame, const char* namesfx=0); // fix parameters
const char* getPureName(const char* name);
const char* getPureName2(const char* name);
double getEne(const char* name);


int main(int argc, char** argv)
{
  //if (argc==1) Ana(argv[1]);
  //if (argc>1)  Ana(argv[1],argv[2]);
  //else Ana();
  TFile *fileout = new TFile("output/output_mcang.root","update");
  for (int i=1; i<argc; i++){
    Ana(argv[i],0,fileout);
  }
  //fileout->Write();
  delete fileout;
  return 0;
}

int Ana(const char *filename, const char* outdir, TFile *fileout)
{
  //TFile *file = new TFile("mc_KPI_22324.root");
  TFile *file;
  if (filename==0) file= new TFile("KK_22324.root");
  else file = new TFile(filename);
  std::cout<<"File name is "<<file->GetName()<<std::endl;
  if (file==0) return -1;
  TTree *tree = (TTree*)file->Get("TwoProng");
  if (tree==0) return -2;

  double ene[22];
  double pcut[22];
  double epcut[22];
  double thecut[22];
  double m_pcut;
  double m_epcut;
  double m_thecut;
  ifstream cutin("cutpar");
  char line[1000];
  if (cutin.is_open()) cutin.getline(line,1000);
  int iene=0;
  // set E/p and p cut, mark: set cut
  while (!cutin.eof()){
    cutin.getline(line,1000);
    istringstream iss;
    iss.str(line);
    iss>>ene[iene]>>pcut[iene]>>epcut[iene];
    iene++;
    if (iene==22) break;
  }
  for (int i=0;i<22;i++) epcut[i] = 1.64295 - 0.629622*ene[i] + 0.104755 *pow(ene[i],2);
  //for (int i=0;i<21;i++) thecut[i] = 173.946 + 1.74736*ene[i];
  for (int i=0;i<22;i++) thecut[i] = 179;

  double kappx,kappy,kappz,kampx,kampy,kampz;
  int nneu;
  int run;
  int idxmc;
  int pdgid[100];
  int motheridx[100];
  int emcstatusInt;
  short emcstatusShort;
  double emctrk1;
  double emctrk2;
  double epratio1;
  double epratio2;
  int ntof1;
  int ntof2;
  int tofl1[5];
  int tofl2[5];
  double tof1[5];
  double tof2[5];
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("indexmc",&idxmc);
  tree->SetBranchAddress("pdgid",pdgid);
  tree->SetBranchAddress("motheridx",motheridx);
  tree->SetBranchAddress("kappx",&kappx);
  tree->SetBranchAddress("kappy",&kappy);
  tree->SetBranchAddress("kappz",&kappz);
  tree->SetBranchAddress("kampx",&kampx);
  tree->SetBranchAddress("kampy",&kampy);
  tree->SetBranchAddress("kampz",&kampz);
  //tree->SetBranchAddress("costheta",&costheta);
  tree->SetBranchAddress("nneu",&nneu);
  if (strncmp(tree->GetBranch("emcstatus")->GetLeaf("emcstatus")->GetTypeName(),"Short_t",7)==0)
  	tree->SetBranchAddress("emcstatus",&emcstatusShort);
  else 
  	tree->SetBranchAddress("emcstatus",&emcstatusInt);
  tree->SetBranchAddress("epratio1",&epratio1);
  tree->SetBranchAddress("epratio2",&epratio2);
  tree->SetBranchAddress("emctrk1",&emctrk1);
  tree->SetBranchAddress("emctrk2",&emctrk2);
  tree->SetBranchAddress("ntof1",&ntof1);
  tree->SetBranchAddress("toflayer1",tofl1);
  tree->SetBranchAddress("tof1",tof1);
  tree->SetBranchAddress("ntof2",&ntof2);
  tree->SetBranchAddress("toflayer2",tofl2);
  tree->SetBranchAddress("tof2",tof2);
  double pi = TMath::Pi();
  double mka = 0.493677;
  //double mka = 0.13957;
  //double mka = 0.1057;
  
  double mckappx, mckappy, mckappz, mckape;
  double mckampx, mckampy, mckampz, mckame;
  tree->SetBranchAddress("mckappx",&mckappx);
  tree->SetBranchAddress("mckappy",&mckappy);
  tree->SetBranchAddress("mckappz",&mckappz);
  tree->SetBranchAddress("mckape" ,&mckape);
  tree->SetBranchAddress("mckampx",&mckampx);
  tree->SetBranchAddress("mckampy",&mckampy);
  tree->SetBranchAddress("mckampz",&mckampz);
  tree->SetBranchAddress("mckame" ,&mckame);


  tree->GetEntry(0);
  double Ebeam = GetEnergy(run);
  if (Ebeam<1.0) Ebeam = getEne(filename);
  //Ebeam = 3.08;

  for (int iene=0;iene<22;iene++){
    if (fabs(Ebeam-ene[iene])<1e-4) {
      m_pcut = pcut[iene];
      m_epcut = epcut[iene];
      m_thecut = thecut[iene];
      break;
    }
  }
  //m_epcut = m_epcut +0.05; // change E/p to determine uncertainty from this cut
  //m_thecut = 179+0.2; // uncertainty from theta cut
  //m_pcut = 0.956409;
  m_thecut = 179;
  m_epcut = 1.64295 - 0.629622*Ebeam + 0.104755 *pow(Ebeam,2);
  double m_tofcut = 3;
  //m_tofcut = 3-1;
  std::cout<<"pcut "<<m_pcut<<", E/p cut "<<m_epcut<<std::endl;

//char name[100];
//if (Ebeam>0){
//  if (outdir==0) sprintf(name,"output_%.5f.root",Ebeam);
//  else sprintf(name,"%s/output_%.5f.root",outdir,Ebeam);
//}
//else {
//  if (outdir==0) sprintf(name,"output_%s",getPureName(filename));
//  else sprintf(name,"%s/output_%s",outdir,getPureName(filename));
//  Ebeam = getEne(filename);
//}
//std::cout<<"output file name is "<<name<<std::endl;
  //TFile *fileout = new TFile(name,"recreate");
  const char *pureName = getPureName2(filename);
  std::cout<<"Pure Name: "<< pureName <<std::endl;
  char name1[1000];
  //sprintf(name1,"output/%s.root",pureName);
  sprintf(name1,"output/%s.root",pureName);
  TFile *dir = new TFile(name1,"recreate");
  TTree *vars = new TTree("vars","vars");

//fileout->cd();
//TDirectory *dir = fileout->GetDirectory(pureName);
//if (dir==0) dir = fileout->mkdir(pureName);;
//dir->cd();
//TTree *vars = new TTree("vars","vars");

  TLorentzVector kap,kam;
  double mass;

  double totp, costheta1, costheta2;
  double costheta;
  double p1,p2;
  double tof11=-1000,tof21=-1000;
  int tofid;
  int isrtag=0;
  int tof1tag=0;
  int tof2tag=0;
  int emctag=0;
  vars->Branch("run",&run,"run/I");
  vars->Branch("indexmc",&idxmc,"indexmc/I");
  vars->Branch("pdgid",pdgid,"pdgid[100]/I");
  vars->Branch("motheridx",motheridx,"motheridx[100]/I");
  vars->Branch("isrtag",&isrtag,"isrtag/I");
  vars->Branch("tof1tag",&tof1tag,"tof1tag/I");
  vars->Branch("tof2tag",&tof2tag,"tof2tag/I");
  vars->Branch("emctag",&emctag,"emctag/I");
  vars->Branch("mass",&mass,"mass/D");
  vars->Branch("costheta",&costheta,"costheta/D");
  vars->Branch("costheta1",&costheta1,"costheta1/D");
  vars->Branch("costheta2",&costheta2,"costheta2/D");
  vars->Branch("totp",&totp,"totp/D");
  vars->Branch("p1",&p1,"p1/D");
  vars->Branch("p2",&p2,"p2/D");
  vars->Branch("ep1",&epratio1,"ep1/D");
  vars->Branch("ep2",&epratio2,"ep2/D");
  vars->Branch("ed1",&emctrk1,"ed1/D");
  vars->Branch("ed2",&emctrk2,"ed2/D");
  vars->Branch("tof1",&tof11,"tof1/D");
  vars->Branch("tof2",&tof21,"tof2/D");

  TTree *datasel1 = new TTree("datasel1","datasel1");
  datasel1->Branch("x1",&p1,"x1/D");
  datasel1->Branch("mass",&mass,"mass/D");
//TTree *datasel2 = new TTree("datasel2","datasel2");
//datasel2->Branch("x1",&p2,"x1/D");
//datasel2->Branch("mass",&mass,"mass/D");
//const int nangle=0;
//TTree* datasel[nangle];
//double thetan[nangle];
//for (int i=0; i<10;i++){
//  char tmpname[100];
//  sprintf(tmpname,"datasel_%02d",i);
//  datasel[i] = new TTree(tmpname,tmpname);
//  datasel[i]->Branch("x1",&p1,"x1/D");
//  datasel[i]->Branch("mass",&mass,"mass/D");
//  thetan[i] = 170+i;
//}
//for (int i=10; i<nangle;i++){
//  char tmpname[100];
//  sprintf(tmpname,"datasel_%02d",i);
//  datasel[i] = new TTree(tmpname,tmpname);
//  datasel[i]->Branch("x1",&p1,"x1/D");
//  datasel[i]->Branch("mass",&mass,"mass/D");
//  thetan[i] = 178.+(180.-178.)/(nangle-10)*(i-10);
//}
//
  fileout->cd();

  int nentries = tree->GetEntries();
  double pexp = sqrt(pow(Ebeam/2,2)-pow(mka,2));
  std::cout<<"Total entries is " << nentries << std::endl;
  std::cout<<"Beam energy is " << Ebeam << std::endl;
  std::cout<<"Expected K momentum is " << pexp << std::endl;
  int count0=0,count1=0,count2=0,count3=0;
  int count[10]={0};
  int tagmc=0;

  TH1D *hmcangpi= (TH1D*)file->Get("thep");
  TH1D *hmcangmi= (TH1D*)file->Get("them");
  if (hmcangpi==0 || hmcangmi==0) return -1;
  TH1D *hmcangp = new TH1D("hmcangp" ,"cos#theta+",200,-1,1);
  TH1D *hmcangp2= new TH1D("hmcangp2","cos#theta+",200,-1,1);
  TH1D *heffp   = new TH1D("heffp"   ,"#epsilon+" ,200,-1,1);
  TH1D *heffp2  = new TH1D("heffp2"  ,"#epsilon+" ,200,-1,1);
  TH1D *hmcangm = new TH1D("hmcangm" ,"cos#theta-",200,-1,1);
  TH1D *hmcangm2= new TH1D("hmcangm2","cos#theta-",200,-1,1);
  TH1D *heffm   = new TH1D("heffm"   ,"#epsilon-" ,200,-1,1);
  TH1D *heffm2  = new TH1D("heffm2"  ,"#epsilon-" ,200,-1,1);


  for (int ien=0;ien<nentries;ien++){
    tree->GetEntry(ien);
    
     if ( run == 39582||run == 39671||run == 39723||run == 39783||run == 40208||run == 40300||run == 40308||run == 40462||run == 40526||run == 40946||run == 41099||run == 41200||run == 41283||run == 41408||run == 41416||run == 41436||run == 41445||run == 41471||run == 41728||run == 41818||run == 41902) continue;
     
    isrtag = 0;
    for (int j=0;j<idxmc;j++){
      if (pdgid[j] == 22 && motheridx[j]==j) { tagmc++; isrtag=1;}
    }

 // if (emcstatusShort!=0x3 && emcstatusInt!=0x3) continue;
 // if (ntof1<1) continue;
 // if (ntof2<1) continue;
 // for (int i=0;i<ntof1;i++){
 //   if (tofl1[i]==1) tof11=tof1[i];
 // }
 // for (int i=0;i<ntof2;i++){
 //   if (tofl2[i]==1) tof21=tof2[i];
 // }
    //if (tof11<-999 || tof21<-999) continue;
    count0++;

    kap.SetVectMag(TVector3(kappx,kappy,kappz),mka);
    kam.SetVectMag(TVector3(kampx,kampy,kampz),mka);
    
    if (isrtag==1) continue;
    TLorentzVector mckp, mckm;
    mckp.SetVectMag(TVector3(mckappx,mckappy,mckappz),mka);
    mckm.SetVectMag(TVector3(mckampx,mckampy,mckampz),mka);
    double mcangp = mckp.Vect().CosTheta();
    double mcangm = mckm.Vect().CosTheta();
    hmcangp->Fill(mcangp);
    hmcangm->Fill(mcangm);

    // angular information in lab coordinate
    costheta1 = kap.CosTheta();
    costheta2 = kam.CosTheta();
    
    // momentum information in cms coordinalte
    kap.Boost(-0.011,0,0);
    kam.Boost(-0.011,0,0);
    costheta = kap.Vect().Dot(kam.Vect())/(kap.Vect().Mag()*kam.Vect().Mag());
    double theta = kap.Vect().Angle(kam.Vect())/TMath::Pi()*180.;
    mass = (kap+kam).M();
    totp = (kap+kam).Rho();
    p1 = kap.Rho();
    p2 = kam.Rho();
 
    //vars->Fill();

    // select candidate, p1 dis
    // if (count0%100000==0) std::cout<<"epratio1 is "<< epratio1<<std::endl;
    //if (fabs(costheta1)>0.93) cout<< "ien: "<< ien << "\t cos1: "<< fabs(costheta1)<<endl;
    if (fabs(costheta1)<0.93 && costheta1 < 0.8) { count[4]++;
    if (fabs(costheta2)<0.93 && costheta2 > -0.8){ count[5]++;
    //if (fabs(costheta1)<0.93 ) { count[4]++;
    //if (fabs(costheta2)<0.93 ){ count[5]++;
    //if (costheta1 > 0.8) continue; // Kp in positive direction
    //if (costheta2 < -0.8) continue; // Km in negetive direction
    if (theta>m_thecut){ count[2]++;
    if (epratio1<m_epcut)
    { count[0]++;
      if (epratio2<m_epcut)
      { count[1]++;
        if (fabs(tof11-tof21)<m_tofcut){
          count3++;
          if (p2<m_pcut) {/*datasel2->Fill();*/ count2++; // p2 dis
            datasel1->Fill();
	if (p1<m_pcut){
	  hmcangp2->Fill(costheta1);
	  hmcangm2->Fill(costheta2);
	}
          }
        }
      }
    }
    }
    
//////////////////////////////
//  if (fabs(tof11-tof21)<3 && p2<m_pcut && epratio1<m_epcut && epratio2<m_epcut){
//        for (int anglei=0; anglei<nangle; anglei++) {
//           if (theta>thetan[anglei]) datasel[anglei]->Fill();
//        }
//  }
//////////////////////////////

    }
    }
  }
  
  // extract efficiency
  for (int i=1;i<201;i++){
    double mcselp = hmcangp->GetBinContent(i);
    double mcselpi= hmcangpi->GetBinContent(i);
    double mcselp2= hmcangp2->GetBinContent(i);
    double mcselm = hmcangm->GetBinContent(i);
    double mcselmi= hmcangmi->GetBinContent(i);
    double mcselm2= hmcangm2->GetBinContent(i);

    if (mcselpi>0) heffp->SetBinContent(i,mcselp/mcselpi);
    else           heffp->SetBinContent(i,0);
    if (mcselp >0) heffp2->SetBinContent(i,mcselp2/mcselp);
    else           heffp2->SetBinContent(i,0);
    if (mcselmi>0) heffm->SetBinContent(i,mcselm/mcselmi);
    else           heffm->SetBinContent(i,0);
    if (mcselm >0) heffm2->SetBinContent(i,mcselm2/mcselm);
    else           heffm2->SetBinContent(i,0);
  }
  
  //dir->Write();
  //
  hmcangpi->Write();
  hmcangmi->Write();
  hmcangp->Write();
  hmcangp2->Write();
  hmcangm->Write();
  hmcangm2->Write();
  heffp->Write();
  heffp2->Write();
  heffm->Write();
  heffm2->Write();
  //return 0;

  ofstream cutflow("cutflow2",std::ios::app);
  //ofstream cutflow("cutflow_cmp665andp01",std::ios::app);
  cutflow<<filename<<std::endl;
  std::cout<<"data selction 1 size is "<<datasel1->GetEntries()<<std::endl;
  //std::cout<<"data selction 2 size is "<<datasel2->GetEntries()<<" "<<count2<<std::endl;
  cutflow<<"Initial size      :"<<nentries<<std::endl;
  cutflow<<"Tagged ISR evt    :"<< tagmc  <<std::endl;
  cutflow<<"Tagged no ISR evt :"<< nentries - tagmc  <<std::endl;
  cutflow<<"Valid tof and emc :"<<count0<<std::endl;
  cutflow<<"After cos1<0.8    :"<<count[4]<<std::endl;
  cutflow<<"After cos2>-0.8   :"<<count[5]<<std::endl;
  cutflow<<"After theta cut   :"<<count[2]<<std::endl;
  cutflow<<"After ep1         :"<<count[0]<<std::endl;
  cutflow<<"After ep2         :"<<count[1]<<std::endl;
  //cutflow<<"After totp<0.05   :"<<count[3]<<std::endl;
  cutflow<<"After dtof<3      :"<<count3<<std::endl;
  //cutflow<<"p1-exp<0.08 for p1:"<<count1<<std::endl;
  cutflow<<"p2-exp<3 sigma    :"<<count2<<std::endl;
  //if (count2==0) return -3;

  return 0;
  //return 0;

  //double sigp1=FitSpectrum5(datasel1, Ebeam, pureName);
  //double sigp1=FitSpectrum4(datasel1, Ebeam, pureName, fileout);
  double sigp1=FitSpectrum3(datasel1, Ebeam, pureName);
  cutflow<<"p1 fit signal     :"<<sigp1 <<std::endl;
  //cutflow<<"m1 fit signal     :"<<sigm1 <<std::endl;
//cutflow<<"p2 fit signal     :"<<sigp2 <<std::endl;
//cutflow<<"m2 fit signal     :"<<sigm2 <<std::endl;

//char nameangle[100];
//for (int anglei=0; anglei<nangle; anglei++){
//  sprintf(nameangle,"%s_%.2f",&pureName[9],thetan[anglei]);
//  //sprintf(nameangle,"%.2f",epn[epi]);
//  std::cout<<nameangle<<" entries is " << datasel[anglei]->GetEntries()<<std::endl;
//  FitSpectrum(datasel[anglei], Ebeam, nameangle);
//}


  // finish program here
  return 0;
 
}

#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include "TPaveText.h"
#include "TMath.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooCBShape.h"
#include "RooGenericPdf.h"
#include "RooChebychev.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooMsgService.h"
#include "TGaxis.h"
using namespace RooFit;


double exp1(double *x, double *par)
{
    return TMath::Exp(par[0]*(x[0]-par[1]))+par[2];
}


double FitSpectrum(TTree *&dataraw, double beame, const char* namesfx)
{
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   std::cout<<"Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;
   double beamlow=peakvalue-0.2;
   double beamup=peakvalue+0.2;
   // try to use roofit
   RooRealVar x("x1","momentum",peakvalue,beamlow,beamup,"GeV");
   //RooRealVar x2("x2","momentum",peakvalue,beamlow,beamup,"GeV");
   //RooRealVar x("x1","momentum",peakvalue-0.1,beame);
   //RooRealVar x("x1","momentum",peakvalue-0.1,peakvalue+0.05);
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.003,peakvalue+0.003);
   RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.013);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   //RooRealVar mean2("mean2","mean of gaussian",beame,peakvalue,beame+0.1);
   RooRealVar mean2("mean2","mean of gaussian",peakvalue+0.2,peakvalue+0.05,3.0);
   RooRealVar sigma2("sigma2","width of gaussian",0.09,0.03,0.15);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean2,sigma2);
   
   RooRealVar co1("co1","coefficient #1",   0 ,-100.,100.);
   //RooRealVar co2("co2","coefficient #1",   0 ,-100.,100.);
 //RooRealVar co3("co3","coefficient #1",0.1,-100.,100.);
   RooChebychev ground("ground","background",x,RooArgList(co1));
   RooRealVar signal("signal"," ",1000,0,10000000);//event number
   RooRealVar background("background"," ",20000,0,100000000);
 //RooRealVar a0("a0","coefficient #0",100,-100000,100000);
 //RooRealVar a1("a1","coefficient #1",-10,-100000,100000);
 //RooRealVar a2("a2","coefficient #2",0,-100000,100000);
   //RooRealVar a3("a3","coefficient #2",0,-100000,100000);
 //RooPolynomial ground("ground","ground",x,RooArgList(a0,a1,a2));
     
   //RooRealVar sigma2("sigma2","width of gaussian",0.01,0.008,0.02);
   RooRealVar alpha1("alpha1","#alpha",1.7,5);
   RooRealVar nnn1("n1","n",100,1,200);
   RooCBShape cbshape("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);
 //RooRealVar sigma2("sigma2","width of gaussian",0.01,0.008,0.02);
 //RooRealVar alpha1("alpha1","#alpha",1.,-5,5);
 //RooRealVar nnn1("n1","n",100,1,200);
 //RooCBShape cbshape1("cbshape1","crystal ball",x,mean2,sigma2,alpha1,nnn1);

   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   //RooDataHist *data_6pi;
   
   // fit bhabha sigma and mean
   RooRealVar meanb("meanb","mean of gaussian",beame,beame-0.015,beame+0.005);
   RooRealVar sigmab("sigmab","width of gaussian",0.005,0.004,0.02);
 //  RooGaussian gausb("gausb","gauss(x,m,s)",x,meanb,sigmab);
   
   RooRealVar alphab("alphab","#alpha",1.0,5);
   RooRealVar nnnb("nnnb","n",100,1,200);
   RooCBShape cbshapeb("cbshapeb","crystal ball",x,meanb,sigmab,alphab,nnnb);
 
   RooRealVar meane("meane","mean of gaussian",beame,beame-0.005,beame+0.003);
   RooRealVar sigmae("sigmae","width of gaussian",0.005,0.004,0.02);
   RooGaussian gause("gause","gauss(x,m,s)",x,meane,sigmae);
   
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1 = new TCanvas("","",800,600);

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title("fit p"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   //sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2),RooArgList(signal,background));
   //sum = new RooAddPdf("sum","sum",RooArgList(cbshape,gaus2),RooArgList(signal,background));
   //sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,gaus3),RooArgList(signal,background,signal1));
   sum = new RooAddPdf("sum","sum",RooArgList(cbshape,cbshapeb),RooArgList(signal,background));
   Npar = 10;
   //sigma.setVal(0.035);
   //signal.setVal(1200);
   //background.setVal(200);
   //co1.setVal(0);
   //sum->fitTo(*dataset,Range(peakvalue-1.0*(beame-peakvalue),peakvalue+0.8*(beame-peakvalue)));
   //sum->fitTo(*dataset,Range(peakvalue-0.07,peakvalue+0.07));
   //x.setRange("sigragi",peakvalue-0.1,peakvalue+0.07);
   //x.setRange("sigragi",peakvalue-0.1,beame+0.01);
   x.setRange("sigragi",peakvalue-0.1,beame+0.012);
   //sum->fitTo(*dataset);
   sum->fitTo(*dataset,Range("sigragi"));
   //sum->fitTo(*dataset,Range(peakvalue-0.1,beame+0.01));
   //sum->fitTo(*dataset,"e",Range(peakvalue-0.1,beame));
   //gause.fitTo(*dataset,Range(beame-0.02, beame+0.01));
   dataset->plotOn(xframe);
   //sum->plotOn(xframe,Components(cbshape1),LineStyle(2),LineColor(2));
   //sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3));
   //sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(cbshape),LineStyle(2),LineColor(2)   ,Range(peakvalue-0.1, peakvalue+0.05)  );
   sum->plotOn(xframe,Components(cbshapeb),LineStyle(2),LineColor(3)  ,Range(peakvalue-0.1, peakvalue+0.05)  );
   //sum->plotOn(xframe,Components(cbshape),LineStyle(2),LineColor(2) );
   //sum->plotOn(xframe,Components(cbshapeb),LineStyle(2),LineColor(3));
   //sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(3));
   //gause.plotOn(xframe,LineStyle(2),LineColor(3));
   sum->plotOn(xframe  ,Range(peakvalue-0.1, peakvalue+0.065)  );
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt = new TPaveText(0.15,0.65,0.45,0.90,"BRNDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(4000);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.035);
   sprintf(tmpchr,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
 //sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
 //pt->AddText(tmpchr);
   sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
   sprintf(tmpchr,"p_spectrum_%s",namesfx);
   c1->SetName(tmpchr);
   c1->Write();

 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   ofstream fitpar("fitpar",std::ios::app);
   fitpar<<" ene = "<< beame*2 <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
   fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
   fitpar<<"\t sigNo = " << signal.getVal() << "\t sigNoE = " << signal.getError();
   fitpar<<"\t bckNo = " << background.getVal() << "\t bckNoE = " << background.getError();
   fitpar<< "\t chi2 = " << xframe->chiSquare(Npar);
   fitpar<< std::endl;;
   //fitpar<<"\t e mean = " << meane.getVal() << "\t e sigma = " << sigmae.getVal()<<std::endl;
   
   double intel3 = mean.getVal()-3*sigma.getVal();
   double inteu3 = mean.getVal()+3*sigma.getVal();
   double intel5 = mean.getVal()-5*sigma.getVal();
   double inteu5 = mean.getVal()+5*sigma.getVal();
   //RooRealVar xint("xint","xint",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame+0.01);
   x.setRange("sigrag3",intel3, inteu3);
   x.setRange("sigrag5",intel5, inteu5);
   RooAbsReal* intsigi = cbshape.createIntegral( x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intbcki = cbshapeb.createIntegral(x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intsig3 = cbshape.createIntegral( x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intbck3 = cbshapeb.createIntegral(x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intsig5 = cbshape.createIntegral( x,  NormSet(x),  Range("sigrag5"));
   RooAbsReal* intbck5 = cbshapeb.createIntegral(x,  NormSet(x),  Range("sigrag5"));
   double signalN3 =   signal.getVal()*intsig3->getVal();
   double backN3 = background.getVal()*intbck3->getVal();
   double signalN5 =   signal.getVal()*intsig5->getVal();
   double backN5 = background.getVal()*intbck5->getVal();
 //double signalN3nom =   signal.getVal()*intsig3->getVal()/intsigi->getVal();
 //double backN3nom = background.getVal()*intbck3->getVal()/intbcki->getVal();
 //double signalN5nom =   signal.getVal()*intsig5->getVal()/intsigi->getVal();
 //double backN5nom = background.getVal()*intbck5->getVal()/intbcki->getVal();
   ofstream angSNR("angSNR.dat",std::ios::app);
   angSNR << namesfx <<"\t"<<signalN3 << "\t" << backN3 <<"\t"<<signalN5 << "\t" << backN5 << std::endl;
 //epSNR << "\t"<<signalN3nom << "\t" << backN3nom <<"\t"<<signalN5nom << "\t" << backN5nom<< std::endl;
   std::cout<< "Total signal int is "<< intsigi->getVal() <<" Total bck int is "<< intbcki->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig3->getVal() <<" Total bck int is "<< intbck3->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig5->getVal() <<" Total bck int is "<< intbck5->getVal()<<std::endl;

   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
}


double FitSpectrum1(TTree *&dataraw, double beame, const char* namesfx)
{
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   double beamlow=peakvalue-0.2;
   double beamup=peakvalue+0.2;
   // try to use roofit
   RooRealVar x("x1","momentum",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,beamlow,beamup);
   RooRealVar sigma("sigma","width of gaussian",0.013,0.004,0.018);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   //RooRealVar mean2("mean2","mean of gaussian",beame,peakvalue,beame+0.1);
   RooRealVar mean2("mean2","mean of gaussian",peakvalue+0.2,peakvalue+0.05,3.0);
   RooRealVar sigma2("sigma2","width of gaussian",0.09,0.03,0.15);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean2,sigma2);
   
   RooRealVar co1("co1","coefficient #1",   0 ,-100.,100.);
   RooRealVar co2("co2","coefficient #1",   0 ,-100.,100.);
   RooRealVar co3("co3","coefficient #1",  0.1,-100.,100.);
   RooChebychev ground("ground","background",x,RooArgList(co1,co2,co3));
   
   RooRealVar signal("signal"," ",1000,0,1000000000);//event number
   RooRealVar background("background"," ",20000,0,100000000);
 //RooRealVar a0("a0","coefficient #0",100,-100000,100000);
 //RooRealVar a1("a1","coefficient #1",-10,-100000,100000);
 //RooRealVar a2("a2","coefficient #2",0,-100000,100000);
   //RooRealVar a3("a3","coefficient #2",0,-100000,100000);
 //RooPolynomial ground("ground","ground",x,RooArgList(a0,a1,a2));
   
   // signal
   RooRealVar alpha1("alpha1","#alpha",1.,-5,5);
   RooRealVar nnn1("n1","n",100,1,200);
   RooCBShape cbshape("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);

 //RooRealVar sigma2("sigma2","width of gaussian",0.01,0.008,0.02);
 //RooRealVar alpha1("alpha1","#alpha",1.,-5,5);
 //RooRealVar nnn1("n1","n",100,1,200);
 //RooCBShape cbshape1("cbshape1","crystal ball",x,mean2,sigma2,alpha1,nnn1);

   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   //RooDataHist *data_6pi;
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1 = new TCanvas("","",800,600);

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title("fit p"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   //sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2),RooArgList(signal,background));
   //sum = new RooAddPdf("sum","sum",RooArgList(cbshape,gaus2),RooArgList(signal,background));
   sum = new RooAddPdf("sum","sum",RooArgList(cbshape,ground),RooArgList(signal,background));
   //sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,gaus3),RooArgList(signal,background,signal1));
   //sum = new RooAddPdf("sum","sum",RooArgList(cbshape1,ground),RooArgList(signal,background));
   Npar = 6;
   //sigma.setVal(0.035);
   //signal.setVal(1200);
   //background.setVal(200);
   //co1.setVal(0);
   //sum->fitTo(*dataset,Range(peakvalue-1.0*(beame-peakvalue),peakvalue+0.8*(beame-peakvalue)));
   sum->fitTo(*dataset,Range(peakvalue-0.08,peakvalue+0.06));
   dataset->plotOn(xframe);
   //sum->plotOn(xframe,Components(cbshape1),LineStyle(2),LineColor(2));
   //sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3));
   //sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(cbshape),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt = new TPaveText(0.60,0.5,0.90,0.90,"BRNDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(4000);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.035);
   sprintf(tmpchr,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
 //sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
 //pt->AddText(tmpchr);
   //sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
   sprintf(tmpchr,"#chi^{2} = %5.6f",xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
   sprintf(tmpchr,"p_spectrum_%s",namesfx);
   c1->SetName(tmpchr);
   c1->Write();

 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   
   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
}

double FitSpectrum2(TTree *&dataraw, double beame, const char* namesfx)
{
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double mmu = 0.1057;
   double mpi = 0.13957;
   double peakvalue = beame;
   double mispeak = 2*sqrt(pow(beame/2,2)-pow(mpi,2)+pow(mka,2));
   double beamlow=peakvalue-0.2;
   double beamup=peakvalue+0.2;

   double sigmamean = 0.005+0.01*(beame-2.0);
   double sigmalow  = sigmamean-0.001;
   // try to use roofit
   RooRealVar x("mass","M(KK)",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.002,peakvalue+0.002);
   RooRealVar sigma("sigma","width of gaussian",sigmamean,sigmalow,0.018);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooRealVar mean2("mean2","mean of gaussian",mispeak+0.1,peakvalue+0.1,mispeak+0.8);
   RooRealVar sigma2("sigma2","width of gaussian",0.09,0.01,0.15);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean2,sigma2);
   
   //double co2up = 0.5+(beame-2.0);
   RooRealVar co1("co1","coefficient #1", 0,-100.,100.);
   RooRealVar co2("co2","coefficient #1", 0,-10,10);
   //RooRealVar co3("co3","coefficient #1",0.1,-100.,100.);
   RooChebychev ground("ground","background",x,RooArgList(co1,co2));
   RooRealVar signal("signal"," ",1000,0,1000000000);//event number
   RooRealVar background("background"," ",20000,0,100000000);
 //RooRealVar a0("a0","coefficient #0",100,-100000,100000);
 //RooRealVar a1("a1","coefficient #1",-10,-100000,100000);
 //RooRealVar a2("a2","coefficient #2",0,-100000,100000);
   //RooRealVar a3("a3","coefficient #2",0,-100000,100000);
 //RooPolynomial ground("ground","ground",x,RooArgList(a0,a1,a2));
     
   //RooRealVar sigma2("sigma2","width of gaussian",0.01,0.008,0.05);
   RooRealVar alpha1("alpha1","#alpha",1.2,0.5,5);
   RooRealVar nnn1("n1","n",100,1,200);
   RooCBShape cbshape1("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);
   //RooCBShape cbshape1("cbshape1","crystal ball",x,mean2,sigma2,alpha1,nnn1);

    
    //TF1 *exp2 = new TF1("exp2",exp1,0,4,3);
   // define user function
     double amean = 0;
     RooRealVar a("a","a",0, -1,1000);
     RooRealVar b("b","b",beamlow, 0,5);
     RooRealVar c("c","c",0, -10000,10000);
   //RooAbsPdf myexp = *(RooFit::bindPdf("exp",exp2,x,a,b,c));
     RooGenericPdf myexp("exp","exp","exp(a*(mass-b))+c",RooArgSet(x,a,b,c));

   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   //RooDataHist *data_6pi;
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1 = new TCanvas("","",800,600);

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title("fit mass"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   //sum = new RooAddPdf("sum","sum",RooArgList(cbshape1,myexp),RooArgList(signal,background));
   //sum = new RooAddPdf("sum","sum",RooArgList(cbshape1,ground),RooArgList(signal,background));
   sum = new RooAddPdf("sum","sum",RooArgList(cbshape1,gaus2),RooArgList(signal,background));
   Npar = 8;
   //sigma.setVal(0.035);
   //signal.setVal(1200);
   //background.setVal(200);
   //co1.setVal(0);
   sum->fitTo(*dataset,Range(peakvalue-0.6*(mispeak-peakvalue),peakvalue+0.6*(mispeak-peakvalue)));
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(cbshape1),LineStyle(2),LineColor(2));
   //sum->plotOn(xframe,Components(myexp),LineStyle(2),LineColor(3));
   //sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3));
   sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
   //sum->plotOn(xframe,Components(cbshape1),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt = new TPaveText(0.60,0.5,0.90,0.90,"BRNDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(4000);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.035);
   sprintf(tmpchr,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
 //sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
 //pt->AddText(tmpchr);
   sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
   sprintf(tmpchr,"mass_spectrum_%s",namesfx);
   c1->SetName(tmpchr);
   c1->Write();
   sprintf(tmpchr,"mass_spectrum_%s.pdf",namesfx);
   c1->Print(tmpchr);

 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   
   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
}

/*
double FitSpectrum3(TTree *&dataraw, double beame, const char* namesfx)
{
   TF1 fsigma("fsigma","-0.00705519+0.00600716*x",2.0,3.1);
   double sigma_ini = fsigma.Eval(beame);
   TCanvas *c1 = new TCanvas("","",800,600);
   
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double mmu = 0.1057;
   double Ecm = beame;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   double pmu = sqrt(pow(beame,2)-pow(mmu,2));
   std::cout<<"Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;
   double beamlow=peakvalue-0.1;
   double beamup=peakvalue+0.2;
   if (Ecm<2.6) {
     beamlow=peakvalue-0.05;
     beamup=peakvalue+0.15;
   }
   // save histogram of pmu
  //TH1D *hp = new TH1D("hp","hp",100,beamlow,beamup);
  //dataraw->Draw("x1>>hp");
  //char name[100];
  //sprintf(name,"hp_mcKK_%.4f",beame*2);
  //hp->Write(name);
  //delete hp;
  //return 0;
   //
   //
   // try to use roofit
   //
   RooRealVar x("x1","momentum",peakvalue,beamlow,beamup,"GeV");
   // signal
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.003,peakvalue+0.003);
   RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.013);
   //RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.014);
   RooRealVar alpha1("alpha1","#alpha",1.5,1.0,5.0);
   RooRealVar nnn1("n1","n",2,1,10);
   RooCBShape cbshape("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);
   RooRealVar sigma2("sigma2","width of gaussian",0.007,0.005,0.02);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma2);
   RooRealVar frac1("frac1","frac1",0.9,0.6,1.0);

   RooAddPdf sig("sig","signal",RooArgList(cbshape,gaus),RooArgList(frac1));
   
   // background, mainly di-mu
   RooRealVar meanb("meanb","mean of gaussian",pmu,pmu-0.002,pmu+0.003);
   //RooRealVar sigmab("sigmab","width of gaussian",sigma_ini-0.001, sigma_ini-0.001,sigma_ini+0.001);
   RooRealVar sigmab("sigmab","width of gaussian",sigma_ini, sigma_ini-0.001,sigma_ini+0.001);
   //RooRealVar sigmab("sigmab","width of gaussian",0.005,0.002,0.02);
 //  RooGaussian gausb("gausb","gauss(x,m,s)",x,meanb,sigmab);
   //RooRealVar alphab("alphab","#alpha",0.8,0.5,1.5);
   RooRealVar alphab("alphab","#alpha",1.2,1.0,1.5);
   //RooRealVar alphab("alphab","#alpha",0.5,0.1,1.5);
   //RooRealVar nnnb("nnnb","n",0.5,0.4,2.0);
   RooRealVar nnnb("nnnb","n",0.5,0.1,2.0);
   //RooRealVar nnnb("nnnb","n",beame);
   RooCBShape cbshapeb("cbshapeb","crystal ball",x,meanb,sigmab,alphab,nnnb);
   //RooRealVar meane("meane","mean of gaussian",beame,beame-0.005,beame+0.003);
   RooRealVar sigmae("sigmae","width of gaussian",0.01,0.005,0.08);
   RooGaussian gause("gause","gauss(x,m,s)",x,meanb,sigmae);
   RooRealVar frac2("frac2","frac2",0.9,0.6,1.0);
   RooAddPdf bck("bck","signal",RooArgList(cbshapeb,gause),RooArgList(frac2));
   
   RooRealVar signal("signal"," ",1000,0,1000000000);//event number
   RooRealVar background("background"," ",20000,0,100000000);
 
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title(""));
   c1->SetMargin(0.13,0.1,0.13,0.1);
   TGaxis::SetMaxDigits(3);
   xframe->SetTitle("");
   xframe->GetXaxis()->SetTitle("p (GeV/c)");
   xframe->GetXaxis()->SetNdivisions(505);
   xframe->GetXaxis()->SetTitleSize(0.05);
   xframe->GetXaxis()->SetLabelSize(0.05);
   xframe->GetXaxis()->SetTitleOffset(1.1);
   xframe->GetYaxis()->SetNdivisions(505);
   xframe->GetYaxis()->SetTitleSize(0.05);
   xframe->GetYaxis()->SetTitleOffset(1.25);
   xframe->GetYaxis()->SetLabelSize(0.05);
   xframe->GetYaxis()->SetNoExponent(0);
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(sig,bck),RooArgList(signal,background));
   Npar = 10;
   //x.setRange("sigragi",peakvalue-0.1,beame+0.012);
   //sum->fitTo(*dataset);
   sum->fitTo(*dataset);
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bck),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt;
   if (Ecm<2.6) pt = new TPaveText(0.5,0.55,0.85,0.90,"BRNDC");
   else pt = new TPaveText(0.15,0.55,0.50,0.90,"BRNDC");
   pt->SetBorderSize(0);
   //pt->SetFillStyle(4000);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.045);
   sprintf(tmpchr,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"signal = %.1f #pm %.1f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
   sprintf(tmpchr,"p_spectrum_%s",namesfx);
   c1->SetName(tmpchr);
   c1->Write();
   sprintf(tmpchr,"p_spectrum_%s.pdf",namesfx);
   c1->Print(tmpchr);

 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   ofstream fitpar("fitpar",std::ios::app);
 //fitpar<<" ene = "<< beame*2 <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
 //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
 //fitpar<<"\t sigNo = " << signal.getVal() << "\t sigNoE = " << signal.getError();
 //fitpar<<"\t bckNo = " << background.getVal() << "\t bckNoE = " << background.getError();
 //fitpar<< "\t chi2 = " << xframe->chiSquare(Npar);
 //fitpar<< std::endl;
   fitpar<< beame*2 <<"\t ";
   
   fitpar<< mean.getVal() <<"\t ";
   fitpar<< sigma.getVal() <<"\t ";
   fitpar<< sigma2.getVal() <<"\t ";
   fitpar<< alpha1.getVal()<<"\t ";
   fitpar<< nnn1.getVal()<<"\t ";
   fitpar<< frac1.getVal()<<"\t ";
   
   fitpar<< meanb.getVal()<<"\t ";
   fitpar<< sigmab.getVal()<<"\t ";
   fitpar<< sigmae.getVal()<<"\t ";
   fitpar<< alphab.getVal()<<"\t ";
   fitpar<< nnnb.getVal()<<"\t ";
   fitpar<< frac2.getVal()<<"\t ";
   
   fitpar<< endl;
   //fitpar<<"\t e mean = " << meane.getVal() << "\t e sigma = " << sigmae.getVal()<<std::endl;
   
   double intel3 = mean.getVal()-3*sigma.getVal();
   double inteu3 = mean.getVal()+3*sigma.getVal();
   double intel5 = mean.getVal()-5*sigma.getVal();
   double inteu5 = mean.getVal()+5*sigma.getVal();
   //RooRealVar xint("xint","xint",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame+0.01);
   x.setRange("sigrag3",intel3, inteu3);
   x.setRange("sigrag5",intel5, inteu5);
   RooAbsReal* intsigi = sig.createIntegral( x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intbcki = bck.createIntegral(x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intsig3 = sig.createIntegral( x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intbck3 = bck.createIntegral(x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intsig5 = sig.createIntegral( x,  NormSet(x),  Range("sigrag5"));
   RooAbsReal* intbck5 = bck.createIntegral(x,  NormSet(x),  Range("sigrag5"));
   double signalN3 =   signal.getVal()*intsig3->getVal();
   double backN3 = background.getVal()*intbck3->getVal();
   double signalN5 =   signal.getVal()*intsig5->getVal();
   double backN5 = background.getVal()*intbck5->getVal();
 //double signalN3nom =   signal.getVal()*intsig3->getVal()/intsigi->getVal();
 //double backN3nom = background.getVal()*intbck3->getVal()/intbcki->getVal();
 //double signalN5nom =   signal.getVal()*intsig5->getVal()/intsigi->getVal();
 //double backN5nom = background.getVal()*intbck5->getVal()/intbcki->getVal();
   ofstream angSNR("angSNR.dat",std::ios::app);
   angSNR << namesfx <<"\t"<<signalN3 << "\t" << backN3 <<"\t"<<signalN5 << "\t" << backN5 << std::endl;
 //epSNR << "\t"<<signalN3nom << "\t" << backN3nom <<"\t"<<signalN5nom << "\t" << backN5nom<< std::endl;
   std::cout<< "Total signal int is "<< intsigi->getVal() <<" Total bck int is "<< intbcki->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig3->getVal() <<" Total bck int is "<< intbck3->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig5->getVal() <<" Total bck int is "<< intbck5->getVal()<<std::endl;

   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
}
*/
/*
double FitSpectrum3(TTree *&dataraw, double beame, const char* namesfx)
{
   TF1 fsigma("fsigma","-0.00705519+0.00600716*x",2.0,3.1);
   double sigma_ini = fsigma.Eval(beame);
   TCanvas *c1 = new TCanvas("","",800,600);
   
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double mmu = 0.1057;
   double Ecm = beame;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   double pmu = sqrt(pow(beame,2)-pow(mmu,2));
   std::cout<<"Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;
   double beamlow=peakvalue-0.1;
   double beamup=peakvalue+0.2;
   if (Ecm<2.6) {
     beamlow=peakvalue-0.05;
     beamup=peakvalue+0.15;
   }
   // save histogram of pmu
  //TH1D *hp = new TH1D("hp","hp",100,beamlow,beamup);
  //dataraw->Draw("x1>>hp");
  //char name[100];
  //sprintf(name,"hp_mcKK_%.4f",beame*2);
  //hp->Write(name);
  //delete hp;
  //return 0;
   //
   //
   // try to use roofit
   //
   RooRealVar x("x1","momentum",peakvalue,beamlow,beamup,"GeV/c");
   // signal
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.003,peakvalue+0.003);
   RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.013);
   //RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.014);
   RooRealVar alpha1("alpha1","#alpha",1.5,1.0,5.0);
   RooRealVar nnn1("n1","n",2,1,10);
   RooCBShape cbshape("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);
   RooRealVar sigma2("sigma2","width of gaussian",0.007,0.005,0.02);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma2);
   RooRealVar frac1("frac1","frac1",0.9,0.1,1.0);

   RooAddPdf sig("sig","signal",RooArgList(cbshape,gaus),RooArgList(frac1));
   
   // background, mainly di-mu
   RooRealVar meanb("meanb","mean of gaussian",pmu,pmu-0.002,pmu+0.003);
   //RooRealVar sigmab("sigmab","width of gaussian",sigma_ini-0.001, sigma_ini-0.001,sigma_ini+0.001);
   RooRealVar sigmab("sigmab","width of gaussian",sigma_ini, sigma_ini-0.001,sigma_ini+0.001);
   //RooRealVar sigmab("sigmab","width of gaussian",0.005,0.002,0.02);
 //  RooGaussian gausb("gausb","gauss(x,m,s)",x,meanb,sigmab);
   RooRealVar alphab("alphab","#alpha",0.5,0.5,5);
   //RooRealVar alphab("alphab","#alpha",1.2,1.0,1.5);
   //RooRealVar alphab("alphab","#alpha",0.5,0.1,1.5);
   //RooRealVar nnnb("nnnb","n",0.5,0.4,2.0);
   RooRealVar nnnb("nnnb","n",0.6,0.1,1.5);
   //RooRealVar nnnb("nnnb","n",beame);
   RooCBShape cbshapeb("cbshapeb","crystal ball",x,meanb,sigmab,alphab,nnnb);
   //RooRealVar meane("meane","mean of gaussian",beame,beame-0.005,beame+0.003);
   RooRealVar sigmae("sigmae","width of gaussian",0.01,0.005,0.08);
   RooGaussian gause("gause","gauss(x,m,s)",x,meanb,sigmae);
   RooRealVar frac2("frac2","frac2",0.9,0.6,1.0);
   RooAddPdf bck("bck","signal",RooArgList(cbshapeb,gause),RooArgList(frac2));
   
   RooRealVar signal("signal"," ",1000,0,1000000000);//event number
   RooRealVar background("background"," ",20000,0,100000000);
 
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title(""));
   c1->SetMargin(0.13,0.1,0.13,0.1);
   TGaxis::SetMaxDigits(3);
   xframe->SetTitle("");
   xframe->GetXaxis()->SetTitle("p (GeV/c)");
   xframe->GetXaxis()->SetNdivisions(505);
   xframe->GetXaxis()->SetTitleSize(0.05);
   xframe->GetXaxis()->SetLabelSize(0.05);
   xframe->GetXaxis()->SetTitleOffset(1.1);
   xframe->GetYaxis()->SetNdivisions(505);
   xframe->GetYaxis()->SetTitleSize(0.05);
   xframe->GetYaxis()->SetTitleOffset(1.25);
   xframe->GetYaxis()->SetLabelSize(0.05);
   xframe->GetYaxis()->SetNoExponent(0);
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(sig,bck),RooArgList(signal,background));
   Npar = 10;
   //x.setRange("sigragi",peakvalue-0.1,beame+0.012);
   //sum->fitTo(*dataset);
   sum->fitTo(*dataset);
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bck),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt;
   if (Ecm<2.65) pt = new TPaveText(0.5,0.55,0.85,0.90,"BRNDC");
   else pt = new TPaveText(0.15,0.55,0.50,0.90,"BRNDC");
   pt->SetBorderSize(0);
   //pt->SetFillStyle(4000);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.045);
   sprintf(tmpchr,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"signal = %.1f #pm %.1f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
   //sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,xframe->chiSquare(Npar));
   sprintf(tmpchr,"#chi^{2} = %5.3f",xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
   sprintf(tmpchr,"p_spectrum_%s",namesfx);
   c1->SetName(tmpchr);
   c1->Write();
   sprintf(tmpchr,"p_spectrum_%s.pdf",namesfx);
   c1->Print(tmpchr);

 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   ofstream fitpar("fitpar",std::ios::app);
 //fitpar<<" ene = "<< beame*2 <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
 //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
 //fitpar<<"\t sigNo = " << signal.getVal() << "\t sigNoE = " << signal.getError();
 //fitpar<<"\t bckNo = " << background.getVal() << "\t bckNoE = " << background.getError();
 //fitpar<< "\t chi2 = " << xframe->chiSquare(Npar);
 //fitpar<< std::endl;
   fitpar<< beame*2 <<"\t ";
   
   fitpar<< mean.getVal() <<"\t ";
   fitpar<< sigma.getVal() <<"\t ";
   fitpar<< sigma2.getVal() <<"\t ";
   fitpar<< alpha1.getVal()<<"\t ";
   fitpar<< nnn1.getVal()<<"\t ";
   fitpar<< frac1.getVal()<<"\t ";
   
   fitpar<< meanb.getVal()<<"\t ";
   fitpar<< sigmab.getVal()<<"\t ";
   fitpar<< sigmae.getVal()<<"\t ";
   fitpar<< alphab.getVal()<<"\t ";
   fitpar<< nnnb.getVal()<<"\t ";
   fitpar<< frac2.getVal()<<"\t ";
   
   fitpar<< endl;
   //fitpar<<"\t e mean = " << meane.getVal() << "\t e sigma = " << sigmae.getVal()<<std::endl;
   
   double intel3 = mean.getVal()-3*sigma.getVal();
   double inteu3 = mean.getVal()+3*sigma.getVal();
   double intel5 = mean.getVal()-5*sigma.getVal();
   double inteu5 = mean.getVal()+5*sigma.getVal();
   //RooRealVar xint("xint","xint",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame+0.01);
   x.setRange("sigrag3",intel3, inteu3);
   x.setRange("sigrag5",intel5, inteu5);
   RooAbsReal* intsigi = sig.createIntegral( x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intbcki = bck.createIntegral(x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intsig3 = sig.createIntegral( x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intbck3 = bck.createIntegral(x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intsig5 = sig.createIntegral( x,  NormSet(x),  Range("sigrag5"));
   RooAbsReal* intbck5 = bck.createIntegral(x,  NormSet(x),  Range("sigrag5"));
   double signalN3 =   signal.getVal()*intsig3->getVal();
   double backN3 = background.getVal()*intbck3->getVal();
   double signalN5 =   signal.getVal()*intsig5->getVal();
   double backN5 = background.getVal()*intbck5->getVal();
 //double signalN3nom =   signal.getVal()*intsig3->getVal()/intsigi->getVal();
 //double backN3nom = background.getVal()*intbck3->getVal()/intbcki->getVal();
 //double signalN5nom =   signal.getVal()*intsig5->getVal()/intsigi->getVal();
 //double backN5nom = background.getVal()*intbck5->getVal()/intbcki->getVal();
   ofstream angSNR("angSNR.dat",std::ios::app);
   angSNR << namesfx <<"\t"<<signalN3 << "\t" << backN3 <<"\t"<<signalN5 << "\t" << backN5 << std::endl;
 //epSNR << "\t"<<signalN3nom << "\t" << backN3nom <<"\t"<<signalN5nom << "\t" << backN5nom<< std::endl;
   std::cout<< "Total signal int is "<< intsigi->getVal() <<" Total bck int is "<< intbcki->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig3->getVal() <<" Total bck int is "<< intbck3->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig5->getVal() <<" Total bck int is "<< intbck5->getVal()<<std::endl;

   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
}
*/

double FitSpectrum3(TTree *&dataraw, double beame, const char* namesfx)
{
   TF1 fsigma("fsigma","-0.00705519+0.00600716*x",2.0,3.1);
   double sigma_ini = fsigma.Eval(beame);
   TCanvas *c1 = new TCanvas("","",800,600);
   
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double mmu = 0.1057;
   double Ecm = beame;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   double pmu = sqrt(pow(beame,2)-pow(mmu,2));
   std::cout<<"Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;
   double beamlow=peakvalue-0.1;
   double beamup=peakvalue+0.2;
   if (Ecm<2.6) {
     beamlow=peakvalue-0.05;
     beamup=peakvalue+0.15;
   }
   // save histogram of pmu or pK
  //TH1D *hp = new TH1D("hp","hp",100,beamlow,beamup);
  //dataraw->Draw("x1>>hp");
  //char name[100];
  //sprintf(name,"hp_mckk_%.4f",beame*2);
  //hp->Write(name);
  //delete hp;
  //return 0;
   //
   //
   // try to use roofit
   //
   RooRealVar x("x1","momentum",peakvalue,beamlow,beamup,"GeV/c");
   // signal
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.003,peakvalue+0.003);
   RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.013);
   //RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.014);
   RooRealVar alpha1("alpha1","#alpha",2,1.0,5.0);
   RooRealVar nnn1("n1","n",2,1,10);
   RooCBShape cbshape("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);
   RooRealVar sigma2("sigma2","width of gaussian",0.007,0.005,0.02);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma2);
   RooRealVar frac1("frac1","frac1",0.8,0.1,1.0);

   RooAddPdf sig("sig","signal",RooArgList(cbshape,gaus),RooArgList(frac1));
   
   // background, mainly di-mu
   RooRealVar meanb("meanb","mean of gaussian",pmu,pmu-0.002,pmu+0.003);
   //RooRealVar sigmab("sigmab","width of gaussian",sigma_ini-0.001, sigma_ini-0.001,sigma_ini+0.001);
   RooRealVar sigmab("sigmab","width of gaussian",sigma_ini, sigma_ini-0.002,sigma_ini+0.002);
   //RooRealVar sigmab("sigmab","width of gaussian",0.005,0.002,0.02);
 //  RooGaussian gausb("gausb","gauss(x,m,s)",x,meanb,sigmab);
   RooRealVar alphab("alphab","#alpha",1.5,0.5,5);
   //RooRealVar alphab("alphab","#alpha",1.2,1.0,1.5);
   //RooRealVar alphab("alphab","#alpha",0.5,0.1,1.5);
   //RooRealVar nnnb("nnnb","n",0.5,0.4,2.0);
   RooRealVar nnnb("nnnb","n",0.6,0.1,1.5);
   //RooRealVar nnnb("nnnb","n",beame);
   RooCBShape cbshapeb("cbshapeb","crystal ball",x,meanb,sigmab,alphab,nnnb);
   //RooRealVar meane("meane","mean of gaussian",beame,beame-0.005,beame+0.003);
   RooRealVar sigmae("sigmae","width of gaussian",0.01,0.005,0.08);
   RooGaussian gause("gause","gauss(x,m,s)",x,meanb,sigmae);
   RooRealVar frac2("frac2","frac2",0.9,0.6,1.0);
   RooAddPdf bck("bck","signal",RooArgList(cbshapeb,gause),RooArgList(frac2));
   
   RooRealVar signal("signal"," ",1000,0,1000000);//event number
   RooRealVar background("background"," ",2000,0,1000000);
 
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title(""));
   c1->SetMargin(0.13,0.1,0.13,0.1);
   TGaxis::SetMaxDigits(3);
   xframe->SetTitle("");
   xframe->GetXaxis()->SetTitle("p (GeV/c)");
   xframe->GetXaxis()->SetNdivisions(505);
   xframe->GetXaxis()->SetTitleSize(0.05);
   xframe->GetXaxis()->SetLabelSize(0.05);
   xframe->GetXaxis()->SetTitleOffset(1.1);
   xframe->GetYaxis()->SetNdivisions(505);
   xframe->GetYaxis()->SetTitleSize(0.05);
   xframe->GetYaxis()->SetTitleOffset(1.25);
   xframe->GetYaxis()->SetLabelSize(0.05);
   xframe->GetYaxis()->SetNoExponent(0);
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(sig,bck),RooArgList(signal,background));
   Npar = 10;
   //x.setRange("sigragi",peakvalue-0.1,beame+0.012);
   //sum->fitTo(*dataset);
   sum->fitTo(*dataset);
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bck),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt;
   if (Ecm<2.65) pt = new TPaveText(0.5,0.55,0.85,0.90,"BRNDC");
   else pt = new TPaveText(0.15,0.55,0.50,0.90,"BRNDC");
   pt->SetBorderSize(0);
   //pt->SetFillStyle(4000);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.045);
   sprintf(tmpchr,"#mu = %1.5f #pm %1.5f",mean.getVal(),mean.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma = %1.5f #pm %1.5f",sigma.getVal(),sigma.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"signal = %.1f #pm %.1f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
   //sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,xframe->chiSquare(Npar));
   sprintf(tmpchr,"#chi^{2} = %5.3f",xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
   sprintf(tmpchr,"p_spectrum_%s",namesfx);
   c1->SetName(tmpchr);
   c1->Write();
   sprintf(tmpchr,"p_spectrum_%s.pdf",namesfx);
   c1->Print(tmpchr);

 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   ofstream fitpar("fitpar",std::ios::app);
 //fitpar<<" ene = "<< beame*2 <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
 //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
 //fitpar<<"\t sigNo = " << signal.getVal() << "\t sigNoE = " << signal.getError();
 //fitpar<<"\t bckNo = " << background.getVal() << "\t bckNoE = " << background.getError();
 //fitpar<< "\t chi2 = " << xframe->chiSquare(Npar);
 //fitpar<< std::endl;
   fitpar<< beame*2 <<"\t ";
   
   fitpar<< mean.getVal() <<"\t ";
   fitpar<< sigma.getVal() <<"\t ";
   fitpar<< sigma2.getVal() <<"\t ";
   fitpar<< alpha1.getVal()<<"\t ";
   fitpar<< nnn1.getVal()<<"\t ";
   fitpar<< frac1.getVal()<<"\t ";
   
   fitpar<< meanb.getVal()<<"\t ";
   fitpar<< sigmab.getVal()<<"\t ";
   fitpar<< sigmae.getVal()<<"\t ";
   fitpar<< alphab.getVal()<<"\t ";
   fitpar<< nnnb.getVal()<<"\t ";
   fitpar<< frac2.getVal()<<"\t ";
   
   fitpar<< endl;
   //fitpar<<"\t e mean = " << meane.getVal() << "\t e sigma = " << sigmae.getVal()<<std::endl;
   
   double intel3 = mean.getVal()-3*sigma.getVal();
   double inteu3 = mean.getVal()+3*sigma.getVal();
   double intel5 = mean.getVal()-5*sigma.getVal();
   double inteu5 = mean.getVal()+5*sigma.getVal();
   //RooRealVar xint("xint","xint",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame+0.01);
   x.setRange("sigrag3",intel3, inteu3);
   x.setRange("sigrag5",intel5, inteu5);
   RooAbsReal* intsigi = sig.createIntegral( x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intbcki = bck.createIntegral(x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intsig3 = sig.createIntegral( x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intbck3 = bck.createIntegral(x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intsig5 = sig.createIntegral( x,  NormSet(x),  Range("sigrag5"));
   RooAbsReal* intbck5 = bck.createIntegral(x,  NormSet(x),  Range("sigrag5"));
   double signalN3 =   signal.getVal()*intsig3->getVal();
   double backN3 = background.getVal()*intbck3->getVal();
   double signalN5 =   signal.getVal()*intsig5->getVal();
   double backN5 = background.getVal()*intbck5->getVal();
 //double signalN3nom =   signal.getVal()*intsig3->getVal()/intsigi->getVal();
 //double backN3nom = background.getVal()*intbck3->getVal()/intbcki->getVal();
 //double signalN5nom =   signal.getVal()*intsig5->getVal()/intsigi->getVal();
 //double backN5nom = background.getVal()*intbck5->getVal()/intbcki->getVal();
   ofstream angSNR("angSNR.dat",std::ios::app);
   angSNR << namesfx <<"\t"<<signalN3 << "\t" << backN3 <<"\t"<<signalN5 << "\t" << backN5 << std::endl;
 //epSNR << "\t"<<signalN3nom << "\t" << backN3nom <<"\t"<<signalN5nom << "\t" << backN5nom<< std::endl;
   std::cout<< "Total signal int is "<< intsigi->getVal() <<" Total bck int is "<< intbcki->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig3->getVal() <<" Total bck int is "<< intbck3->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig5->getVal() <<" Total bck int is "<< intbck5->getVal()<<std::endl;

   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
}


#include "RooHistPdf.h"
#include "RooFFTConvPdf.h"
double FitSpectrum4(TTree *&dataraw, double beame, const char* namesfx, TFile* fileout)
{
   TCanvas *c1 = new TCanvas("","",800,600);
   
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double mmu = 0.1057;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   double pmu = sqrt(pow(beame,2)-pow(mmu,2));
   std::cout<<"Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;
   double beamlow=peakvalue-0.2;
   double beamup=peakvalue+0.2;
   //double beamlow=peakvalue-0.10;
   //double beamup=peakvalue+0.10;
   
   char hname[1000];
   TFile file("output/output_noF.root");
   sprintf(hname,"hp_mcKK_%.4f",beame*2);
   TH1D hk = *(TH1D*)(((TH1D*)file.Get(hname))->Clone());
   sprintf(hname,"hp_mumu_%.4f",beame*2);
   TH1D hmu = *(TH1D*)(((TH1D*)file.Get(hname))->Clone());
   file.Close();
   //if (hk==0 || hmu==0) {cout<<"Error : spare pointer!!!!"<< endl; exit(-1);}
   hk.Draw();
   c1->Print("hktmp.pdf");
   hmu.Draw();
   c1->Print("hmutmp.pdf");
   //
   fileout->cd();
   //
   // try to use roofit
   //
   RooRealVar x("x1","momentum",peakvalue,beamlow,beamup,"GeV");
    // signal
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.003,peakvalue+0.003);
   RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.014);
   RooRealVar alpha1("alpha1","#alpha",1.5,1.0,5.0);
   RooRealVar nnn1("n1","n",2,1,10);
   RooCBShape cbshape("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);
   RooRealVar sigma2("sigma2","width of gaussian",0.007,0.005,0.02);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma2);
   RooRealVar frac1("frac1","frac1",0.9,0.8,1.0);

   RooAddPdf sig("sig","signal",RooArgList(cbshape,gaus),RooArgList(frac1));
 
 
 
 //RooDataHist mckhist("mcKhist","mc K spec",x,&hk);
   RooDataHist mcmuhist("mcmuhist","mc mu spec",x,&hmu);
 //RooHistPdf kpdf("kpdf","signal",x,mckhist,2); 
   RooHistPdf mupdf("mupdf","background",x,mcmuhist,2);

 //RooRealVar mean("mean","mean of gaussian",0.0007,-0.001, 0.001);
 //RooRealVar sigma("sigma","width of gaussian",0.0024,0.002,0.003);
 //RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);

   // signal
 //RooFFTConvPdf sig("sig","signal",x,kpdf,gaus); 
   //RooHistPdf sig("sig","signal",x,mckhist,2); 
   // background, mainly di-mu
   //RooFFTConvPdf bck("bck","background",x,mupdf,gaus); 
   RooHistPdf bck("bck","background",x,mcmuhist,2);
   
   RooRealVar signal("signal"," ",2000,0,100000);//event number
   RooRealVar background("background"," ",200,0,100000);
 
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title("fit p"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(sig,bck),RooArgList(signal,background));
   //sum = new RooAddPdf("sum","sum",RooArgList(sig),RooArgList(signal));
   Npar = 4;
   //x.setRange("sigragi",peakvalue-0.1,beame+0.012);
   //sum->fitTo(*dataset,Range("sigragi"));
   //sum->fitTo(*dataset,Minimizer("Minuit","simplex"),Hesse(false));
   //sum->fitTo(*dataset,Minimizer("GSLSimAn"));
   //sum->fitTo(*dataset,Minimizer("GSLMultiMin"),Hesse(false) );
   sum->chi2FitTo(*dataset);
   dataset->plotOn(xframe,Binning(nBins));
   //sum->fitTo(mckhist);
   //mckhist.plotOn(xframe);
   sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bck),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt = new TPaveText(0.15,0.65,0.45,0.90,"BRNDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(4000);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.035);
// sprintf(tmpchr,"#mu = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
// pt->AddText(tmpchr);
// sprintf(tmpchr,"#sigma = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
// pt->AddText(tmpchr);
   sprintf(tmpchr,"N_{sig} = %.2f #pm %.2f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"N_{bck} = %.2f #pm %.2f",background.getVal(),background.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.3f",nBins,Npar,xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
   sprintf(tmpchr,"p_spectrum_%s",namesfx);
   c1->SetName(tmpchr);
   c1->Write();
   sprintf(tmpchr,"p_spectrum_%s.pdf",namesfx);
   c1->Print(tmpchr);
   
   
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
}

double m_par[12];

bool getPars(double energy)
{
 if (fabs(energy-2)<1e-4)	{m_par[0]= 0.870367	;m_par[1]=0.00460545	;m_par[2]= 0.0117588	;m_par[3]= 1.32605	;m_par[4]= 9.99997	;m_par[5]= 0.930997	;m_par[6]= 0.995603	;m_par[7]= 0.00595913	;m_par[8]= 0.0128092	;m_par[9]= 1.99593	;m_par[10]= 0.1	;m_par[11]= 0.635921	 ; return 1;}
 if (fabs(energy-2.05)<1e-4)	{m_par[0]= 0.898002	;m_par[1]=0.00485945	;m_par[2]= 0.0167493	;m_par[3]= 4.31001	;m_par[4]= 2.42602	;m_par[5]= 0.890877	;m_par[6]= 1.01999	;m_par[7]= 0.00625541	;m_par[8]= 0.0799418	;m_par[9]= 1.59312	;m_par[10]= 0.649966	;m_par[11]= 0.956378	 ; return 1;}
 if (fabs(energy-2.1)<1e-4)	{m_par[0]= 0.926588	;m_par[1]=0.00433223	;m_par[2]= 0.00866258	;m_par[3]= 1.98911	;m_par[4]= 1.7437	;m_par[5]= 0.719703	;m_par[6]= 1.04453	;m_par[7]= 0.00651913	;m_par[8]= 0.037212	;m_par[9]= 1.80297	;m_par[10]= 0.521477	;m_par[11]= 0.933227	 ; return 1;}
 if (fabs(energy-2.15)<1e-4)	{m_par[0]= 0.955803	;m_par[1]=0.00465133	;m_par[2]= 0.0163199	;m_par[3]= 1.09182	;m_par[4]= 9.99761	;m_par[5]= 0.939027	;m_par[6]= 1.07041	;m_par[7]= 0.00650751	;m_par[8]= 0.0260943	;m_par[9]= 2	;m_par[10]= 0.316106	;m_par[11]= 0.878499	 ; return 1;}
 if (fabs(energy-2.175)<1e-4)	{m_par[0]= 0.968996	;m_par[1]=0.00432316	;m_par[2]= 0.00732269	;m_par[3]= 1.36759	;m_par[4]= 3.62356	;m_par[5]= 0.6	;m_par[6]= 1.08268	;m_par[7]= 0.0064486	;m_par[8]= 0.0799998	;m_par[9]= 1.99999	;m_par[10]= 0.323123	;m_par[11]= 0.901922	 ; return 1;}
 if (fabs(energy-2.2)<1e-4)	{m_par[0]= 0.983023	;m_par[1]=0.00464484	;m_par[2]= 0.0101405	;m_par[3]= 2.66973	;m_par[4]= 1	;m_par[5]= 0.693302	;m_par[6]= 1.09523	;m_par[7]= 0.00692449	;m_par[8]= 0.0404309	;m_par[9]= 1.84652	;m_par[10]= 0.641975	;m_par[11]= 0.918789	 ; return 1;}
 if (fabs(energy-2.2324)<1e-4)	{m_par[0]= 1.00066	;m_par[1]=0.00795524	;m_par[2]= 0.005	;m_par[3]= 1.52041	;m_par[4]= 9.99945	;m_par[5]= 0.6	;m_par[6]= 1.11116	;m_par[7]= 0.0053552	;m_par[8]= 0.011222	;m_par[9]= 1.3118	;m_par[10]= 0.429974	;m_par[11]= 0.651126	 ; return 1;}
 if (fabs(energy-2.3094)<1e-4)	{m_par[0]= 1.04347	;m_par[1]=0.00860883	;m_par[2]= 0.005	;m_par[3]= 1.80503	;m_par[4]= 2.11545	;m_par[5]= 0.6	;m_par[6]= 1.14876	;m_par[7]= 0.00771785	;m_par[8]= 0.0267276	;m_par[9]= 1.94156	;m_par[10]= 0.453565	;m_par[11]= 0.931807	 ; return 1;}
 if (fabs(energy-2.3864)<1e-4)	{m_par[0]= 1.08654	;m_par[1]=0.00597808	;m_par[2]= 0.0112714	;m_par[3]= 1.23222	;m_par[4]= 9.25633	;m_par[5]= 0.694457	;m_par[6]= 1.18824	;m_par[7]= 0.00790824	;m_par[8]= 0.0304545	;m_par[9]= 1.89918	;m_par[10]= 0.521754	;m_par[11]= 0.924465	 ; return 1;}
 if (fabs(energy-2.396)<1e-4)	{m_par[0]= 1.09157	;m_par[1]=0.0062201	;m_par[2]= 0.0116202	;m_par[3]= 2.15054	;m_par[4]= 1.31319	;m_par[5]= 0.733822	;m_par[6]= 1.19304	;m_par[7]= 0.00767405	;m_par[8]= 0.0296603	;m_par[9]= 1.53355	;m_par[10]= 0.97687	;m_par[11]= 0.923091	 ; return 1;}
 if (fabs(energy-2.5)<1e-4)	{m_par[0]= 1.14785	;m_par[1]=0.00590093	;m_par[2]= 0.0159552	;m_par[3]= 4.67433	;m_par[4]= 4.31463	;m_par[5]= 0.6	;m_par[6]= 1.24852	;m_par[7]= 0.00853351	;m_par[8]= 0.0139589	;m_par[9]= 1.12254	;m_par[10]= 0.999999	;m_par[11]= 0.99997	 ; return 1;}
 if (fabs(energy-2.6444)<1e-4)	{m_par[0]= 1.22666	;m_par[1]=0.0075946	;m_par[2]= 0.013568	;m_par[3]= 2.50572	;m_par[4]= 1	;m_par[5]= 0.6	;m_par[6]= 1.31725	;m_par[7]= 0.00840415	;m_par[8]= 0.0405968	;m_par[9]= 1.53972	;m_par[10]= 0.985816	;m_par[11]= 0.927726	 ; return 1;}
 if (fabs(energy-2.6464)<1e-4)	{m_par[0]= 1.22784	;m_par[1]=0.00704253	;m_par[2]= 0.0163882	;m_par[3]= 2.36911	;m_par[4]= 1.00181	;m_par[5]= 0.660468	;m_par[6]= 1.31801	;m_par[7]= 0.00791839	;m_par[8]= 0.0262927	;m_par[9]= 1.50979	;m_par[10]= 0.99999	;m_par[11]= 0.818353	 ; return 1;}
 if (fabs(energy-2.7)<1e-4)	{m_par[0]= 1.25556	;m_par[1]=0.00844267	;m_par[2]= 0.00843719	;m_par[3]= 4.88274	;m_par[4]= 3.29246	;m_par[5]= 0.649839	;m_par[6]= 1.34404	;m_par[7]= 0.0101641	;m_par[8]= 0.005	;m_par[9]= 1.23024	;m_par[10]= 1	;m_par[11]= 0.683616	 ; return 1;}
 if (fabs(energy-2.8)<1e-4)	{m_par[0]= 1.3093	;m_par[1]=0.00649012	;m_par[2]= 0.02	;m_par[3]= 2.39573	;m_par[4]= 1.85898	;m_par[5]= 0.600002	;m_par[6]= 1.394	;m_par[7]= 0.00876493	;m_par[8]= 0.027531	;m_par[9]= 1.24628	;m_par[10]= 0.999996	;m_par[11]= 0.810064	 ; return 1;}
 if (fabs(energy-2.9)<1e-4)	{m_par[0]= 1.36311	;m_par[1]=0.0137545	;m_par[2]= 0.00749657	;m_par[3]= 1.52293	;m_par[4]= 6.95113	;m_par[5]= 0.600002	;m_par[6]= 1.44441	;m_par[7]= 0.0103404	;m_par[8]= 0.0377723	;m_par[9]= 1.70034	;m_par[10]= 0.99999	;m_par[11]= 0.87429	 ; return 1;}
 if (fabs(energy-2.95)<1e-4)	{m_par[0]= 1.39036	;m_par[1]=0.014	;m_par[2]= 0.0062766	;m_par[3]= 1.32192	;m_par[4]= 9.99948	;m_par[5]= 0.600001	;m_par[6]= 1.46969	;m_par[7]= 0.0116659	;m_par[8]= 0.005	;m_par[9]= 1.67524	;m_par[10]= 1	;m_par[11]= 0.940405	 ; return 1;}
 if (fabs(energy-2.981)<1e-4)	{m_par[0]= 1.406	;m_par[1]=0.00928873	;m_par[2]= 0.00513155	;m_par[3]= 1.00056	;m_par[4]= 7.27885	;m_par[5]= 0.945569	;m_par[6]= 1.48475	;m_par[7]= 0.0105223	;m_par[8]= 0.0338782	;m_par[9]= 1.71679	;m_par[10]= 0.999992	;m_par[11]= 0.865187	 ; return 1;}
 if (fabs(energy-3)<1e-4)	{m_par[0]= 1.4168	;m_par[1]=0.0139919	;m_par[2]= 0.0102106	;m_par[3]= 3.75874	;m_par[4]= 4.66468	;m_par[5]= 0.603467	;m_par[6]= 1.49554	;m_par[7]= 0.0119663	;m_par[8]= 0.005	;m_par[9]= 1.58379	;m_par[10]= 0.99999	;m_par[11]= 0.947467	 ; return 1;}
 if (fabs(energy-3.02)<1e-4)	{m_par[0]= 1.42553	;m_par[1]=0.0127171	;m_par[2]= 0.00967	;m_par[3]= 1.53159	;m_par[4]= 9.85947	;m_par[5]= 0.610111	;m_par[6]= 1.5043	;m_par[7]= 0.0100886	;m_par[8]= 0.0251589	;m_par[9]= 1.5	;m_par[10]= 0.999991	;m_par[11]= 0.759624	 ; return 1;}
 if (fabs(energy-3.08)<1e-4)	{m_par[0]= 1.45739	;m_par[1]=0.0094642	;m_par[2]= 0.02	;m_par[3]= 4.02917	;m_par[4]= 3.67045	;m_par[5]= 0.6	;m_par[6]= 1.53437	;m_par[7]= 0.011389	;m_par[8]= 0.0423672	;m_par[9]= 1.17153	;m_par[10]= 1.88202	;m_par[11]= 0.921112	 ; return 1;}
 return false;
}
bool getParsMC(double energy)
{
if (fabs(energy-2)<1e-4)	{m_par[0]= 0.869729	 ;m_par[1]=0.0040562	 ;m_par[2]=0.0087734	 ;m_par[3]=3.10275	;m_par[4]= 1.59408	;m_par[5]= 0.873228	 ;m_par[6]=0.992404	;m_par[7]= 0.00395914	 ;m_par[8]=0.0050307	 ;m_par[9]=0.1	 ;m_par[10]=0.749329	 ;m_par[11]=1	 ; return 1;}
 if (fabs(energy-2.05)<1e-4)	{m_par[0]= 0.898391	 ;m_par[1]=0.0043809	 ;m_par[2]=0.0107161	 ;m_par[3]=3.29088	;m_par[4]= 1.35412	;m_par[5]= 0.923564	 ;m_par[6]=1.01754	;m_par[7]= 0.00425949	 ;m_par[8]=0.0798512	 ;m_par[9]=0.113476	 ;m_par[10]=0.111696	 ;m_par[11]=0.778929	 ; return 1;}
 if (fabs(energy-2.1)<1e-4)	{m_par[0]= 0.92691	 ;m_par[1]=0.0045242	 ;m_par[2]=0.0105272	 ;m_par[3]=3.58752	;m_par[4]= 1.00783	;m_par[5]= 0.908488	 ;m_par[6]=1.04273	;m_par[7]= 0.00455985	 ;m_par[8]=0.0799988	 ;m_par[9]=0.100005	 ;m_par[10]=0.948016	 ;m_par[11]=0.600111	 ; return 1;}
 if (fabs(energy-2.15)<1e-4)	{m_par[0]= 0.955126	 ;m_par[1]=0.0048575	 ;m_par[2]=0.0125467	 ;m_par[3]=2.95272	;m_par[4]= 1.59553	;m_par[5]= 0.932985	 ;m_par[6]=1.06779	;m_par[7]= 0.0048602 	 ;m_par[8]=0.0799997	 ;m_par[9]=0.509777	 ;m_par[10]=0.10004	 ;m_par[11]=0.798717	 ; return 1;}
 if (fabs(energy-2.175)<1e-4)	{m_par[0]= 0.969164	 ;m_par[1]=0.0049292	 ;m_par[2]=0.0130111	 ;m_par[3]=3.19578	;m_par[4]= 1.47553	;m_par[5]= 0.937991	 ;m_par[6]=1.08418	;m_par[7]= 0.00501039	 ;m_par[8]=0.0799997	 ;m_par[9]=0.1	 ;m_par[10]=0.99999	 ;m_par[11]=0.603054	 ; return 1;}
 if (fabs(energy-2.2)<1e-4)	{m_par[0]= 0.983198	 ;m_par[1]=0.0050458	 ;m_par[2]=0.0127866	 ;m_par[3]=3.08916	;m_par[4]= 1.48143	;m_par[5]= 0.941396	 ;m_par[6]=1.09291	;m_par[7]= 0.00516056	 ;m_par[8]=0.0799919	 ;m_par[9]=0.100007	 ;m_par[10]=0.999996	 ;m_par[11]=0.807618	 ; return 1;}
 if (fabs(energy-2.2324)<1e-4)	{m_par[0]= 1.00143	 ;m_par[1]=0.0062354	 ;m_par[2]=0.005	 ;m_par[3]=2.06359	;m_par[4]= 3.77784	;m_par[5]= 0.6	 ;m_par[6]=1.10918	;m_par[7]= 0.0053552 	 ;m_par[8]=0.0799991	 ;m_par[9]=0.1	 ;m_par[10]=0.999997	 ;m_par[11]=0.889524	 ; return 1;}
 if (fabs(energy-2.3094)<1e-4)	{m_par[0]= 1.04414	 ;m_par[1]=0.0067412	 ;m_par[2]=0.005	 ;m_par[3]=2.01489	;m_par[4]= 4.43378	;m_par[5]= 0.6	 ;m_par[6]=1.14785	;m_par[7]= 0.00581775	 ;m_par[8]=0.0799998	 ;m_par[9]=0.1	 ;m_par[10]=0.999926	 ;m_par[11]=0.843378	 ; return 1;}
 if (fabs(energy-2.3864)<1e-4)	{m_par[0]= 1.08657	 ;m_par[1]=0.0071719	 ;m_par[2]=0.005	 ;m_par[3]=2.3555	;m_par[4]= 2.45342	;m_par[5]= 0.6	 ;m_par[6]=1.18651	;m_par[7]= 0.0062803 	 ;m_par[8]=0.0799987	 ;m_par[9]=0.1	 ;m_par[10]=0.999915	 ;m_par[11]=0.855787	 ; return 1;}
 if (fabs(energy-2.396)<1e-4)	{m_par[0]= 1.09184	 ;m_par[1]=0.0074491	 ;m_par[2]=0.005	 ;m_par[3]=2.30098	;m_par[4]= 3.46626	;m_par[5]= 0.6	 ;m_par[6]=1.19133	;m_par[7]= 0.00633797	 ;m_par[8]=0.079988	 ;m_par[9]=0.1	 ;m_par[10]=0.722438	 ;m_par[11]=0.915832	 ; return 1;}
 if (fabs(energy-2.5)<1e-4)	{m_par[0]= 1.14888	 ;m_par[1]=0.0083017	 ;m_par[2]=0.005	 ;m_par[3]=2.17074	;m_par[4]= 4.131	;m_par[5]= 0.6	 ;m_par[6]=1.24352	;m_par[7]= 0.00696271	 ;m_par[8]=0.0765445	 ;m_par[9]=0.100002	 ;m_par[10]=0.999997	 ;m_par[11]=0.600002	 ; return 1;}
 if (fabs(energy-2.6444)<1e-4)	{m_par[0]= 1.22704	 ;m_par[1]=0.0092900	 ;m_par[2]=0.0055570	 ;m_par[3]=2.50932	;m_par[4]= 2.90496	;m_par[5]= 0.6	 ;m_par[6]=1.31597	;m_par[7]= 0.00783015	 ;m_par[8]=0.08	 ;m_par[9]=0.1	 ;m_par[10]=0.999993	 ;m_par[11]=0.722241	 ; return 1;}
 if (fabs(energy-2.6464)<1e-4)	{m_par[0]= 1.22814	 ;m_par[1]=0.0092702	 ;m_par[2]=0.0054301	 ;m_par[3]=2.0624	;m_par[4]= 6.1071	;m_par[5]= 0.6	 ;m_par[6]=1.31697	;m_par[7]= 0.00784216	 ;m_par[8]=0.0799993	 ;m_par[9]=0.131213	 ;m_par[10]=0.999999	 ;m_par[11]=0.904941	 ; return 1;}
 if (fabs(energy-2.7)<1e-4)	{m_par[0]= 1.25689	 ;m_par[1]=0.0071972	 ;m_par[2]=0.0142876	 ;m_par[3]=2.74333	;m_par[4]= 1.93963	;m_par[5]= 0.860784	 ;m_par[6]=1.34386	;m_par[7]= 0.00816415	 ;m_par[8]=0.0433022	 ;m_par[9]=0.189595	 ;m_par[10]=0.980516	 ;m_par[11]=0.658781	 ; return 1;}
 if (fabs(energy-2.8)<1e-4)	{m_par[0]= 1.31082	 ;m_par[1]=0.0081207	 ;m_par[2]=0.02	 ;m_par[3]=4.54691	;m_par[4]= 3.57218	;m_par[5]= 0.933149	 ;m_par[6]=1.394	;m_par[7]= 0.00876486	 ;m_par[8]=0.0347896	 ;m_par[9]=0.486668	 ;m_par[10]=0.99999	 ;m_par[11]=0.600001	 ; return 1;}
 if (fabs(energy-2.9)<1e-4)	{m_par[0]= 1.36399	 ;m_par[1]=0.0085246	 ;m_par[2]=0.0182617	 ;m_par[3]=4.78015	;m_par[4]= 4.20569	;m_par[5]= 0.915264	 ;m_par[6]=1.44414	;m_par[7]= 0.00936558	 ;m_par[8]=0.0754732	 ;m_par[9]=0.182487	 ;m_par[10]=0.999999	 ;m_par[11]=1	 ; return 1;}
 if (fabs(energy-2.95)<1e-4)	{m_par[0]= 1.39074	 ;m_par[1]=0.0088054	 ;m_par[2]=0.0184032	 ;m_par[3]=2.62066	;m_par[4]= 2.35007	;m_par[5]= 0.901543	 ;m_par[6]=1.46953	;m_par[7]= 0.00966593	 ;m_par[8]=0.0799999	 ;m_par[9]=0.197325	 ;m_par[10]=0.108751	 ;m_par[11]=0.600002	 ; return 1;}
 if (fabs(energy-2.981)<1e-4)	{m_par[0]= 1.40688	 ;m_par[1]=0.0083764	 ;m_par[2]=0.0150285	 ;m_par[3]=3.17317	;m_par[4]= 1.13186	;m_par[5]= 0.8053	 ;m_par[6]=1.48475	;m_par[7]= 0.00985215	 ;m_par[8]=0.08	 ;m_par[9]=0.16804	 ;m_par[10]=0.999995	 ;m_par[11]=0.728436	 ; return 1;}
 if (fabs(energy-3)<1e-4)	{m_par[0]= 1.41706	 ;m_par[1]=0.0084445	 ;m_par[2]=0.013986	 ;m_par[3]=2.8561	;m_par[4]= 1.22928	;m_par[5]= 0.739855	 ;m_par[6]=1.49427	;m_par[7]= 0.00996629	 ;m_par[8]=0.0481069	 ;m_par[9]=0.100024	 ;m_par[10]=0.999987	 ;m_par[11]=0.600005	 ; return 1;}
 if (fabs(energy-3.02)<1e-4)	{m_par[0]= 1.4279	 ;m_par[1]=0.0092258	 ;m_par[2]=0.0179922	 ;m_par[3]=2.80752	;m_par[4]= 3.11775	;m_par[5]= 0.908431	 ;m_par[6]=1.5093	;m_par[7]= 0.0100864 	 ;m_par[8]=0.0799874	 ;m_par[9]=0.168894	 ;m_par[10]=0.102016	 ;m_par[11]=1	 ; return 1;}
 if (fabs(energy-3.08)<1e-4)	{m_par[0]= 1.45944	 ;m_par[1]=0.0124944	 ;m_par[2]=0.0076822	 ;m_par[3]=2.12209	;m_par[4]= 7.8597	;m_par[5]= 0.6	 ;m_par[6]=1.53437	;m_par[7]= 0.0104471 	 ;m_par[8]=0.08	 ;m_par[9]=0.445917	 ;m_par[10]=1	 ;m_par[11]=0.6	 ; return 1;}
return false;
}

double FitSpectrum5(TTree *&dataraw, double beame, const char* namesfx) // fix par
{
 //TF1 fsigma("fsigma","-0.00705519+0.00600716*x",2.0,3.1);
 //double sigma_ini = fsigma.Eval(beame);
   if (!getPars(beame)) return -1;
   //if (!getParsMC(beame)) return -1;
   TCanvas *c1 = new TCanvas("","",800,600);
   
   int nBins=100;
   int Npar;
   double mka = 0.493677;
   double mmu = 0.1057;
   beame = beame/2;
   double peakvalue = sqrt(pow(beame,2)-pow(mka,2));
   double pmu = sqrt(pow(beame,2)-pow(mmu,2));
   std::cout<<"Fitspectrum "<< peakvalue <<"GeV/c" <<std::endl;
   double beamlow=peakvalue-0.2;
   double beamup=peakvalue+0.2;
   
   // save histogram of pmu
  //TH1D *hp = new TH1D("hp","hp",100,beamlow,beamup);
  //dataraw->Draw("x1>>hp");
  //char name[100];
  //sprintf(name,"hp_mcKK_%.4f",beame*2);
  //hp->Write(name);
  //delete hp;
  //return 0;
   //
   // try to use roofit
   //
   RooRealVar x("x1","momentum",peakvalue,beamlow,beamup,"GeV");
   // signal
   RooRealVar mean("mean","mean of gaussian",m_par[0]);
   RooRealVar sigma("sigma","width of gaussian",m_par[1]);
   RooRealVar alpha1("alpha1","#alpha",m_par[3]);
   RooRealVar nnn1("n1","n",m_par[4]);
   RooCBShape cbshape("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);
   RooRealVar sigma2("sigma2","width of gaussian",m_par[2]);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma2);
   RooRealVar frac1("frac1","frac1",m_par[5]);

   RooAddPdf sig("sig","signal",RooArgList(cbshape,gaus),RooArgList(frac1));
   
   // background, mainly di-mu
   RooRealVar meanb("meanb","mean of gaussian",m_par[6]);
   //RooRealVar sigmab("sigmab","width of gaussian",m_par[7]);
   RooRealVar sigmab("sigmab","width of gaussian",m_par[7],m_par[7]-0.001,m_par[7]+0.001);
   //RooRealVar sigmab("sigmab","width of gaussian",sigma_ini, sigma_ini-0.001,sigma_ini+0.001);
   //RooRealVar sigmab("sigmab","width of gaussian",0.005,0.002,0.02);
 //  RooGaussian gausb("gausb","gauss(x,m,s)",x,meanb,sigmab);
   //RooRealVar alphab("alphab","#alpha",m_par[9]);
   RooRealVar alphab("alphab","#alpha",m_par[9]+1,m_par[9]+0.05,m_par[9]+10);
   //RooRealVar alphab("alphab","#alpha",0.5,0.1,1.5);
   //RooRealVar nnnb("nnnb","n",m_par[10]);
   RooRealVar nnnb("nnnb","n",m_par[10]+2,m_par[10]+0.2,m_par[10]+10);
   RooCBShape cbshapeb("cbshapeb","crystal ball",x,meanb,sigmab,alphab,nnnb);
   //RooRealVar meane("meane","mean of gaussian",beame,beame-0.005,beame+0.003);
   RooRealVar sigmae("sigmae","width of gaussian",m_par[8]);
   RooGaussian gause("gause","gauss(x,m,s)",x,meanb,sigmae);
   RooRealVar frac2("frac2","frac2",m_par[11],m_par[11]-0.3,1);
   RooAddPdf bck("bck","signal",RooArgList(cbshapeb,gause),RooArgList(frac2));
   
   RooRealVar signal("signal"," ",1000,0,1000000000);//event number
   RooRealVar background("background"," ",20000,0,100000000);
 
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   
  
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit

   char tmpchr[100];
   sprintf(tmpchr,"data_2trk_%s",namesfx);
   xframe = x.frame(Title(""));
   c1->SetMargin(0.13,0.1,0.13,0.1);
   TGaxis::SetMaxDigits(3);
   xframe->SetTitle("");
   xframe->GetXaxis()->SetTitle("p (GeV/c)");
   xframe->GetXaxis()->SetNdivisions(505);
   xframe->GetXaxis()->SetTitleSize(0.05);
   xframe->GetXaxis()->SetLabelSize(0.05);
   xframe->GetXaxis()->SetTitleOffset(1.1);
   xframe->GetYaxis()->SetNdivisions(505);
   xframe->GetYaxis()->SetTitleSize(0.05);
   xframe->GetYaxis()->SetTitleOffset(1.25);
   xframe->GetYaxis()->SetLabelSize(0.05);
   xframe->GetYaxis()->SetNoExponent(0);  
   
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(sig,bck),RooArgList(signal,background));
   Npar = 2;
   x.setRange("sigragi",peakvalue-0.1,beame+0.012);
   //sum->fitTo(*dataset);
   sum->fitTo(*dataset);
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bck),LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
   //sum->plotOn(xframe);
   xframe->Draw();
   TPaveText *pt = new TPaveText(0.15,0.65,0.45,0.90,"BRNDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(4000);
   pt->SetTextAlign(12);
   pt->SetTextFont(42);
   pt->SetTextSize(0.035);
   sprintf(tmpchr,"#mu = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#sigma = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
   pt->AddText(tmpchr);
   sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
   pt->AddText(tmpchr);
   pt->Draw();
   sprintf(tmpchr,"p_spectrum_%s",namesfx);
   c1->SetName(tmpchr);
   c1->Write();
   sprintf(tmpchr,"p_spectrum_%s.pdf",namesfx);
   c1->Print(tmpchr);

 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   ofstream fitpar("fitpar",std::ios::app);
 //fitpar<<" ene = "<< beame*2 <<"\t sig mean = "<< mean.getVal() << "\t sig sigma = "<< sigma.getVal();
 //fitpar<<"\t bkg mean = " << meanb.getVal() << "\t bkg sigma = " << sigmab.getVal();
 //fitpar<<"\t sigNo = " << signal.getVal() << "\t sigNoE = " << signal.getError();
 //fitpar<<"\t bckNo = " << background.getVal() << "\t bckNoE = " << background.getError();
 //fitpar<< "\t chi2 = " << xframe->chiSquare(Npar);
 //fitpar<< std::endl;
   fitpar<< beame*2 <<"\t";
   fitpar<< signal.getVal() <<" +/- "<< signal.getError()<<"\t";
   fitpar<< background.getVal() <<" +/- "<< background.getError();
   
 //fitpar<< mean.getVal() <<" ";
 //fitpar<< sigma.getVal() <<" ";
 //fitpar<< sigma2.getVal() <<" ";
 //fitpar<< alpha1.getVal()<<"  ";
 //fitpar<< nnn1.getVal()<<"  ";
 //fitpar<< frac1.getVal()<<"  ";
 //
 //fitpar<< meanb.getVal()<<"  ";
 //fitpar<< sigmab.getVal()<<"  ";
 //fitpar<< sigmae.getVal()<<"  ";
 //fitpar<< alphab.getVal()<<"  ";
 //fitpar<< nnnb.getVal()<<"  ";
 //fitpar<< frac2.getVal()<<"  ";
   
   fitpar<< endl;
   //fitpar<<"\t e mean = " << meane.getVal() << "\t e sigma = " << sigmae.getVal()<<std::endl;
   
   double intel3 = mean.getVal()-3*sigma.getVal();
   double inteu3 = mean.getVal()+3*sigma.getVal();
   double intel5 = mean.getVal()-5*sigma.getVal();
   double inteu5 = mean.getVal()+5*sigma.getVal();
   //RooRealVar xint("xint","xint",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame);
   //x.setRange("sigragi",peakvalue-0.1,beame+0.01);
   x.setRange("sigrag3",intel3, inteu3);
   x.setRange("sigrag5",intel5, inteu5);
   RooAbsReal* intsigi = sig.createIntegral( x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intbcki = bck.createIntegral(x,  NormSet(x),  Range("sigragi"));
   RooAbsReal* intsig3 = sig.createIntegral( x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intbck3 = bck.createIntegral(x,  NormSet(x),  Range("sigrag3"));
   RooAbsReal* intsig5 = sig.createIntegral( x,  NormSet(x),  Range("sigrag5"));
   RooAbsReal* intbck5 = bck.createIntegral(x,  NormSet(x),  Range("sigrag5"));
   double signalN3 =   signal.getVal()*intsig3->getVal();
   double backN3 = background.getVal()*intbck3->getVal();
   double signalN5 =   signal.getVal()*intsig5->getVal();
   double backN5 = background.getVal()*intbck5->getVal();
 //double signalN3nom =   signal.getVal()*intsig3->getVal()/intsigi->getVal();
 //double backN3nom = background.getVal()*intbck3->getVal()/intbcki->getVal();
 //double signalN5nom =   signal.getVal()*intsig5->getVal()/intsigi->getVal();
 //double backN5nom = background.getVal()*intbck5->getVal()/intbcki->getVal();
   ofstream angSNR("angSNR.dat",std::ios::app);
   angSNR << namesfx <<"\t"<<signalN3 << "\t" << backN3 <<"\t"<<signalN5 << "\t" << backN5 << std::endl;
 //epSNR << "\t"<<signalN3nom << "\t" << backN3nom <<"\t"<<signalN5nom << "\t" << backN5nom<< std::endl;
   std::cout<< "Total signal int is "<< intsigi->getVal() <<" Total bck int is "<< intbcki->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig3->getVal() <<" Total bck int is "<< intbck3->getVal()<<std::endl;
   std::cout<< "Total signal int is "<< intsig5->getVal() <<" Total bck int is "<< intbck5->getVal()<<std::endl;

   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
}




double GetEnergy(int run)
{
  run = abs(run);
  if (run>=39335 && run<=39618) return 3.08;
  if (run>=39711 && run<=39738) return 3.02;
  if (run>=39680 && run<=39710) return 3.00;
  if (run>=39651 && run<=39679) return 2.981;
  if (run>=39619 && run<=39650) return 2.95;
  
  if (run>=39775 && run<=40069) return 2.90;
  if (run>=40128 && run<=40296) return 2.6444;
  if (run>=40300 && run<=40435) return 2.6464;
  if (run>=40436 && run<=40439) return 2.70;
  if (run>=40440 && run<=40443) return 2.80;
  if (run>=40459 && run<=40769) return 2.396;
  if (run>=40771 && run<=40776) return 2.5;
  if (run>=40777 && run<=40804) return 2.6444;//separated beam
  if (run>=40806 && run<=40951) return 2.3864;
  if (run>=40989 && run<=41121) return 2.2;
  if (run>=41122 && run<=41239) return 2.2324;
  if (run>=41240 && run<=41411) return 2.3094;
  if (run>=41416 && run<=41532) return 2.175;
  if (run>=41533 && run<=41570) return 2.15;
  if (run>=41588 && run<=41727) return 2.1;

  if (run>=41729 && run<=41909) return 2.0;
  if (run>=41911 && run<=41958) return 2.05;
  if (run>=41959 && run<=41999) return 2.2324; // separated beam
  
  if (run>=42004) return 2.125;
  return -1;
}


const char* getPureName(const char* name)
{
  int pos=0;
  for (int i=0;;i++){
    if (name[i]=='\0') return &name[pos];
    if (name[i]=='/') pos = i+1;
    if (i>10000) return NULL;
  } 
}

const char* getPureName2(const char* name)
{
  char *name1 = new char[1000];
  int pos=0;
  int ppos=0;

  for (int i=0;;i++){
    if (name[i]=='\0') break; //return &name[pos];
    if (name[i]=='/') pos = i+1;
    if (name[i]=='.') ppos = i;
    if (i>1000) return NULL;
  }
  //std::cout<<"/ pos is "<< pos <<", . pos is "<< ppos <<std::endl;
  for (int i=0;i<ppos-pos;i++){
    name1[i] = name[pos+i];
  }
  name1[ppos-pos]='\0';

  return name1;
}


double getEne(const char* name)
{
  name = getPureName(name);
  std::string str = name;
  int _pos=0;
  int ppos=0;
  for (int i=0;;i++){
    if (name[i]=='\0') break;
    if (name[i]=='_') _pos=i;
    if (name[i]=='.') ppos=i;
  }
  int len = ppos-_pos-1;
  string vstr=str.substr(_pos+1,len);
  if (atof(vstr.c_str())/10000>1) return atof(vstr.c_str())/10000;

  else {// find index number in _*_
    int _pos2[2]={0,0};
    int j=0;
    for (int i=0;;i++){
      if (name[i]=='_') {_pos2[j]=i; j++;}
      if (j>=2) break;
    }
    if (j!=2) return -1;
    string vstr = str.substr(_pos2[0]+1,3);
    int idx = atoi(vstr.c_str());
    std::cout<<vstr<<'\t'<<idx<<std::endl;
    switch (idx) {
      case 1:    return  2.0;   
      case 2:    return  2.05;  
      case 3:    return  2.1;   
      case 4:    return  2.15;  
      case 5:    return  2.175; 
                                
      case 6:    return  2.2;   
      case 7:    return  2.2324;
      case 8:    return  2.3094;
      case 9:    return  2.3864;
      case 10:   return  2.396; 
      case 11:   return  2.5;   
      case 12:   return  2.6444;
      case 13:   return  2.6464;
      case 14:   return  2.7;   
      case 15:   return  2.8;   
      case 16:   return  2.9;   
      case 17:   return  2.95;  
      case 18:   return  2.981; 
      case 19:   return  3.0;   
      case 20:   return  3.02;  
                                
      case 21:   return  3.08;  
    }
  }
  return -2;
}

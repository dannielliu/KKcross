#include "TFile.h"
//#include "TFolder.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TLorentzVector.h"
#include "TH2D.h"
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
const char* getPureName(const char* name);
const char* getPureName2(const char* name);
double getEne(const char* name);


int main(int argc, char** argv)
{
  //if (argc==1) Ana(argv[1]);
  //if (argc>1)  Ana(argv[1],argv[2]);
  //else Ana();
  TFile *fileout = new TFile("output/output_mc.root","update");
  for (int i=1; i<argc; i++){
    Ana(argv[i],0,fileout);
  }
  fileout->Write();
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

  double ene[21];
  double pcut[21];
  double epcut[21];
  double m_pcut;
  double m_epcut;
  ifstream cutin("cutpar");
  char line[1000];
  if (cutin.is_open()) cutin.getline(line,1000);
  int iene=0;
  while (!cutin.eof()){
    cutin.getline(line,1000);
    istringstream iss;
    iss.str(line);
    iss>>ene[iene]>>pcut[iene]>>epcut[iene];
    iene++;
    if (iene==21) break;
  }
  
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
  double pid1[5];
  double pid2[5];
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
  tree->SetBranchAddress("pid1",pid1);
  tree->SetBranchAddress("pid2",pid2);
  double pi = TMath::Pi();
  double mka = 0.493677;
  //double mka = 0.13957;
  //double mka = 0.1057;

  tree->GetEntry(0);
  double Ebeam = GetEnergy(run);
  if (Ebeam<1.0) Ebeam = getEne(filename);

  for (int iene=0;iene<21;iene++){
    if (Ebeam == ene[iene]) {
      m_pcut = pcut[iene];
      m_epcut = epcut[iene];
      break;
    }
  }
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
  
  fileout->cd();
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
  vars->Branch("run",&run,"run/I");
  vars->Branch("indexmc",&idxmc,"indexmc/I");
  vars->Branch("pdgid",pdgid,"pdgid[100]/I");
  vars->Branch("motheridx",motheridx,"motheridx[100]/I");
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
  TTree *datasel2 = new TTree("datasel2","datasel2");
  datasel2->Branch("x1",&p2,"x1/D");
  datasel2->Branch("mass",&mass,"mass/D");

  int nentries = tree->GetEntries();
  double pexp = sqrt(pow(Ebeam/2,2)-pow(mka,2));
  std::cout<<"Total entries is " << nentries << std::endl;
  std::cout<<"Beam energy is " << Ebeam << std::endl;
  std::cout<<"Expected K momentum is " << pexp << std::endl;
  int count0=0,count1=0,count2=0,count3=0;
  int count[10]={0};
  int tagmc=0;

  TH1D *ldis = new TH1D("ldis","ldis",100,-50,50);
  TH1D *lcos1 = new TH1D("lcos1","lcos1",3,0,3);
  TH2D *hlcos = new TH2D("hlcos","hlcos",1000,-1,1,3,0,3);
  TH2D *hlcos1 = new TH2D("hlcos1","hlcos1",1000,-1,1,30,0,30);
  TH2D *hlcos2 = new TH2D("hlcos2","hlcos2",1000,-1,1,30,0,30);
  TH2D *hlcos3 = new TH2D("hlcos3","hlcos3",1000,-1,1,300,0,300);

  for (int ien=0;ien<nentries;ien++){
    tree->GetEntry(ien);
    
    int isrtag=0;
    for (int j=0;j<idxmc;j++){
      if (pdgid[j] == 22 && motheridx[j]==j) { tagmc++; isrtag=1;}
    }
    //if (isrtag==1) continue;
    //if (run<41122) continue;
    //if (costheta>-0.995) continue;
    //if (nneu>0) continue;
    if (emcstatusShort!=0x3 && emcstatusInt!=0x3) continue;
    if (ntof1<1) continue;
    if (ntof2<1) continue;
    for (int i=0;i<ntof1;i++){
      if (tofl1[i]==1) tof11=tof1[i];
    }
    for (int i=0;i<ntof2;i++){
      if (tofl2[i]==1) tof21=tof2[i];
    }
    if (pid1[3]<pid1[0] || pid1[3]<pid1[1] || pid1[3]<pid1[2] || pid1[3]<pid1[4]) continue;
    if (pid2[3]<pid2[0] || pid2[3]<pid2[1] || pid2[3]<pid2[2] || pid2[3]<pid2[4]) continue;
    //if (tof11<-999 || tof21<-999) continue;
    count0++;

    kap.SetVectMag(TVector3(kappx,kappy,kappz),mka);
    kam.SetVectMag(TVector3(kampx,kampy,kampz),mka);
    
    // angular information in lab coordinate
    costheta1 = kap.CosTheta();
    costheta2 = kam.CosTheta();
    
    for (int i=0;i<ntof1;i++){
      //if (tofl1[i]==1) tof11=tof1[i];
      hlcos->Fill(costheta1,tofl1[i]);
    }
    if (ntof1==1) hlcos1->Fill(costheta1,tofl1[0]);
    if (ntof1==2) hlcos2->Fill(costheta1,tofl1[0]*10+tofl1[1]);
    if (ntof1==3) hlcos3->Fill(costheta1,tofl1[0]*100+tofl1[1]*10+tofl1[2]);

    // momentum information in cms coordinalte
    kap.Boost(-0.011,0,0);
    kam.Boost(-0.011,0,0);
    costheta = kap.Vect().Dot(kam.Vect())/(kap.Vect().Mag()*kam.Vect().Mag());
    mass = (kap+kam).M();
    totp = (kap+kam).Rho();
    p1 = kap.Rho();
    p2 = kam.Rho();
 
    vars->Fill();

    // select candidate, p1 dis
    // if (count0%100000==0) std::cout<<"epratio1 is "<< epratio1<<std::endl;
    if (fabs(costheta1)<0.93){ count[4]++;
    if (fabs(costheta2)<0.93){ count[5]++;
    if (costheta<cos(pi*179/180)){ count[2]++;
    if (epratio1<m_epcut){ count[0]++;
    if (epratio2<m_epcut){ count[1]++;
    //if (epratio2<0.85){ count[1]++;
    //if (totp<0.05){ count[3]++;
    if (fabs(tof11-tof21)<3){
       count3++;
       //if (fabs(p1-pexp)<0.08) {/*datasel1->Fill();*/ count1++; // p1 dis
       //if (fabs(p2-pexp)<0.03) {/*datasel2->Fill();*/ count2++; // p2 dis
       if (p2<m_pcut) {/*datasel2->Fill();*/ count2++; // p2 dis
         datasel1->Fill();
       }
       //}
    }
    }
   // }
    }
    }
    }
    }
  }

  ofstream cutflow("cutflow",std::ios::app);
  cutflow<<filename<<std::endl;
  std::cout<<"data selction 1 size is "<<datasel1->GetEntries()<<std::endl;
  //std::cout<<"data selction 2 size is "<<datasel2->GetEntries()<<" "<<count2<<std::endl;
  cutflow<<"Initial size      :"<<nentries<<std::endl;
  cutflow<<"Tagged ISR evt    :"<< tagmc  <<std::endl;
  cutflow<<"After cos1<0.93   :"<<count[4]<<std::endl;
  cutflow<<"After cos2<0.93   :"<<count[5]<<std::endl;
  cutflow<<"After cos<179     :"<<count[2]<<std::endl;
  cutflow<<"After ep1         :"<<count[0]<<std::endl;
  cutflow<<"After ep2         :"<<count[1]<<std::endl;
  //cutflow<<"After totp<0.05   :"<<count[3]<<std::endl;
  cutflow<<"After dtof<3      :"<<count3<<std::endl;
  //cutflow<<"p1-exp<0.08 for p1:"<<count1<<std::endl;
  cutflow<<"p2-exp<3 sigma    :"<<count2<<std::endl;
  //if (count2==0) return -3;

  double sigp1=FitSpectrum(datasel1, Ebeam, pureName);
  //double sigm1=FitSpectrum2(datasel1, Ebeam, pureName);
//double sigp2=FitSpectrum(datasel2, Ebeam,"p2");
//double sigm2=FitSpectrum2(datasel2, Ebeam,"p2");
  cutflow<<"p1 fit signal     :"<<sigp1 <<std::endl;
  //cutflow<<"m1 fit signal     :"<<sigm1 <<std::endl;
//cutflow<<"p2 fit signal     :"<<sigp2 <<std::endl;
//cutflow<<"m2 fit signal     :"<<sigm2 <<std::endl;

  char name2[100];
  sprintf(name2,"hlcos");
  dir->WriteTObject(hlcos,name2);
  sprintf(name2,"hlcos1");
  dir->WriteTObject(hlcos1,name2);
  sprintf(name2,"hlcos2");
  dir->WriteTObject(hlcos2,name2);
  sprintf(name2,"hlcos3");
  dir->WriteTObject(hlcos3,name2);
  //fileout->Close();
  dir->WriteTObject(vars);
  //tmpfile->WriteTObject(vars);

  
  
  // cut study
  //sprintf(name2,"%s.root",name);
  //TFile *file2 = new TFile(name2,"recreate");
  TH1D* hpdis = new TH1D("hpdis","p distribution",1000,0,2.0);
  TH2D* h2pdis = new TH2D("h2pdis","p distribution",1000,0,2.0,1000,0,2.0);
  char cut[1000];
  
  sprintf(cut,"");
  vars->Draw("p1>>hpdis",cut);
  vars->Draw("p1:p2>>h2pdis",cut);
  sprintf(name2,"hpNocut");
  dir->WriteTObject(hpdis,name2);
  sprintf(name2,"h2pNocut");
  dir->WriteTObject(h2pdis,name2);
  
  sprintf(cut,"abs(p2-%f)<0.03",pexp);
  vars->Draw("p1>>hpdis",cut);
  //vars->Draw("p1:p2>>h2pdis",cut);
  sprintf(name2,"hpcut1p");
  dir->WriteTObject(hpdis,name2);
  //sprintf(name2,"h2pcut1p");
  //dir->WriteTObject(h2pdis,name2);
  
  sprintf(cut,"abs(p2-%f)<0.03 & abs(costheta1)<0.8 & abs(costheta2)<0.8",pexp);
  vars->Draw("p1>>hpdis",cut);
  sprintf(cut,"abs(costheta1)<0.8 & abs(costheta2)<0.8");
  vars->Draw("p1:p2>>h2pdis",cut);
  sprintf(name2,"hpcutcosi");
  dir->WriteTObject(hpdis,name2);
  sprintf(name2,"h2pcutcosi");
  dir->WriteTObject(h2pdis,name2);
  
  sprintf(cut,"abs(p2-%f)<0.03 & abs(costheta1)<0.8 & abs(costheta2)<0.8 & costheta<-0.99",pexp);
  vars->Draw("p1>>hpdis",cut);
  sprintf(cut,"abs(costheta1)<0.8 & abs(costheta2)<0.8 & costheta<-0.99");
  vars->Draw("p1:p2>>h2pdis",cut);
  sprintf(name2,"hpcutcos");
  dir->WriteTObject(hpdis,name2);
  sprintf(name2,"h2pcutcos");
  dir->WriteTObject(h2pdis,name2);
  
  sprintf(cut,"abs(p2-%f)<0.03 & abs(costheta1)<0.8 & abs(costheta2)<0.8  & costheta<-0.99 & ep1<0.85 & ep2<0.85",pexp);
  vars->Draw("p1>>hpdis",cut);
  sprintf(cut,"abs(costheta1)<0.8 & abs(costheta2)<0.8  & costheta<-0.99& ep1<0.85 & ep2<0.85");
  vars->Draw("p1:p2>>h2pdis",cut);
  sprintf(name2,"hpcutep");
  dir->WriteTObject(hpdis,name2);
  sprintf(name2,"h2pcutep");
  dir->WriteTObject(h2pdis,name2);
  
  sprintf(cut,"abs(p2-%f)<0.03 & abs(costheta1)<0.8 & abs(costheta2)<0.8 & costheta<-0.99 & ep1<0.85 & ep2<0.85 & abs(tof1-tof2)<4",pexp);
  vars->Draw("p1>>hpdis",cut);
  sprintf(cut,"abs(costheta1)<0.8 & abs(costheta2)<0.8 & costheta<-0.99 & ep1<0.85 & ep2<0.85 & abs(tof1-tof2)<0.7");
  vars->Draw("p1:p2>>h2pdis",cut);
  sprintf(name2,"hpcutTOF");
  dir->WriteTObject(hpdis,name2);
  sprintf(name2,"h2pcutTOF");
  dir->WriteTObject(h2pdis,name2);
  
  dir->Write();
  //file2->Delete("vars;*");
  //fileout->Delete("vars;*");
  

  delete ldis;
  delete lcos1;
  delete hlcos;
  delete hlcos1;
  delete hlcos2;
  delete hlcos3;
  delete hpdis;
  delete h2pdis;

  delete datasel1;
  delete datasel2;
  delete []pureName;
  delete vars;
  //delete tmpfile;
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
   RooRealVar mean("mean","mean of gaussian",peakvalue,peakvalue-0.003,peakvalue+0.003);
   RooRealVar sigma("sigma","width of gaussian",0.005,0.004,0.012);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   //RooRealVar mean2("mean2","mean of gaussian",beame,peakvalue,beame+0.1);
   RooRealVar mean2("mean2","mean of gaussian",peakvalue+0.2,peakvalue+0.05,3.0);
   RooRealVar sigma2("sigma2","width of gaussian",0.09,0.03,0.15);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean2,sigma2);
   
   RooRealVar co1("co1","coefficient #1",   0 ,-100.,100.);
   //RooRealVar co2("co2","coefficient #1",   0 ,-100.,100.);
 //RooRealVar co3("co3","coefficient #1",0.1,-100.,100.);
   RooChebychev ground("ground","background",x,RooArgList(co1));
   RooRealVar signal("signal"," ",1000,0,1000000000);//event number
   RooRealVar background("background"," ",20000,0,100000000);
 //RooRealVar a0("a0","coefficient #0",100,-100000,100000);
 //RooRealVar a1("a1","coefficient #1",-10,-100000,100000);
 //RooRealVar a2("a2","coefficient #2",0,-100000,100000);
   //RooRealVar a3("a3","coefficient #2",0,-100000,100000);
 //RooPolynomial ground("ground","ground",x,RooArgList(a0,a1,a2));
     
   //RooRealVar sigma2("sigma2","width of gaussian",0.01,0.008,0.02);
   RooRealVar alpha1("alpha1","#alpha",1.7,1,5);
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
   RooRealVar meanb("meanb","mean of gaussian",beame,beame-0.003,beame+0.003);
   RooRealVar sigmab("sigmab","width of gaussian",0.005,0.004,0.02);
 //  RooGaussian gausb("gausb","gauss(x,m,s)",x,meanb,sigmab);
   
   RooRealVar alphab("alphab","#alpha",1.7,1,5);
   RooRealVar nnnb("n1","n",100,1,200);
   RooCBShape cbshapeb("cbshapeb","crystal ball",x,meanb,sigmab,alphab,nnnb);
 
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
   sum->fitTo(*dataset,Range(peakvalue-0.1,beame+0.01));
   //gausb.fitTo(*dataset,Range(beame-0.01, beame+0.01));
   dataset->plotOn(xframe);
   //sum->plotOn(xframe,Components(cbshape1),LineStyle(2),LineColor(2));
   //sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3));
   //sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(cbshape),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(cbshapeb),LineStyle(2),LineColor(3));
   //sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(3));
   //gausb.plotOn(xframe,LineStyle(2),LineColor(3));
   sum->plotOn(xframe);
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
   fitpar<<" ene = "<< beame*2 <<"\t signal mean = "<< mean.getVal() << "sigbal sigma = "<< sigma.getVal();
   fitpar<<" background mean = " << meanb.getVal() << " background sigma = " << sigmab.getVal()<<std::endl;

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
   sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
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

 //ofstream outf("f6pi",std::ios::app);
 //outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   
   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   return signal.getVal();
}


double GetEnergy(int run)
{
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

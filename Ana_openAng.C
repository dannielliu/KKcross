//##############################
//### this code is modified ####
//### dliu13@mail.ustc.edu.cn###
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
const char* getPureName(const char* name);
const char* getPureName2(const char* name);
double getEne(const char* name);


int main(int argc, char** argv)
{
  //if (argc==1) Ana(argv[1]);
  //if (argc>1)  Ana(argv[1],argv[2]);
  //else Ana();
  TFile *fileout = new TFile("output/output_openang2.root","update");
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
  //for (int i=0;i<22;i++) thecut[i] = 179;

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


//std::cout<<"output file name is "<<name<<std::endl;
  //TFile *fileout = new TFile(name,"recreate");
  const char *pureName = getPureName2(filename);
  std::cout<<"Pure Name: "<< pureName <<std::endl;
  char name1[1000];
  //sprintf(name1,"output/%s.root",pureName);
//sprintf(name1,"output/%s.root",pureName);
//TFile *dir = new TFile(name1,"recreate");
//TTree *vars = new TTree("vars","vars");

  fileout->cd();
  TDirectory *dir = fileout->GetDirectory(pureName);
  if (dir==0) dir = fileout->mkdir(pureName);;
  dir->cd();
  TTree *vars = new TTree("vars","vars");

  TLorentzVector kap,kam;
  double mass;

  double totp, costheta1, costheta2;
  double costheta, theta;
  double p1,p2;
  double pt1,pt2;
  double tof11=-1000,tof21=-1000;
  int tofid;
  int isrtag=0;
  int tof1tag=0;
  int tof2tag=0;
  int emctag=0;
  double dtof;
  double eop1,eop2;
//vars->Branch("run",&run,"run/I");
//vars->Branch("indexmc",&idxmc,"indexmc/I");
//vars->Branch("pdgid",pdgid,"pdgid[indexmc]/I");
//vars->Branch("motheridx",motheridx,"motheridx[indexmc]/I");
//vars->Branch("costheta",&costheta,"costheta/D");
//vars->Branch("theta",&theta,"theta/D");
  vars->Branch("costheta1",&costheta1,"costheta1/D");
  vars->Branch("costheta2",&costheta2,"costheta2/D");
  vars->Branch("p1",&p1,"p1/D");
  vars->Branch("p2",&p2,"p2/D");
  vars->Branch("pt1",&pt1,"pt1/D");
  vars->Branch("pt2",&pt2,"pt2/D");
  vars->Branch("dtof",&dtof,"dtof/D");
  vars->Branch("eop1",&eop1,"eop1/D");
  vars->Branch("eop2",&eop2,"eop2/D");
  
  //fileout->cd();

  int nentries = tree->GetEntries();
  double pexp = sqrt(pow(Ebeam/2,2)-pow(mka,2));
  std::cout<<"Total entries is " << nentries << std::endl;
  std::cout<<"Beam energy is " << Ebeam << std::endl;
  std::cout<<"Expected K momentum is " << pexp << std::endl;
  int count0=0,count1=0,count2=0,count3=0;
  int count[10]={0};
  int tagmc=0;

  for (int ien=0;ien<nentries;ien++){
    tree->GetEntry(ien);
    
    isrtag = 0;
    for (int j=0;j<idxmc;j++){
      if (pdgid[j] == 22 && motheridx[j]==j) { tagmc++; isrtag=1;}
    }
    if (emcstatusShort!=0x3 && emcstatusInt!=0x3) continue;
    if (ntof1<1) continue;
    if (ntof2<1) continue;
    for (int i=0;i<ntof1;i++){
      if (tofl1[i]==1) tof11=tof1[i];
    }
    for (int i=0;i<ntof2;i++){
      if (tofl2[i]==1) tof21=tof2[i];
    }
 
    count0++;

    kap.SetVectMag(TVector3(kappx,kappy,kappz),mka);
    kam.SetVectMag(TVector3(kampx,kampy,kampz),mka);
    
    // angular information in lab coordinate
    costheta1 = kap.CosTheta();
    costheta2 = kam.CosTheta();
    
    // momentum information in cms coordinalte
    pt1 = kap.Perp();
    pt2 = kam.Perp();
    kap.Boost(-0.011,0,0);
    kam.Boost(-0.011,0,0);
    costheta = kap.Vect().Dot(kam.Vect())/(kap.Vect().Mag()*kam.Vect().Mag());
    theta = kap.Vect().Angle(kam.Vect())/TMath::Pi()*180.;
    mass = (kap+kam).M();
    totp = (kap+kam).Rho();
    p1 = kap.Rho();
    p2 = kam.Rho();
    dtof = tof11 - tof21;
    eop1 = epratio1;
    eop2 = epratio2;

    //vars->Fill();

    if (fabs(costheta1)<0.93 && costheta1 < 0.8 ) { count[4]++;
    if (fabs(costheta2)<0.93 && costheta2 > -0.8){ count[5]++;
    if (theta>m_thecut)
    { count[2]++;
    //if (epratio1<m_epcut)
    { count[0]++;
      //if (epratio2<m_epcut)
      { count[1]++;
        //if (fabs(tof11-tof21)<m_tofcut)
        if (p1<m_pcut && p2<m_pcut)
        {
          count3++;
          vars->Fill();
        //if (p2<m_pcut) {/*datasel2->Fill();*/ count2++; // p2 dis
        //  //datasel1->Fill();
        //  vars->Fill();
        //}
        }
      }
    }
    }
    
    }
    }
  }
  dir->Write();

  ofstream cutflow("cutflow2",std::ios::app);
  //ofstream cutflow("cutflow_cmp665andp01",std::ios::app);
  cutflow<<filename<<std::endl;
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
  cutflow<<"p1-exp<3 sigma    :"<<count2<<std::endl;

  return 0;
 
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
  if (run>=27147 && run<=27288) return 3.08;
  if (run>=28241 && run<=28266) return 3.08;
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

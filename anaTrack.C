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
const char* getPureName(const char* name);
const char* getPureName2(const char* name);
double getEne(const char* name);
/*
double getDeffp(double pt)
{
  if (fabs(pt-0.025)<=0.25 )   return    -0.00152233; // 0.0106329 uncertainty
  if (fabs(pt-0.075)<=0.25 )   return    -0.00809036; // 0.0136784
  if (fabs(pt-0.125)<=0.25 )   return    -0.0459078 ; // 0.023104
  if (fabs(pt-0.175)<=0.25 )   return    -0.0153731 ; // 0.0190038
  if (fabs(pt-0.225)<=0.25 )   return    0.00989558 ; // 0.0149169
  if (fabs(pt-0.275)<=0.25 )   return    0.011055   ; // 0.0112094
  if (fabs(pt-0.325)<=0.25 )   return    0.0257912  ; // 0.0125118
  if (fabs(pt-0.375)<=0.25 )   return    0.0130194  ; // 0.0102408
  if (fabs(pt-0.425)<=0.25 )   return    0.018983   ; // 0.00836595
  if (fabs(pt-0.475)<=0.25 )   return    -0.00927574; // 0.0114852
  if (fabs(pt-0.525)<=0.25 )   return    0.0137483  ; // 0.011598
  if (fabs(pt-0.575)<=0.25 )   return    0.0246526  ; // 0.0163572
  if (fabs(pt-0.625)<=0.25 )   return    -0.00092793; // 0.0183032
  if (fabs(pt-0.675)<=0.25 )   return    -0.00672899; // 0.0171582
  if (fabs(pt-0.725)<=0.25 )   return    -0.00935493; // 0.0129323
  if (fabs(pt-0.775)<=0.25 )   return    0.00330426 ; // 0.0141946
  if (fabs(pt-0.825)<=0.25 )   return    -0.00998153; // 0.0127833
  if (fabs(pt-0.875)<=0.25 )   return    -0.0412788 ; // 0.00988909
  if (fabs(pt-0.925)<=0.25 )   return    0.00461584 ; // 0.0215712
  if (fabs(pt-0.975)<=0.25 )   return    0.00966137 ; //  0.0293682
  if (fabs(pt-1.025)<=0.25 )   return    0.00516112 ; //  0.0226896
  if (fabs(pt-1.075)<=0.25 )   return    0.0148174  ; //  0.020542
  if (fabs(pt-1.125)<=0.25 )   return    -0.0001399 ; //  0.0147543
  if (fabs(pt-1.175)<=0.25 )   return    5.2645e-05 ; //  0.038999
  if (fabs(pt-1.25 )<=0.5  )   return    -0.0746708 ; //  0.0318817
  return 0.0006; // 0.0028
}
double getDeffm(double pt)
{
  if (fabs(pt-0.025)<=0.25 )   return    0          ; //  0.0132117
  if (fabs(pt-0.075)<=0.25 )   return    0.00015375 ; //  0.0126108
  if (fabs(pt-0.125)<=0.25 )   return    -0.0205695 ; //  0.0184559
  if (fabs(pt-0.175)<=0.25 )   return    0.0342148  ; //  0.02635
  if (fabs(pt-0.225)<=0.25 )   return    -0.0127009 ; //  0.0181818
  if (fabs(pt-0.275)<=0.25 )   return    0.00775857 ; //  0.0126351
  if (fabs(pt-0.325)<=0.25 )   return    0.0118153  ; //  0.0116439
  if (fabs(pt-0.375)<=0.25 )   return    0.00931635 ; //  0.0097009
  if (fabs(pt-0.425)<=0.25 )   return    0.0044542  ; //  0.0115362
  if (fabs(pt-0.475)<=0.25 )   return    -0.00624288; //  0.0115809
  if (fabs(pt-0.525)<=0.25 )   return    0.0191736  ; //  0.0163926
  if (fabs(pt-0.575)<=0.25 )   return    0.0144993  ; //  0.00876636
  if (fabs(pt-0.625)<=0.25 )   return    0.00392102 ; //  0.00836115
  if (fabs(pt-0.675)<=0.25 )   return    -0.00648427; //  0.00977012
  if (fabs(pt-0.725)<=0.25 )   return    -0.0286454 ; //  0.00771639
  if (fabs(pt-0.775)<=0.25 )   return    -0.00888085; //  0.0099398
  if (fabs(pt-0.825)<=0.25 )   return    0.0124098  ; //  0.0102919
  if (fabs(pt-0.875)<=0.25 )   return    -0.0285338 ; //  0.0122111
  if (fabs(pt-0.925)<=0.25 )   return    -0.0513527 ; //  0.0143669
  if (fabs(pt-0.975)<=0.25 )   return    -0.008623  ; //  0.0168785
  if (fabs(pt-1.025)<=0.25 )   return    -0.0229544 ; //  0.0140932
  if (fabs(pt-1.075)<=0.25 )   return    -0.0333059 ; //  0.0156714
  if (fabs(pt-1.125)<=0.25 )   return    0.00702537 ; //  0.0178509
  if (fabs(pt-1.175)<=0.25 )   return    -0.0402473 ; //  0.0373916
  if (fabs(pt-1.25 )<=0.5  )   return    -0.0130793 ; //  0.0285311

  return -0.0047; // 0.0025
}
*/

double getDeffp(double pt)
{
  if (fabs(pt-0.05)<=0.05 )   return  -0.114615	  ; // 1.05382   uncertainty
  if (fabs(pt-0.15)<=0.05 )   return  0.0146269	  ; // 0.0459368
  if (fabs(pt-0.25)<=0.05 )   return  0.00801038  ; // 0.0156752
  if (fabs(pt-0.35)<=0.05 )   return  0.00709051  ; // 0.00823205
  if (fabs(pt-0.45)<=0.05 )   return  0.00026086  ; // 0.00899845
  if (fabs(pt-0.55)<=0.05 )   return  -0.0107742  ; // 0.00620605
  if (fabs(pt-0.65)<=0.05 )   return  -0.0148236  ; // 0.0101031
  if (fabs(pt-0.75)<=0.05 )   return  -0.0145307  ; // 0.00725206
  if (fabs(pt-0.85)<=0.05 )   return  -0.0148492  ; // 0.00571412
  if (fabs(pt-0.95)<=0.05 )   return  -0.0083579  ; // 0.00495995
  if (fabs(pt-1.05)<=0.05 )   return  -0.0160191  ; // 0.00908389
  if (fabs(pt-1.15)<=0.05 )   return  -0.0452576  ; // 0.0186137
  if (fabs(pt-1.25)<=0.05 )   return  -0.0151194  ; // 0.0223524
  return -0.0101372;// 0.00239245
}
double getDeffpErr(double pt)
{
  if (fabs(pt-0.05)<=0.05 )   return  1.05382    ; // 1.05382   uncertainty
  if (fabs(pt-0.15)<=0.05 )   return  0.0459368  ; // 0.0459368
  if (fabs(pt-0.25)<=0.05 )   return  0.0156752  ; // 0.0156752
  if (fabs(pt-0.35)<=0.05 )   return  0.00823205 ; // 0.00823205
  if (fabs(pt-0.45)<=0.05 )   return  0.00899845 ; // 0.00899845
  if (fabs(pt-0.55)<=0.05 )   return  0.00620605 ; // 0.00620605
  if (fabs(pt-0.65)<=0.05 )   return  0.0101031  ; // 0.0101031
  if (fabs(pt-0.75)<=0.05 )   return  0.00725206 ; // 0.00725206
  if (fabs(pt-0.85)<=0.05 )   return  0.00571412 ; // 0.00571412
  if (fabs(pt-0.95)<=0.05 )   return  0.00495995 ; // 0.00495995
  if (fabs(pt-1.05)<=0.05 )   return  0.00908389 ; // 0.00908389
  if (fabs(pt-1.15)<=0.05 )   return  0.0186137  ; // 0.0186137
  if (fabs(pt-1.25)<=0.05 )   return  0.0223524  ; // 0.0223524
  return 0.00239245;
}

double getDeffm(double pt)
{
  if (fabs(pt-0.05)<=0.05 )   return    -0.0307692	 ;//0.463593  //  
  if (fabs(pt-0.15)<=0.05 )   return    -0.0580717	 ;//0.056870  //  
  if (fabs(pt-0.25)<=0.05 )   return    -0.0118475	 ;//0.010829  //  
  if (fabs(pt-0.35)<=0.05 )   return    -0.00388159	 ;//0.009016  //  
  if (fabs(pt-0.45)<=0.05 )   return    -0.00263652	 ;//0.018656  //  
  if (fabs(pt-0.55)<=0.05 )   return    -0.00821807	 ;//0.006484  //  
  if (fabs(pt-0.65)<=0.05 )   return    -0.00346634	 ;//0.008547  //  
  if (fabs(pt-0.75)<=0.05 )   return    -0.00167316	 ;//0.005131  //  
  if (fabs(pt-0.85)<=0.05 )   return    0.00740485	 ;//0.010778  //  
  if (fabs(pt-0.95)<=0.05 )   return    -0.021221	 ;//0.008313  //  
  if (fabs(pt-1.05)<=0.05 )   return    -0.0137069	 ;//0.011612  //  
  if (fabs(pt-1.15)<=0.05 )   return    0.018318	 ;//0.021540  //  
  if (fabs(pt-1.25)<=0.05 )   return    0.00197621	 ;//0.041395  //  
  return  -0.00658; //  0.00274522
}
double getDeffmErr(double pt)
{
  if (fabs(pt-0.05)<=0.05 )   return   0.463593  ; //  
  if (fabs(pt-0.15)<=0.05 )   return   0.056870  ; //  
  if (fabs(pt-0.25)<=0.05 )   return   0.010829  ; //  
  if (fabs(pt-0.35)<=0.05 )   return   0.009016  ; //  
  if (fabs(pt-0.45)<=0.05 )   return   0.018656  ; //  
  if (fabs(pt-0.55)<=0.05 )   return   0.006484  ; //  
  if (fabs(pt-0.65)<=0.05 )   return   0.008547  ; //  
  if (fabs(pt-0.75)<=0.05 )   return   0.005131  ; //  
  if (fabs(pt-0.85)<=0.05 )   return   0.010778  ; //  
  if (fabs(pt-0.95)<=0.05 )   return   0.008313  ; //  
  if (fabs(pt-1.05)<=0.05 )   return   0.011612  ; //  
  if (fabs(pt-1.15)<=0.05 )   return   0.021540  ; //  
  if (fabs(pt-1.25)<=0.05 )   return   0.041395  ; //  
  return   0.00274522;
}



int main(int argc, char** argv)
{
  //if (argc==1) Ana(argv[1]);
  //if (argc>1)  Ana(argv[1],argv[2]);
  //else Ana();
  TFile *fileout = new TFile("output/output_cutuncertainty.root","update");
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

  //TFile *fileout = new TFile(name,"recreate");
  const char *pureName = getPureName2(filename);
  std::cout<<"Pure Name: "<< pureName <<std::endl;
  char name1[1000];
  //sprintf(name1,"output/%s.root",pureName);
  sprintf(name1,"output/%s.root",pureName);
  
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
  
  fileout->cd();

  int nentries = tree->GetEntries();
  double pexp = sqrt(pow(Ebeam/2,2)-pow(mka,2));
  std::cout<<"Total entries is " << nentries << std::endl;
  std::cout<<"Beam energy is " << Ebeam << std::endl;
  std::cout<<"Expected K momentum is " << pexp << std::endl;
  int count0=0,count1=0,count2=0,count3=0;
  int count[10]={0};
  int tagmc=0;

  double sumdeffp=0;
  double sumdeffm=0;
  double sumcntp=0;
  double sumcntm=0;

  for (int ien=0;ien<nentries;ien++){
    tree->GetEntry(ien);
    
     //if ( run == 39582||run == 39671||run == 39723||run == 39783||run == 40208||run == 40300||run == 40308||run == 40462||run == 40526||run == 40946||run == 41099||run == 41200||run == 41283||run == 41408||run == 41416||run == 41436||run == 41445||run == 41471||run == 41728||run == 41818||run == 41902) continue;
     
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
    //if (tof11<-999 || tof21<-999) continue;
    count0++;

    kap.SetVectMag(TVector3(kappx,kappy,kappz),mka);
    kam.SetVectMag(TVector3(kampx,kampy,kampz),mka);
    
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
 
    // select candidate, p1 dis
    // if (count0%100000==0) std::cout<<"epratio1 is "<< epratio1<<std::endl;
    if (fabs(costheta1)<0.93 && costheta1 < 0.8) { count[4]++;
    if (fabs(costheta2)<0.93 && costheta2 > -0.8){ count[5]++;
    if (theta>m_thecut){ count[2]++;
    if (epratio1<m_epcut)
    { count[0]++;
      if (epratio2<m_epcut)
      { count[1]++;
        if (fabs(tof11-tof21)<m_tofcut){
          count3++;
          if (p2<m_pcut && p1<m_pcut) {/*datasel2->Fill();*/ count2++; // p2 dis
            //datasel1->Fill();
	sumdeffp += getDeffp(kap.Pt())/pow(getDeffpErr(kap.Pt()),2);
	sumdeffm += getDeffm(kam.Pt())/pow(getDeffmErr(kam.Pt()),2);
	sumcntp += 1.0/pow(getDeffpErr(kap.Pt()),2);
	sumcntm += 1.0/pow(getDeffmErr(kam.Pt()),2);
          }
        }
      }
    }
    }
    
    }
    }
  }
  double avedeffp = sumdeffp/sumcntp;
  double avedeffm = sumdeffm/sumcntm;
  double deff = sqrt(avedeffp*avedeffp+avedeffm*avedeffm);
  cout<< "for file: "<< filename<<endl;
  cout<<"+ tracking error: "<< avedeffp<<endl;
  cout<<"- tracking error: "<< avedeffm<<endl;
  cout<<"tracking error: sqrt    "<< deff<<endl;
  cout<<"tracking error: e+ + e- "<< avedeffp+avedeffm<<endl;
  cout<<"tracking error: abs+ + abs- "<<  fabs(avedeffp)+fabs(avedeffm)<<endl;
  

  ofstream cutflow("cutflow2",std::ios::app);
  cutflow<<filename<<std::endl;
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

  // finish program here
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

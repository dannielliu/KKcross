#include <sstream>
#include <fstream>
#include <iostream>

using namespace std;

//int main(int argc, char** argv)
int cmp_ps_polars()
{
        ifstream infile("cmp_ps_polars.dat");
        char line[1000];
        infile.getline(line,1000);
        int n=0;
        double enes[21];
        double eenes[21];
        double n_p[21];
        double ne_p[21];
        double n_pol[21];
        double ne_pol[21];

        double n_rela[21];
        double ne_rela[21];
         while (!infile.eof()){
	    infile.getline(line,1000);
	    if (line[0]=='\0') continue;
	    istringstream iss;
	    iss.clear();
	    iss.str(line);
	    //iss.seekg(30,ios::beg);
	    iss>>enes[n]>>n_p[n]>>ne_p[n]>>n_pol[n]>>ne_pol[n];
	  //n_pol[n] = n_pol[n];
	  //ne_pol[n] = ne_pol[n];
	    eenes[n] = 0;

	    // relative difference
	    n_rela[n]  = n_pol[n]/n_p[n];
	    ne_rela[n] = n_pol[n]/n_p[n]*sqrt(pow(ne_pol[n]/n_pol[n],2)+pow(ne_p[n]/n_p[n],2));
	    
	    cout<<enes[n]<<"\t"<<(n_rela[n]-1)*100<<endl;
	    
	    n++;
        }

        TGraphErrors* g_p = new TGraphErrors(n,enes,n_p,eenes,ne_p);
        TGraphErrors* g_pol = new TGraphErrors(n,enes,n_pol,eenes,ne_pol);
        g_p->SetMarkerStyle(20);
        g_p->SetMarkerColor(2);
        g_p->SetFillColor(0);
        g_p->SetTitle("");
        g_p->SetLineWidth(2);
        g_p->GetXaxis()->SetNdivisions(505);
        g_p->GetYaxis()->SetNdivisions(505);
        g_p->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
        g_p->GetYaxis()->SetTitle("events");
        g_p->GetXaxis()->SetLabelSize(0.06);
        g_p->GetXaxis()->SetTitleSize(0.06);
        g_p->GetYaxis()->SetLabelSize(0.06);
        g_p->GetYaxis()->SetTitleSize(0.07);
        g_p->GetXaxis()->SetLabelOffset(1.5);
        g_p->GetYaxis()->SetTitleOffset(0.75);
        g_pol->SetMarkerStyle(21);
        g_pol->SetMarkerColor(3);
        g_pol->SetFillColor(0);
        g_pol->SetLineWidth(2);
        
        TLegend *legendf = new TLegend(0.42,0.65,0.95,0.85);
        legendf->SetBorderSize(0); 
        legendf->SetMargin(0.25); 
        legendf->SetFillStyle(0); 
        legendf->AddEntry(g_p,"N from p spectrum");
        legendf->AddEntry(g_pol,"N from polar angle spectrum");

        TCanvas *c1 = new TCanvas();
        TPad *pad1 = new TPad("pad1","pad1",0.05,0.3,1.0,0.98);
        TPad *pad2 = new TPad("pad2","pad2",0.05,0.02,1.0,0.3);
        //c1->SetMargin(0.15,0.1,0.15,0.1);
        pad1->SetTopMargin(0.05);
        pad1->SetBottomMargin(0.02);
        pad2->SetTopMargin(0.05);
        pad2->SetBottomMargin(0.3);
        pad1->Draw();
        pad2->Draw();
        pad1->cd();
 
        g_p->Draw("AP");
        g_pol->Draw("P");
        legendf->Draw();

        pad2->cd();
        pad2->SetGridy();
        // relative comparison
      //TCanvas *c2 = new TCanvas();
      //c2->SetMargin(0.15,0.1,0.15,0.1);
        TGraphErrors* g_rela = new TGraphErrors(n,enes,n_rela,eenes,ne_rela);
        g_rela->SetMarkerStyle(20);
        g_rela->SetMarkerColor(2);
        g_rela->SetFillColor(0);
        g_rela->SetTitle("");
        g_rela->SetLineWidth(2);
        g_rela->GetXaxis()->SetNdivisions(505);
        g_rela->GetYaxis()->SetNdivisions(502);
        g_rela->GetXaxis()->SetTitle("#sqrt{s} (GeV)");
        g_rela->GetYaxis()->SetTitle("N_{pol}/N_{p}");
        g_rela->GetXaxis()->SetLabelSize(0.15);
        g_rela->GetXaxis()->SetTitleSize(0.15);
        g_rela->GetYaxis()->SetLabelSize(0.15);
        g_rela->GetYaxis()->SetTitleSize(0.15);
        g_rela->GetYaxis()->SetTitleOffset(0.3);
        g_rela->GetYaxis()->SetRangeUser(0.5,1.5);
 
        g_rela->Draw("AP");

        c1->cd();
        c1->Print("cmp_ps_pols.pdf");
        return 0;
}

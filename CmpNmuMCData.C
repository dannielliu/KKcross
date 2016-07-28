int CmpNmuMCData()
{
  double ene[21]    = {2.0,     2.05,    2.1,     2.15,   2.175, 
                       2.2,     2.2324,  2.3094,  2.3864, 2.396,  
	           2.5,     2.6444,  2.6464,  2.7,    2.8,
	           2.9,     2.95,    2.981,   3.0,    3.02,     3.08};
  double enee[21]   = {0.0};
  double Nmu_mc[21] = {73.4868, 24.4878, 106.242, 28.242, 123.445, 
                       171.389, 155.009, 336.687, 412.597, 1314.02,
	           26.2072, 1117.27, 1111.97, 37.838, 41.3881,
	           5159.19, 899.304, 924.003, 928.963, 1074.39, 8685.27};
  double Nmue_mc[21]; for (int i=0; i<21;i++) Nmue_mc[i] = sqrt(Nmu_mc[i]);
  double Nmu_da[21] = {91.9,    39.1,    143.0,   25.8,    129.5,
                       163.7,   159.6,   361.1,   443.5,   1240.1,
	           22.0,    1081.5,  1073.7,  53.1,    41.9,
	           5441.4,  894.8,   974.4,   977.4,   1172.0,  8805.7};
  double Nmue_da[21]= {11.9,    6.4,     13.5,    5.,      12.7,
                       13.2,    17.1,    24.9,    24.8,    50.5,
	           4.1,     35.2,    42.7,    7.5,     7.7,
	           124.6,   32.5,    32.4,    31.6,    35.6,    102.7};

  double Nmu_mc_s[21];
  double Nmue_mc_s[21];
  double Nmu_da_s[21];
  double Nmue_da_s[21];
  for (int i=0; i<21; i++){
    Nmu_mc_s[i] = Nmu_mc[i]/Nmu_da[i];
    Nmue_mc_s[i] = Nmue_mc[i]/Nmu_da[i];
    Nmu_da_s[i] = Nmu_da[i]/Nmu_da[i];
    Nmue_da_s[i] = Nmue_da[i]/Nmu_da[i];
    cout << "At "<< ene[i]<<" GeV, Nmu from data is "<< Nmu_da[i]<<", Nmu from MC is "<< Nmu_mc[i]<<endl;
  //cout << ene[i]<<"\t" << Nmu_da[i];
  //cout << "\t" << Nmu_mc[i] << " "<< Nmu_mc_s[i] ;
  //cout << "\t "<< Nmue_mc[i] << " "<< Nmue_mc_s[i] ;
  //cout << "\t "<< Nmu_da[i] <<" "<< Nmu_da_s[i];
  //cout <<"\t " << Nmue_da[i] <<" "<<  Nmue_da_s[i] << endl;
  }

  TGraphErrors *gdata = new TGraphErrors(21, ene, Nmu_da_s, enee, Nmue_da_s);
  TGraphErrors *gmcmu = new TGraphErrors(21, ene, Nmu_mc_s, enee, Nmue_mc_s);
  gdata->Draw("AP");
  gdata->SetFillColor(0);
  gmcmu->SetFillColor(0);
  gmcmu->SetLineColor(3);
  gmcmu->SetMarkerColor(3);
  gmcmu->Draw("P");
  TLegend *lg = new TLegend(0.5,0.5,0.9,0.9);
  lg->SetFillStyle(0);
  lg->AddEntry(gdata,"Nmu from data");
  lg->AddEntry(gmcmu,"Nmu from MC");
  lg->Draw();




}

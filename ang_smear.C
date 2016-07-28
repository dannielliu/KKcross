int ang_smear()
{
  ifstream ifile("ang_smear.txt");
  double ene[22];
  double enee[22];
  double mean[22], meane[22],sigma[22],sigmae[22];
  for (int i=0; i<22;i++){
    ifile >> ene[i] >> mean[i]>>meane[i]>>sigma[i]>>sigmae[i];
    ene[i] *=2;
    enee[i] =0;
  }

  //TGraphErrors *graph_m = new TGraphErrors(22,ene,mean,enee,meane);
  TGraphErrors *graph_m = new TGraphErrors(22,ene,sigma,enee,sigmae);
  graph_m->Draw();
}

int tmpisr()
{
  ifstream infile("tmp.txt");
  istringstream iss;
  char line[1000];

  while (!infile.eof()){
    iss.clear();
    infile.getline(line,1000); 
    iss.str(line);
    double energy, mceff, isr, lum, tag, ep, ang, tof, trk;;
    iss>>energy>>mceff>>isr>>lum>>tag>>ep>>ang>>tof>>trk;
    double sum = sqrt(mceff*mceff+isr*isr+lum*lum+tag*tag+ep*ep+ang*ang+tof*tof+trk*trk);
    cout<<energy<<"\t"<<sum<<endl;
  }
//while (!infile.eof()){
//  iss.clear();
//  infile.getline(line,1000); 
//  iss.str(line);
//  double isrnew, isrold;
//  iss>>isrnew>>isrold;
//  cout<<(isrnew-isrold)/isrold*100<<endl;
//}

}

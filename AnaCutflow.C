#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

int main(int argc, char **argv)
{
  char filename[1000];
  if (argc>1) sprintf(filename,"%s",argv[1]);
  ifstream infile(filename);
  char linec[1000];
  ofstream outflow("cut_out");
  while (!infile.eof()){
    infile.getline(linec,1000);
    //if (strncmp(linec,"###",3)==0) break;
    if (linec[0]=='#' || linec[0]=='\0' || linec[0] == '\n') continue;
    if (linec[0]=='/') {
      outflow<<std::endl;
      outflow << linec<<"\t";
    }
    else 
      outflow << &linec[19]<<"\t";

  //infile.getline(linec,1000);
  //outflow << &linec[19]<<"\t";
  //infile.getline(linec,1000);
  //outflow << &linec[19]<<"\t";
  //
  //infile.getline(linec,1000);
  //outflow << &linec[19]<<"\t";
  //infile.getline(linec,1000);
  //outflow << &linec[19]<<"\t";
  //infile.getline(linec,1000);
  //outflow << &linec[19]<<"\t";
  //
  //infile.getline(linec,1000);
  //outflow << &linec[19]<<"\t";
  //infile.getline(linec,1000);
  //outflow << &linec[19]<<"\t";
  //infile.getline(linec,1000);
  //outflow << &linec[19]<<"\t";

  }
}



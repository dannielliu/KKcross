int PrintPSpec()
{
  TFile *file = new TFile("output_mc.root");
  int n = file->GetListOfKeys()->GetSize();
  for (int i=0; i<n; i++){
    char *name;
    TCanvas *c1;
    name = file->GetListOfKeys()->At(i)->GetName();
    if (name[0] != 'p') continue;
    c1 = (TCanvas*)file->Get(name);
    char outfile[1000];
    sprintf(outfile,"%s.pdf",name);
    c1->Draw();
    c1->Print(outfile);
  }
  return 0;
}

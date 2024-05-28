 
void cla::upclass(cla A){

  TFile* output = new TFile("myclass.root", "recreate");

  TTree* const_tree = new TTree("const_tree", "Tree with constants");
  TBranch* b = const_tree->Branch("dabranch", &A);

  const_tree->Fill();

  output->WriteObject(const_tree, "const_tree");
  output->Close();
}  

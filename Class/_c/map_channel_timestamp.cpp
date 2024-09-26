
void cla::TooAnnoying(){
  int len = 1024;
  int wf_counter = 0;

  std::map<Short_t,std::vector<ULong64_t>> m;

  wffunctions bs;

  map<TString, TString> filename = {{"file1", wf_file}};

  for (auto f : filename){

     TChain *t[] = {NULL};

     t[0] = new TChain();
     t[0]->Add(Form("%s?#raw_waveforms", f.second.Data()));
     // t[0]->SetImplicitMT(true);
     Long64_t nentries = t[0]->GetEntries();
     hdf5torootclass event(t[0]);

     cout << "\nFile open -> " << f.second << "\tentries: " << nentries << endl;

     for (size_t ievt=0; ievt<nentries && wf_counter<n_wf; ievt++){ // loop over entries in root file

      event.GetEntry(ievt);

      m[event.channel].push_back(event.timestamp);
      wf_counter++;
     }
  
    std::map<Short_t,std::vector<ULong64_t>>::iterator it;
    for (it = m.begin(); it != m.end(); it++){
      std::sort(it->second.begin(), it->second.end());
      std::cout << it->first << " ";
      for(size_t i=0; i<5; i++) std::cout << it->second[i] << " ";
      std::cout << " " << std::endl;

    }

  }

}

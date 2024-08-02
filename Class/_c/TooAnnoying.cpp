
void cla::TooAnnoying(){
  
  CompleteWF_Binary_Swap(wf_file, wfs, n_wf, memorydepth); 
  n_wf = wfs.size();
  vector<vector<double>> wfs2 = wfs;

  SubBaseline(wfs, prepulse_ticks, true);
  SubBaseline2(wfs2, 3., true);

  DisplayWFs (wfs2, wfs, 1., 20);
}

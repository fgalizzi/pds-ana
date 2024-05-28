void cla::read(){
  if(wfs.size()!=n_wf || oldwf_file!=wf_file || oldprepulse_ticks != prepulse_ticks){
    oldwf_file = wf_file;
    oldprepulse_ticks = prepulse_ticks;
    wfs.erase(wfs.begin(),wfs.end());
    if(data_format == "caen")   CAEN_WF_Binary(wf_file, wfs, n_wf, memorydepth);
    if(data_format == "daphne") CompleteWF_Binary(wf_file, wfs, n_wf, memorydepth);
    if(data_format == "esteban") CompleteWF_Binary_Swap(wf_file, wfs, n_wf, memorydepth);
    if(data_format == "csv")    CSV_WF_Binary(wf_file, wfs, n_wf, memorydepth);
    
    if(invert == false) SubBaseline(wfs, prepulse_ticks);
    if(invert == true ) SubBaseline_Invert(wfs, prepulse_ticks);
  }

}

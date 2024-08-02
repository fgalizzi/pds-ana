void cla::read(){
  if(wfs.size()!=n_wf || oldwf_file!=wf_file || oldprepulse_ticks != prepulse_ticks ||
      oldchannel!=channel){
    //Store old values to decide whether to re-read the file 
    oldwf_file = wf_file;
    oldprepulse_ticks = prepulse_ticks;
    oldchannel = channel;

    wfs.erase(wfs.begin(),wfs.end());
    
    if(data_format == "caen")    CAEN_WF_Binary(wf_file, wfs, n_wf, memorydepth);
    if(data_format == "daphne")  CompleteWF_Binary(wf_file, wfs, n_wf, memorydepth);
    if(data_format == "esteban") CompleteWF_Binary_Swap(wf_file, wfs, n_wf, memorydepth);
    if(data_format == "csv")     CSV_WF_Binary(wf_file, wfs, n_wf, memorydepth);
    
    //ProtoDUNE-HD: update n_wf because the reading function stops automatically
    if(data_format == "pdhd"){
      PDHD_ch_wfs(wf_file, wfs, channel, n_wf); 
      n_wf = wfs.size();
    }
   
    //Subtract the baseline and invert the wfs according to "invert"
    if(sub_bsl == true) SubBaseline(wfs, prepulse_ticks, invert);
  }

}

void cla::read(){
  if(WFS.size()!=MEMORYDEPTH*N_WF){
    WFS.erase(WFS.begin(),WFS.end());
    if(DATA_FORMAT == "caen")   CAEN_WF_Binary(WF_FILE, WFS, N_WF, MEMORYDEPTH);
    if(DATA_FORMAT == "daphne") CompleteWF_Binary(WF_FILE, WFS, N_WF, MEMORYDEPTH);
    if(DATA_FORMAT == "csv")    CSV_WF_Binary(WF_FILE, WFS, N_WF, MEMORYDEPTH);
    
    if(INVERT == false) SubBaseline(WFS, MEMORYDEPTH, PREPULSE_TICKS);
    if(INVERT == true ) SubBaseline_Invert(WFS, MEMORYDEPTH, PREPULSE_TICKS);
  }

}

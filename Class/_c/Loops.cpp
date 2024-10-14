#include "../classe.hpp"

void cla::Loop_ST_Analysis(){
  //////////////////////////////////////////////////////////////////////
  // HARD CODE HERE
  //
  ep = 112;
  mask = 1400;
  //////////////////////////////////////////////////////////////////////

  vector<string> signal_files = read_chs("ch_sipm.txt");
  vector<string> self_files   = read_chs("ch_self.txt");
  vector<size_t> channels = read_pdhd_ch_map(mask);
  print = 1;
  
  for(size_t i=0; i<signal_files.size(); i++){
    wf_file = signal_files[i];
    trg_f = self_files[i];
    //Check wether the LED mask was good for this channel
    size_t this_channel = (size_t)extract_channel_from_filename(wf_file);
    int cnt = std::count(channels.begin(), channels.end(), ep*100+this_channel);
    if (cnt == 0) continue;
    
    //Calibrate and perform the self-trigger analisys
    LED_Analysis();
    LoadFitParameters(fgaus);
    ST_Analysis();
  }

}

#include "../classe.hpp"

void cla::Loop_ST_Analysis(){
  //////////////////////////////////////////////////////////////////////
  // HARD CODE HERE
  //
  calibration_run = 9999; //Needed create the root file with the results
  ep = 112;
  mask = 1400;
  channel_low = ep*100;     //Lower channel to look at (included)
  channel_up  = (ep+1)*100;     //Upper " "
  //////////////////////////////////////////////////////////////////////

  vector<string> signal_files = read_chs("ch_sipm.txt");
  vector<string> self_files   = read_chs("ch_self.txt");
  vector<size_t> channels = read_pdhd_ch_map(mask);
  print = 1;
  
  for(size_t i=0; i<signal_files.size(); i++){
    wf_file = signal_files[i];
    trg_f = self_files[i];
    //Check wether the LED mask was good for this channel
    channel = size_t(extract_channel_from_filename(wf_file))+ep*100;
    int cnt = std::count(channels.begin(), channels.end(), channel);
    std::cout << cnt << std::endl;
    if (cnt == 0 || channel < channel_low || channel > channel_up) continue;
    std::cout << i << std::endl;
    
    //Calibrate and perform the self-trigger analisys
    LED_Analysis();
    LoadFitParameters(fgaus);
    ST_Analysis();
  }

}

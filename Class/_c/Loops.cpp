#include "../classe.hpp"
#include <cstddef>
#include <utility>

// ****************************************************************
// Description
// ****************************************************************

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
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

void cla::Loop_RMS_Analysis(){
  size_t run = 27877; //Needed create the root file with the results
  channel_low = 10400;     //Lower channel to look at (included)
  channel_up  = 11050;     //Upper " "
  vector<size_t> channels = read_pdhd_ch_map(mask);
  print = 1;
  int nwf_here = 20000;
  prepulse_ticks = 1023;
  vector<vector<double>> y2;
  vector<pair<size_t,double>> ch_rms;


  for(auto ch : channels){
    if(ch<channel_low || ch>channel_up) continue;

    this->channel = ch;
    n_wf = nwf_here;
    read();
    SelCalib_WF(wfs, y2, prepulse_ticks, sat_low, sat_up, bsl);


    TH2D* h2 = new TH2D("h2", Form("%s;%s;%s", "Persistence", "Ticks", "ADC Counts"),
                      memorydepth, 0., memorydepth, int(bsl*4), -bsl*2, bsl*2);
  
    for (auto& wf : y2) for (int j=0; j<memorydepth; j=j+2) h2->Fill(j, wf[j]);
  

    TH1D* h_bsl = h2->ProjectionY("h_bsl", 0, prepulse_ticks);
    ch_rms.push_back(make_pair(channel, h_bsl->GetRMS()));
    h_bsl->Delete();
    h2->Delete();

  }
 
    std::ofstream outfile(Form("Ch_rms_bsl_%f_run_%zu", bsl, run));

    for (const auto& pair : ch_rms) {
        outfile << pair.first << "\t" << pair.second << std::endl;
    }
    outfile.close();

}












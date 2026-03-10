#include "../classe.hpp"

// ****************************************************************
// Change this marco when you can easily loop over runs and you
// don't want to repeat the same commands every time!
// If the loop could be useful for the future, store it in Loops.cpp
// ****************************************************************

///////////////////////////////////////////////////////////////////
//////// HARD CODE ////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////
std::vector<std::pair<int, std::vector<int>>> vgains_offsets = {
  {1500, {2450}}
};

std::vector<int> channels_to_analyze = {4,5,6,7,12,13,14,15,18,19,21};

std::string base_folder = "/Users/federico/CERN/M1/cb_mar_26/setting_up/";

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void cla::TooAnnoying(){
  // set batch mode to avoid opening canvases
  gROOT->SetBatch(kTRUE);

  memorydepth = 1024;
  data_format = "esteban";
  sub_bsl = false;
  peak_to_peak = 2;

  double y0 = 0;
  double y1 = 16500;

  for(auto & vgain_offset : vgains_offsets){
    int vgain = vgain_offset.first;
    std::vector<int> offsets = vgain_offset.second;

    for(auto & offset : offsets){

    for(size_t i=0; i<channels_to_analyze.size(); i++){
      n_wf = -1;
      int channel = channels_to_analyze[i];
      wf_file = base_folder+"vgain_"+std::to_string(vgain)+"_offset_"+std::to_string(offset)+"/raw/channel_"+std::to_string(channel)+".dat";
      std::cout << wf_file << std::endl;
      read();
      std::cout << vgain << " " << offset << " " << channel << std::endl;

      TCanvas* cTime = new TCanvas(Form("wavedec_ch%i_vgain%i_offset%i", channel, vgain, offset), Form("wavedec_ch%i_vgain%i_offset%i", channel, vgain, offset));
      TH2D* h2 = new TH2D("h2", Form("%s;%s;%s", "DecoWF", "Ticks", "ADC Counts"), memorydepth/2, 0., memorydepth, 200, y0, y1); 
      for (auto  &wf : wfs) for (int j=0; j<memorydepth; j=j+2) h2->Fill(j, wf[j]);
      h2->Draw("COLZ");
      gPad->SetLogz();
        cTime->SaveAs(Form("/Users/federico/CERN/M1/cb_mar_26/setting_up/plots/vgain_%i_offset_%i_channel_%i.png", vgain, offsets[0], channel));
        delete cTime;
      }
    }
  }
   
  return;   
}

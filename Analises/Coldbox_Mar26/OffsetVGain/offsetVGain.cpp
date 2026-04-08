#include "../../../Class/_c/class_include.hpp"
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"

std::vector<std::pair<int, std::vector<int>>> vgains_offsets = {
  {500, {2000, 2050, 2100, 2150, 2200, 2250, 2300, 2400, 2500, 2550}},
  {750, {2000, 2050, 2100, 2150, 2200, 2250, 2300, 2400, 2500, 2550}},
  {1000, {2000, 2050, 2100, 2150, 2200, 2250, 2300, 2400, 2500, 2550}},
  {1250, {2000, 2050, 2100, 2150, 2200, 2250, 2300, 2400, 2500, 2550}},
  {1500, {2000, 2050, 2100, 2150, 2200, 2250, 2300, 2400, 2500, 2550}},
  {1750, {2000, 2050, 2100, 2150, 2200, 2250, 2300, 2400, 2500, 2550}},
  {2000, {2000, 2050, 2100, 2150, 2200, 2250, 2300, 2400, 2500, 2550}},
  {2250, {2000, 2050, 2100, 2150, 2200, 2250, 2300, 2400, 2500, 2550}},
  {2500, {2000, 2050, 2100, 2150, 2200, 2250, 2300, 2400, 2500, 2550}},
  {2750, {2000, 2050, 2100, 2150, 2200, 2250, 2300, 2400, 2500, 2550}},
  {3000, {2000, 2050, 2100, 2150, 2200, 2250, 2300, 2400, 2500, 2550}},
};

std::vector<int> channels_to_analyze = {4,5,6,7,12,13,14,15,18,19,20,21};

std::string base_folder = "/eos/experiment/neutplatform/protodune/experiments/ColdBoxVD/March2026run/spybuffer/offset_vgain_scan_SiPM_on/";

void offsetVGain(){
  gROOT->SetBatch(kTRUE);
  auto a = cla();
  a.memorydepth = 1024;
  a.data_format = "esteban";
  a.sub_bsl = false;
  a.peak_to_peak = 2;

  double y0 = 0;
  double y1 = 16500;

  std::string data_folder = base_folder+"data/";
  std::string ana_folder = base_folder+"analysis/";

  for(auto & vgain_offset : vgains_offsets){
    int vgain = vgain_offset.first;
    // if ana_folder/vgain/  does not exist, create it
    std::string vgain_folder = ana_folder+"vgain_"+std::to_string(vgain);
    if (vgain < 1000) vgain_folder = ana_folder+"vgain_0"+std::to_string(vgain);
    if (!std::filesystem::exists(vgain_folder)) {
      std::filesystem::create_directory(vgain_folder);
    }
    vgain_folder = vgain_folder+"/";

    std::vector<int> offsets = vgain_offset.second;
    for(auto & offset : offsets){
    std::string offset_folder = vgain_folder+"offset_"+std::to_string(offset);
    if (!std::filesystem::exists(offset_folder)) {
      std::filesystem::create_directory(offset_folder);
    }
    offset_folder = offset_folder+"/";

    for(size_t i=0; i<channels_to_analyze.size(); i++){
      a.n_wf = -1;
      int channel = channels_to_analyze[i];
      a.wf_file = data_folder+"vgain_"+std::to_string(vgain)+"_offset_"+std::to_string(offset)+"/raw/channel_"+std::to_string(channel)+".dat";
      if (vgain < 1000) a.wf_file = data_folder+"vgain_0"+std::to_string(vgain)+"_offset_"+std::to_string(offset)+"/raw/channel_"+std::to_string(channel)+".dat";
      std::cout << "Opening " << a.wf_file << std::endl;
      a.read();
      std::cout << vgain << " " << offset << " " << channel << std::endl;

      if (a.wfs.size() == 0){
          std::cout << "No waveforms found for vgain " << vgain << " offset " << offset << " channel " << channel << std::endl;
          continue;
        }

      TCanvas* cTime = new TCanvas(Form("wavedec_ch%i_vgain%i_offset%i", channel, vgain, offset), Form("wavedec_ch%i_vgain%i_offset%i", channel, vgain, offset), 800, 600);
      TH2D* h2 = new TH2D("h2", Form("%s;%s;%s", "WF", "Ticks", "ADC Counts"), a.memorydepth/2, 0., a.memorydepth, 200, y0, y1); 
      for (auto  &wf : a.wfs) for (int j=0; j<a.memorydepth; j=j+2) h2->Fill(j, wf[j]);
      h2->Draw("COLZ");
      gPad->SetLogz();
        cTime->SaveAs(Form("%svgain_%i_offset_%i_channel_%i.png", offset_folder.c_str(), vgain, offset, channel));
        delete cTime;
      }
    }
  }
   
  return;   
}

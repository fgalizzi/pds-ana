#include "../classe.hpp"

// ****************************************************************
// Description
// ****************************************************************

size_t noise_run = 27877;

void cla::Pdhd_FFT(){
  //vector<size_t> channels = {11145, 11147};//, 11145, 11147};
  vector<size_t> channels = read_pdhd_ch_map();
  vector<vector<double>> sel_wf;
  double t;
  int x_point = memorydepth/2+1;
  double y_points[600] = {-70.};

  TFile hf(Form("NoiseRun_%zu.root", noise_run), "recreate"); 
  vector<TGraph*> g_fft;

  // To better select the waferoms !
  prepulse_ticks = memorydepth-1;

  for(size_t ch_index=0; ch_index<channels.size(); ch_index++) {
    
    channel = channels[ch_index];
    if (channel > 11000-1) continue;
    // if (channel < 11100 || channel > 11400) continue;
    // Read the wfs for this channel and subtract the baseline
    n_wf = 10000;
    read();
    n_wf = 10000;

    SelCalib_WF(wfs, sel_wf, prepulse_ticks, -bsl, bsl, bsl);  
    
    TGraph* gr;
    if (sel_wf.size()>5) {
      gr = build_ch_fft(memorydepth, tick_len*memorydepth, tick_len, sel_wf, res);
    }
    else{
      gr = new TGraph(x_point, y_points);
    }

    std::ofstream OutFile  (Form("FFT_ch_%lu.dat", channels[ch_index]), ios::binary);
    for(int i = 0; i < gr->GetN(); i++){
      t = gr->GetPointY(i);
      OutFile.write(reinterpret_cast<char*> (&t), sizeof(t));
      t = gr->GetPointX(i);
      OutFile.write(reinterpret_cast<char*> (&t), sizeof(t));
    }   
    std::cout << "Vector saved in ---> " << std::endl;
    gr->SetTitle(Form("Channel_%lu", channels[ch_index]));
    gr->SetName(Form("Channel_%lu", channels[ch_index]));
    g_fft.push_back(gr);
    OutFile.close();
    
    std::cout << "\n\nsel / wfs = " << sel_wf.size() << "/" << wfs.size() << std::endl;
    sel_wf.erase(sel_wf.begin(), sel_wf.end());
  } 

  // Save in file.root
  for(auto g : g_fft) g->Write();
  
  hf.Close();
}

size_t noise_run = 27875;

void cla::Pdhd_FFT(){
  gStyle->SetOptFit(1111); gStyle->SetOptTitle(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit"); 
  
  vector<size_t> channels = {11145, 11147};//, 11145, 11147};
  //vector<size_t> channels = read_ch_map();
  vector<vector<double>> sel_wf;
  double t;

  TFile hf(Form("NoiseRun_%zu.root", noise_run), "recreate"); 
  vector<TGraph*> g_fft;

  for(size_t ch_index=0; ch_index<channels.size(); ch_index++) {
    
    channel = channels[ch_index];
    if (channel < 11100 || channel > 11400) continue;
    // Read the wfs for this channel and subtract the baseline
    read();

    SelCalib_WF(wfs, sel_wf, prepulse_ticks, -bsl, bsl, bsl);  
    
    TGraph* gr = build_ch_fft(memorydepth, tick_len*memorydepth, tick_len,
        sel_wf, res);

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
   
    sel_wf.erase(sel_wf.begin(), sel_wf.end());
  } 

  // Save in file.root
  for(auto g : g_fft) g->Write();
  
  hf.Close();
}

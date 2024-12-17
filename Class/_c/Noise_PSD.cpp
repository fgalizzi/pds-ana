#include "../classe.hpp"

// ****************************************************************
// Description
// ****************************************************************

///////////////////////////////////////////////////////////////////
//////// HARD CODE ////////////////////////////////////////////////

string ofilename_noise_fft    = "";
string ofilename_avgnoise_fft = "";
///////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void cla::Noise_PSD(){
  double t;
  vector<double> x, avg, int_wf, noise_td;
  vector<vector<double>> noise, noise2, avg_wf;
    
  for (size_t i = 0; i < memorydepth; i++) x.push_back( (double) i);
  string ifile_noisetd          = noise_f;
  string ofilename_noise_td     = trg_f;

  // Read and subtract the baseline
  read();
  // Subtract the coherent noise of the digitiser (only once, if you re-run the macro)
  if (ifile_noisetd!=""){
    CompleteWF_Binary(ifile_noisetd, noise_td, memorydepth);
    if(ite==0){
      SubVec_to_WFs(wfs, noise_td);
      ite++;
    }
  }
  
  SelCalib_WF(wfs, noise, prepulse_ticks, -bsl, bsl, bsl);
  TH1D* hI = BuildRawChargeHisto(noise, int_wf, int_low, int_up, nbins);
  Avg_Sel_WF (noise, noise2, avg, int_wf, mu0_low, mu0_up);

  avg_wf.push_back(avg);
  
  TGraph* gAvg= build_avg_spectral_density(memorydepth, tick_len*memorydepth, tick_len, avg_wf, res);
  gAvg->SetLineColor(2);
  TGraph* gNoise_spectral_density = build_avg_spectral_density(memorydepth,
      tick_len*memorydepth, tick_len, noise2, res);

  gNoise_spectral_density->SaveAs(muon_f.c_str());

  TCanvas *c2 = new TCanvas("c2","c2",20,20,1000,900);
  c2->cd();
  gNoise_spectral_density->Draw();
  //gAvg->Draw();
  c2->Modified();
  c2->Update();
  
  TGraph *g1 = new TGraph(avg.size(), &x[0], &avg[0]);
  g1->GetXaxis()->SetTitle("Time [#mus]");
  g1->GetYaxis()->SetTitle("ADC counts");
  g1->SetLineColor(2);

  TCanvas *c3 = new TCanvas("c3","c3",30,30,1000,900);
  c3->cd();
  g1->Draw();
  c3->Modified();
  c3->Update();

  if(print==true){
    if (ofilename_noise_fft!=""){
      std::ofstream OutFile  (ofilename_noise_fft, ios::binary);
      for(int i = 0; i < gNoise_spectral_density->GetN(); i++){
        // y-axis
        t = gNoise_spectral_density->GetPointY(i);
        OutFile.write(reinterpret_cast<char*> (&t), sizeof(t));
        // x-axis
        t = gNoise_spectral_density->GetPointX(i);
        OutFile.write(reinterpret_cast<char*> (&t), sizeof(t));
        }
      OutFile.close();
    }
    if (ofilename_avgnoise_fft!=""){
      std::ofstream OutFile2 (ofilename_avgnoise_fft, ios::binary);
      for(int i = 0; i < gAvg->GetN(); i++){
        // y-axis
        t = gAvg->GetPointY(i);
        OutFile2.write(reinterpret_cast<char*> (&t), sizeof(t));
        // x-axis
        t = gAvg->GetPointX(i);
        OutFile2.write(reinterpret_cast<char*> (&t), sizeof(t));
      } 
      OutFile2.close();
    }
    if (ofilename_noise_td!="") VecDouble_in_Binary("Noise_files/"+ofilename_noise_td+".dat", avg);
    std::cout << "\n\n Saved the files you requested : )\n\n" << std::endl;
  }
  print=false;
}
 

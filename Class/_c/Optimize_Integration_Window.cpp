#include "../classe.hpp"

void cla::Optimise_Integration_Window(size_t min_int_up, size_t max_int_up, size_t increment){
  vector<pair<int, double>> tuple_intup_snr;
  double best_snr = 0;
  int best_int_up = 0;

  for(size_t try_int_up = min_int_up; try_int_up < max_int_up; try_int_up+=increment){
    int_up = try_int_up;
    LED_Analysis();
    tuple_intup_snr.push_back({int_up, SNR});
    
    if(SNR > best_snr){
      best_snr = SNR;
      best_int_up = int_up;
    }
  }

  TGraph* graph = new TGraph();

  for(auto& tuple : tuple_intup_snr){
    cout << "int_up = " << tuple.first << " SNR = " << tuple.second << endl;
    graph->SetPoint(graph->GetN(), tuple.first, tuple.second);
  }

  cout << "Best int_up and SNR" << endl;
  cout << "int_up = " << best_int_up << " SNR = " << best_snr << endl;

  TCanvas* snr_vs_intup = new TCanvas("snr_vs_intup ","snr_vs_intup ",0,0,800,600);
  snr_vs_intup->cd();
  graph->SetTitle("SNR vs int_up;int_up;SNR");
  graph->Draw();
  snr_vs_intup->Modified(); snr_vs_intup->Update();
}

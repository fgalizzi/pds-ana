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

  for(auto& tuple : tuple_intup_snr){
    cout << "int_up = " << tuple.first << " SNR = " << tuple.second << endl;
  }

  cout << "Best int_up and SNR" << endl;
  cout << "int_up = " << best_int_up << " SNR = " << best_snr << endl;

}

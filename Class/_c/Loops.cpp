#include "../classe.hpp"

void cla::Loop_ST_Analysis(){
  vector<string> signal_files = read_chs("ch_sipm.txt");
  vector<string> self_files   = read_chs("ch_self.txt");
  print = 1;
  
  for(size_t i=0; i<signal_files.size(); i++){
    wf_file = signal_files[i];
    trg_f = self_files[i];
    LED_Analysis();
    LoadFitParameters(fgaus);
    ST_Analysis();
  }

}

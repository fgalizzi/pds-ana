#include "../classe.hpp"
#include <cstdio>
#include <numeric>
#include <string>

// ****************************************************************
// This macro is used to create an average muon waveform and save it
// in a binary file, setting print=true, once you are happy with the
// selection cuts. The idea is to select the muon waveforms on an
// f_prompt basis and remove the ones that have light in coincidence.
//
// Parameters you can tune:
// - int_low, int_prompt, int_up <-> integration limits
// - f_prompt (integral[int_low;int_prompt]/integral[int_low;int_up])
//      We select waveforms with f_prompt lower than this value.
// - sat_low <-> discard wfs with samples below this level
// - amp_low, amp_up <-> the max element in the [int_low;int_prompt]
//      range must fall in [amp_low;amp_up]
// - bsl <-> accept only wfs with sample whithin [-bsl;+bsl] in the
//      pre trigger [0;prepulse_ticks]
// - rms <-> the selected wfs can differ from the average only by random
//      fluctuation, we must ensure there are no pulses in the tail of
//      our signal.
// ****************************************************************

////////////////////////////////////////////////////////
/////// HARD CODE //////////////////////////////////////

std::string muon_files_path = "/Users/federico/PhD/PDE/Muon_files/";


////////////////////////////////////////////////////////

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void cla::Avg_Muon(){
  std::vector<double> avg_mu;
  std::vector<std::vector<double>> mu_wfs;

  read();
  
  // TH1D* h_prompt = AllFpromptHisto(wfs, int_low, int_up, int_prompt);
  TH2D* h2_fprompt_charge_wfs = BuildChargeFpromptHisto(wfs, int_low, int_up, int_prompt);
  TH1D* h_fprompt_wfs = h2_fprompt_charge_wfs->ProjectionY("h_prompt_wfs");

  // SelPDE_WF(wfs, sel_wf, prepulse_ticks, int_prompt,
  //           sat_low, amp_low, amp_up, bsl, rms); 
  

  //-----------------------------------------------------------------
  //----- SELECTION LOOPS -------------------------------------------
  double max_el, min_el, t;
  size_t len = wfs[0].size();
  size_t nwfs = wfs.size();
  int bsl_counter = 0;
  int amp_counter = 0;
  int fpr_counter = 0;
  int sel_counter = 0;

  vector<double> first_avg(len, 0.0);
  vector<bool> preliminary_selection(nwfs, false);
  vector<bool> final_selection(nwfs, false);

  mu_wfs.reserve(nwfs);
    
  for (size_t i=0; i<nwfs; i++) {
    max_el = *max_element(wfs[i].begin(), wfs[i].begin()+prepulse_ticks);
    min_el = *min_element(wfs[i].begin(), wfs[i].begin()+prepulse_ticks);
   
    // Select wfs with no pulses in the pre-trigger
    if (max_el<bsl && min_el > -bsl) {
      bsl_counter++; 
      max_el = *max_element(wfs[i].begin()+prepulse_ticks, wfs[i].begin()+int_prompt);
      min_el = *min_element(wfs[i].begin()+prepulse_ticks, wfs[i].begin()+int_prompt);

      // Select wfs with amplitude in given range and with no lower saturation
      if (max_el<amp_up && max_el>amp_low && min_el > sat_low){
        amp_counter++;
        double integral_prompt = accumulate(wfs[i].begin()+int_low, wfs[i].begin()+int_prompt, 0.);
        double integral        = accumulate(wfs[i].begin()+int_low, wfs[i].begin()+int_up, 0.);
        double fprompt = integral_prompt/integral;
        
        // Select wfs with frompt < f_prompt
        if (fprompt<f_prompt){
          fpr_counter++;
          t = max_el;
          max_el = *max_element(wfs[i].begin()+int_prompt, wfs[i].end());
          min_el = *min_element(wfs[i].begin()+int_prompt, wfs[i].end());

          // Select wfs with no lower-saturation and other big-pulses in the tail
          if (max_el < t && min_el > sat_low){
            sel_counter++;
            for(size_t j=0; j<len; j++) first_avg[j] += wfs[i][j];
            preliminary_selection[i]=true;
          } 
        }     
      }
    }
  }  
 
  double norm = 1./ *max_element(first_avg.begin(), first_avg.end());
  for (auto& e : first_avg) e*=norm;

  for (size_t i=0; i<nwfs; i++) {
    if(preliminary_selection[i]==true){
      vector wf = wfs[i];
      norm = 1./ *max_element(wf.begin(), wf.end());
      double norm_bsl = rms*norm*5.;
      for (auto& e : wf) e *= norm;

      bool select_shape = true;
      for (size_t j=0; j<len; j++){
        if (wf[j] > first_avg[j]+norm_bsl || wf[j] < first_avg[j]-norm_bsl){
          select_shape = false;
          j = len+1;
        }
      }
      final_selection[i]=select_shape;
    }

    if(final_selection[i]==true) mu_wfs.push_back(wfs[i]);
  }
  
  cout << "WFs selected on the baseline "  << bsl_counter << "/" << nwfs << endl;
  cout << "WFs selected on the amplitude " << amp_counter << "/" << nwfs << endl;
  cout << "WFs selected on the fprompt " << fpr_counter << "/" << nwfs << endl;
  cout << "WFs preliminary selected " << sel_counter << "/" << nwfs << endl;
  cout << "WFs selected " << mu_wfs.size() << "/" << nwfs << endl;

 
  //----- SELECTION LOOPS END ---------------------------------------
  //-----------------------------------------------------------------

  TH2D* h2_fprompt_charge_muwfs = BuildChargeFpromptHisto(mu_wfs, int_low, int_up, int_prompt);
  TH1D* h_fprompt_muwfs = h2_fprompt_charge_muwfs->ProjectionY("h_prompt_muwfs");



  cout << "\nMuon canidates: " << mu_wfs.size() << endl;
  avgWF(mu_wfs, avg_mu);

  if(print==true) {
    string outfile_name;
    cout << "Name of the muon file (without .dat)" << endl;
    cin >> outfile_name;
    outfile_name = muon_files_path+outfile_name;
    VecDouble_in_Binary(outfile_name, avg_mu);
    print = false;
  }

  double ymin, ymax;
  min_max_element(mu_wfs, ymin, ymax);

  if (plot == true){
    TH2D* h2_mu_pers = new TH2D("h2_mu_pers", Form("%s;%s;%s", "Persistence", "Ticks", "ADC Counts"),
        memorydepth/2, 0., memorydepth, 
        120, ymin, ymax);
    
    for (int i=0; i< mu_wfs.size(); i++) 
      for (int j=0; j<memorydepth; j=j+2) h2_mu_pers->Fill(j, mu_wfs[i][j]);

    
    TCanvas *c_mu_pers = new TCanvas("c_mu_pers","c_mu_pers",0,0,500,450);  
    c_mu_pers->cd();
    h2_mu_pers->Draw("COLZ");
    c_mu_pers->SetLogz();
    c_mu_pers->Modified(); c_mu_pers->Update();


    TCanvas *c_fprompt_charge_wfs = new TCanvas("c_fprompt_charge_wfs","c_fprompt_charge_wfs",500,0,500,450);
    c_fprompt_charge_wfs->cd();
    h2_fprompt_charge_wfs->Draw();
    c_fprompt_charge_wfs->Modified(); c_fprompt_charge_wfs->Update();

    TCanvas *c_fprompt_charge_muwfs = new TCanvas("c_fprompt_charge_muwfs","c_fprompt_charge_muwfs",500,0,500,450);
    c_fprompt_charge_muwfs->cd();
    h2_fprompt_charge_muwfs->Draw();
    c_fprompt_charge_muwfs->Modified(); c_fprompt_charge_muwfs->Update();

    TGraph *g_mu_avg = new TGraph(avg_mu.size(), &avg_mu[0]);
    //g_mu_avg->GetXaxis()->SetTitle("Time [#mus]");
    g_mu_avg->GetYaxis()->SetTitle("ADC counts");
    g_mu_avg->SetLineColor(2);
    TCanvas *c_mu_avg = new TCanvas("c_mu_avg","c_mu_avg",0,550,500,450);
    c_mu_avg->cd();
    g_mu_avg->Draw();
    c_mu_avg->Modified(); c_mu_avg->Update();

    TCanvas *c_fprompt = new TCanvas("c_fprompt","c_fprompt",500,550,500,450);
    c_fprompt->cd();
    h_fprompt_wfs->SetLineColor(kRed);
    h_fprompt_wfs->Draw();
    h_fprompt_muwfs->Draw("SAME");
    c_fprompt->Modified(); c_fprompt->Update();
  }
}

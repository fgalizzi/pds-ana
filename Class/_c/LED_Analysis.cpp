#include "../classe.hpp"

// ****************************************************************
// Description
// ****************************************************************

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void cla::LED_Analysis(){
  if (!keep_hcharge){
    vector<vector<double>> sel_wf; // Wavefroms I select for the analysis
    vector<double> int_wf;         // Integrals of the selected waveforms

    // Read and subtract the baseline
    read();

    // Select the with a baseline within the [-bsl; bsl] range and 
    // fully contained in the [sat_low; sat_up] range
    SelCalib_WF(wfs, sel_wf, prepulse_ticks, sat_low, sat_up, bsl);  
    if (sel_wf.size() == 0){
      class_skip = 1;
      return;
    }
    h_charge = new TH1D();

    if (manual == 0) {
      if(mov_win == true){
        MovingAverageWF(sel_wf, sel_wf, win);   
        h_charge = BuildRawChargeHisto(sel_wf, int_wf, int_low+win, int_up+win, nbins);
      } else h_charge = BuildRawChargeHisto(sel_wf, int_wf, int_low, int_up, nbins);
    }
    else {
      if(mov_win == true){
        MovingAverageWF(sel_wf, sel_wf, win);   
        h_charge = BuildRawChargeHisto(sel_wf, int_wf, int_low+win, int_up+win,
            hmin, hmax, nbins);
      } else h_charge = BuildRawChargeHisto(sel_wf, int_wf, int_low, int_up,
          hmin, hmax, nbins);
    }
  }
  
  TH1D* hFind = new TH1D(*h_charge); //Delcared to plot it later
  
  // See private_methods.hpp
  // "fgaus" is the function to fit the chargehistogram, it is initialized 
  // by looking for peaks in the histo or manually by setting manual=true
  fgaus = set_charge_fit_function(h_charge, hFind);
  cout << "\n\n----------- FIT CHARGE HISTO ----------------------\n" << endl;
  h_charge->Fit(fgaus, "QNR");
  // Re-fit in case a peak is fallen into a valley
  if (fgaus->GetParameter(5) < 0.5*fgaus->GetParameter(4) && fgaus->GetParameter(5) < 0.5*fgaus->GetParameter(6)){
    fgaus->SetParameter(5, fgaus->GetParameter(6));
    fgaus->SetParameter(1, 2*fgaus->GetParameter(1));
    h_charge->Fit(fgaus, "QNR");
  }
  // Compare the x axis upper limit of the histo h_charge and 6*fgaus->GetParameter(1)
  // set the upper limit of fgaus range to the minimum of the two
  spe_charge  = fgaus->GetParameter(1);
  sigma_zero  = fgaus->GetParameter(2);
  double last_fitted_peak_pe; // last fitted peak in photoelectrons
  if (h_charge->GetXaxis()->GetXmax() < (nmaxpeaks-1)*spe_charge+sigma_zero){
    std::cout << "\n\n\n HERE \n\n\n" << std::endl;
    fgaus->SetRange(-2.5*sigma_zero, h_charge->GetXaxis()->GetXmax());
    last_fitted_peak_pe = h_charge->GetXaxis()->GetXmax()/(spe_charge);
  }
  else{
    std::cout << "\n\n\n HpE \n\n\n" << std::endl;
    fgaus->SetRange(-2.5*sigma_zero, (nmaxpeaks-1)*spe_charge+sigma_zero);
    last_fitted_peak_pe = nmaxpeaks;
  }

  TFitResultPtr FitRes = h_charge->Fit(fgaus, "RMES");
  fgaus->SetRange(fgaus->GetXmin()*4, fgaus->GetXmax()); // Extend the range for the CX fit
  cout << "\n---------------------------------------------------  \n" << endl;

  if(display==true){
    h_charge->Draw();fgaus->Draw("same");gPad->Update();//gPad->WaitPrimitive();
    std::cout << "Want to repeat with manual mode?: 0=no - 1=yes" << endl;
    cin >> manual;
    if(manual==1){
      cout << "mu0_low" << endl;
      cin  >> mu0_low;
      cout << "mu0_up" << endl;
      cin  >> mu0_up;
      cout << "spe_low" << endl;
      cin  >> spe_low;
      cout << "spe_up" << endl;
      cin  >> spe_up;
      cout << "s0_low" << endl;
      cin  >> s0_low;
      cout << "s0_up" << endl;
      cin  >> s0_up;
      cout << "sc_low" << endl;
      cin  >> sc_low;
      cout << "sc_up" << endl;
      cin  >> sc_up;
      fgaus->SetParameter(0,(mu0_low+mu0_up)/2.);fgaus->SetParLimits(0,mu0_low,mu0_up);
      fgaus->SetParameter(1,(spe_low+spe_up/2.));fgaus->SetParLimits(1,spe_low,spe_up);
      fgaus->SetParameter(2,(s0_low+s0_up)/2.);fgaus->SetParLimits(2,s0_low,s0_up);
      fgaus->SetParameter(3,(sc_low+sc_up)/2.);fgaus->SetParLimits(3,sc_low,sc_up);

      FitRes = h_charge->Fit("fgaus", "RS");
    }
  }

  // --- CROSS TALK ----------------------------------------------
  auto g_CX = Build_CX_Graph(fgaus, h_charge, avg_n_photons, nmaxpeaks);
  // auto g_CX = Build_CX_Graph_Cov(fgaus, h_charge, FitRes, avg_n_photons); // With covariance matrix
                                                                             // seems to overestimate the errors
  TF1* f_CX = new TF1("f_CX", fCX, -0.5, last_fitted_peak_pe+0.5, 3);
  fCX_set(f_CX);
  cout << "\n\n----------- FIT CROSS-TALK ------------------------\n" << endl;
  g_CX->Fit("f_CX", "R");
  cout << "\n---------------------------------------------------  \n" << endl;


  // --- CLASS UPDATE --------------------------------------------
  pedestal    = fgaus->GetParameter(0);
  spe_charge  = fgaus->GetParameter(1); err_spe_charge = fgaus->GetParError(1);
  sigma_zero  = fgaus->GetParameter(2); err_sigma_zero = fgaus->GetParError(2);
  SNR         = spe_charge/sigma_zero;  err_SNR = error_propagation(FitRes, fgaus, 1, 2, "div");
  avg_n_ph_cx = f_CX->GetParameter(0);  err_avg_n_ph_cx = f_CX->GetParError(0);
  cx          = f_CX->GetParameter(1);  err_cx = f_CX->GetParError(1);
  avg_n_photoelectrons = (h_charge->GetMean()-fgaus->GetParameter(0)) / fgaus->GetParameter(1);
  
  

  // --- STD::COUT -----------------------------------------------
  cout << "\n---------------------------------------------------  \n" << endl;
  cout << "Results:"  << endl;
  cout << "SNR = " << SNR << endl;
  cout << "CX  = " << cx << "+/-"  << err_cx << endl;
  cout << "#Avg photons CX fit = " << avg_n_ph_cx << " +/- " << err_avg_n_ph_cx << endl;
  cout << "#Avg photons        = " << avg_n_photons << endl;
  cout << "#Avg photoelectrons = " << avg_n_photoelectrons << endl;
  cout << "\n---------------------------------------------------  \n" << endl;

  // --- PLOT ----------------------------------------------------
  if (plot == true) {
    TCanvas *c1 = new TCanvas("c1","c1",0,0,1000,900);
    c1->cd();
    g_CX->GetXaxis()->SetTitle("Peaks");
    g_CX->GetYaxis()->SetTitle("Probability");
    g_CX->Draw("ap");
    c1->Modified();c1->Update();
    
    TCanvas *c2 = new TCanvas("c2","c2",20,20,1000,900);
    c2->cd();
    h_charge->Draw();
    c2->Modified();c2->Update();
  }
  
}

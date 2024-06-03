
void cla::update(){
  std::ofstream outputFile("const.hpp");

  if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return;
    }
	outputFile << "std::string WF_FILE = \"" << wf_file  << "\";" << std::endl;
	outputFile << "std::string TEMPL_F = \"" << templ_f << "\";" << std::endl;
	outputFile << "std::string NOISE_F = \"" << noise_f << "\";" << std::endl;
	outputFile << "std::string MUON_F  = \"" << muon_f  << "\";" << std::endl;
	outputFile << "std::string TRG_F  = \"" << trg_f  << "\";" << std::endl;
	outputFile << "// daphne caen csv " << std::endl;
	outputFile << "std::string DATA_FORMAT  = \"" << data_format << "\";" << std::endl;
	outputFile << "bool INVERT              = " << invert << ";" << std::endl;
  outputFile << "const int MEMORYDEPTH    = " << memorydepth << ";" << std::endl;
	outputFile << "const int N_WF           = " << n_wf << ";" << std::endl;
	outputFile << "const int RES            = " << res << ";" << std::endl;
	outputFile << "const int PREPULSE_TICKS = " << prepulse_ticks << ";" << std::endl;
	outputFile << "const int PRETRG         = " << pretrg << ";" << std::endl;
	outputFile << "const int AFTTRG         = " << afttrg << ";" << std::endl;
	outputFile << "const double TICK_LEN    = " << tick_len << ";" << std::endl;
	outputFile << " " << std::endl;
	outputFile << "const double SAT_UP      = " << sat_up << ";" << std::endl;
	outputFile << "const double SAT_LOW     = " << sat_low << ";" << std::endl;
	outputFile << "const double BSL         = " << bsl << ";" <<  std::endl;
	outputFile << " " << std::endl;
	outputFile << " " << std::endl;
	outputFile << "//calibration"  << std::endl;
	outputFile << "const int INT_LOW        = " << int_low << ";" << std::endl;
	outputFile << "const int INT_UP         = " << int_up << ";" << std::endl;
	outputFile << "const int NBINS          = " << nbins << ";" << std::endl;
	outputFile << "const int NMAXPEAKS      = " << nmaxpeaks << ";" << std::endl;
	outputFile << "const double MU0_LOW     = " << mu0_low << ";" << std::endl;
	outputFile << "const double MU0_UP      = " << mu0_up << ";" << std::endl;
	outputFile << "const double SPE_LOW     = " << spe_low << ";" << std::endl;
	outputFile << "const double SPE_UP      = " << spe_up << ";" << std::endl;
	outputFile << "const double S0_LOW      = " << s0_low << ";" << std::endl;
	outputFile << "const double S0_UP       = " << s0_up << ";" << std::endl;
	outputFile << "const double SC_LOW      = " << sc_low << ";" << std::endl;
	outputFile << "const double SC_UP       = " << sc_up << ";" << std::endl;
	outputFile << "const double FIT_LOW     = " << fit_low << ";" << std::endl;
	outputFile << "const double FIT_UP      = " << fit_up << ";" << std::endl;
	outputFile << "const double HMIN        = " << hmin << ";" << std::endl;
	outputFile << "const double HMAX        = " << hmax << ";" << std::endl;
	outputFile << "const double SPE_AMPL    = " << spe_ampl << ";" << std::endl;
	outputFile << "const double SPE_CHARGE  = " << spe_charge << ";" << std::endl;
	outputFile << "const double PEDESTAL    = " << pedestal << ";" << std::endl;
	outputFile << "//deconvolution"  << std::endl;
	outputFile << "const int INT_PROMPT     = " << int_prompt << ";"  << std::endl;
	outputFile << "const double F_PROMPT    = " << f_prompt << ";"  << std::endl;
	outputFile << "const double AMP_LOW     = " << amp_low<< ";"  << std::endl;
	outputFile << "const double AMP_UP      = " << amp_up << ";"  << std::endl;
	outputFile << "const double DECO_SM     = " << deco_sm << ";"  << std::endl;
	outputFile << "const double N2_         = " << n2_ << ";"  << std::endl;
	outputFile << "const double FIT_L       = " << fit_l << ";"  << std::endl;
	outputFile << "const double FIT_U       = " << fit_u << ";"  << std::endl;
	outputFile << "const double A_FAST      = " << a_fast << ";"  << std::endl;
	outputFile << "const double TAU_FAST    = " << tau_fast << ";"  << std::endl;
	outputFile << "const double A_SLOW      = " << a_slow << ";"  << std::endl;
	outputFile << "const double TAU_SLOW    = " << tau_slow << ";"  << std::endl;
	outputFile << "const double SIGMA       = " << sigma << ";"  << std::endl;
	outputFile << "const double T_0         = " << t_0 << ";"  << std::endl;

  
  outputFile.close();
}

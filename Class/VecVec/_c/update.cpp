
void cla::update(){
  std::ofstream outputFile("_c/const.hpp");

  if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return;
    }
	outputFile << "std::string wf_file = \"" << WF_FILE  << "\";" << std::endl;
	outputFile << "// daphne caen csv " << std::endl;
	outputFile << "std::string data_format  = \"" << DATA_FORMAT << "\";" << std::endl;
	outputFile << "bool invert              = " << INVERT << ";" << std::endl;
  outputFile << "const int memorydepth    = " << MEMORYDEPTH << ";" << std::endl;
	outputFile << "const int n_wf           = " << N_WF << ";" << std::endl;
	outputFile << "const int res            = " << RES << ";" << std::endl;
	outputFile << "const int prepulse_ticks = " << PREPULSE_TICKS << ";" << std::endl;
	outputFile << "const double tick_len    = " << TICK_LEN << ";" << std::endl;
	outputFile << " " << std::endl;
	outputFile << "const double sat_up      = " << SAT_UP << ";" << std::endl;
	outputFile << "const double sat_low     = " << SAT_LOW << ";" << std::endl;
	outputFile << "const double bsl         = " << BSL << ";" <<  std::endl;
	outputFile << " " << std::endl;
	outputFile << " " << std::endl;
	outputFile << "//calibration"  << std::endl;
	outputFile << "const int int_low        = " << INT_LOW << ";" << std::endl;
	outputFile << "const int int_up         = " << INT_UP << ";" << std::endl;
	outputFile << "const int nbins          = " << NBINS << ";" << std::endl;
	outputFile << "const int nmaxpeaks      = " << NMAXPEAKS << ";" << std::endl;
	outputFile << "const double mu0_low     = " << MU0_LOW << ";" << std::endl;
	outputFile << "const double mu0_up      = " << MU0_UP << ";" << std::endl;
	outputFile << "const double spe_low     = " << SPE_LOW << ";" << std::endl;
	outputFile << "const double spe_up      = " << SPE_UP << ";" << std::endl;
	outputFile << "const double s0_low      = " << S0_LOW << ";" << std::endl;
	outputFile << "const double s0_up       = " << S0_UP << ";" << std::endl;
	outputFile << "const double sc_low      = " << SC_LOW << ";" << std::endl;
	outputFile << "const double sc_up       = " << SC_UP << ";" << std::endl;
	outputFile << "const double fit_low     = " << FIT_LOW << ";" << std::endl;
	outputFile << "const double fit_up      = " << FIT_UP << ";" << std::endl;
	outputFile << "const double hmin        = " << HMIN << ";" << std::endl;
	outputFile << "const double hmax        = " << HMAX << ";" << std::endl;
	outputFile << "//deconvolution"  << std::endl;
	outputFile << "const int int_prompt     = 0;"  << std::endl;
	outputFile << "const double f_prompt    = 0;"  << std::endl;
  
  outputFile.close();
}

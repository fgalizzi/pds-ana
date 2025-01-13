#include "../classe.hpp"

// ****************************************************************
// Update the const.hpp file containing the class member initial
// values with the current ones
// ****************************************************************

//-----------------------------------------------------------------
//------- Macro ---------------------------------------------------
void cla::update(string out_file_name){
  out_file_name = out_file_name + ".hpp";
  std::ofstream outputFile(out_file_name);

  if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return;
    }
	outputFile << "WF_FILE = \"" << wf_file  << "\";" << std::endl;
	outputFile << "TEMPL_F = \"" << templ_f << "\";" << std::endl;
	outputFile << "NOISE_F = \"" << noise_f << "\";" << std::endl;
	outputFile << "MUON_F  = \"" << muon_f  << "\";" << std::endl;
	outputFile << "TRG_F   = \"" << trg_f  << "\";" << std::endl;
	outputFile << "// daphne caen csv " << std::endl;
	outputFile << "DATA_FORMAT    = \"" << data_format << "\";" << std::endl;
	outputFile << "INVERT         = " << invert << ";" << std::endl;
  outputFile << "MEMORYDEPTH    = " << memorydepth << ";" << std::endl;
	outputFile << "N_WF           = " << n_wf << ";" << std::endl;
	outputFile << "RES            = " << res << ";" << std::endl;
	outputFile << "PREPULSE_TICKS = " << prepulse_ticks << ";" << std::endl;
	outputFile << "PRETRG         = " << pretrg << ";" << std::endl;
	outputFile << "AFTTRG         = " << afttrg << ";" << std::endl;
	outputFile << "TICK_LEN       = " << tick_len << ";" << std::endl;
	outputFile << " " << std::endl;
	outputFile << "SAT_UP      = " << sat_up << ";" << std::endl;
	outputFile << "SAT_LOW     = " << sat_low << ";" << std::endl;
	outputFile << "BSL         = " << bsl << ";" <<  std::endl;
	outputFile << "RMS         = " << rms << ";" <<  std::endl;
	outputFile << " " << std::endl;
	outputFile << " " << std::endl;
	outputFile << "//calibration"  << std::endl;
	outputFile << "INT_LOW        = " << int_low << ";" << std::endl;
	outputFile << "INT_UP         = " << int_up << ";" << std::endl;
	outputFile << "NBINS          = " << nbins << ";" << std::endl;
	outputFile << "NMAXPEAKS      = " << nmaxpeaks << ";" << std::endl;
	outputFile << "MU0_LOW        = " << mu0_low << ";" << std::endl;
	outputFile << "MU0_UP         = " << mu0_up << ";" << std::endl;
	outputFile << "SPE_LOW        = " << spe_low << ";" << std::endl;
	outputFile << "SPE_UP         = " << spe_up << ";" << std::endl;
	outputFile << "S0_LOW         = " << s0_low << ";" << std::endl;
	outputFile << "S0_UP          = " << s0_up << ";" << std::endl;
	outputFile << "SC_LOW         = " << sc_low << ";" << std::endl;
	outputFile << "SC_UP          = " << sc_up << ";" << std::endl;
	outputFile << "FIT_LOW        = " << fit_low << ";" << std::endl;
	outputFile << "FIT_UP         = " << fit_up << ";" << std::endl;
	outputFile << "HMIN           = " << hmin << ";" << std::endl;
	outputFile << "HMAX           = " << hmax << ";" << std::endl;
	outputFile << "SPE_AMPL       = " << spe_ampl << ";" << std::endl;
	outputFile << "SPE_CHARGE     = " << spe_charge << ";" << std::endl;
	outputFile << "PEDESTAL       = " << pedestal << ";" << std::endl;
	outputFile << "//deconvolution"  << std::endl;
	outputFile << "INT_PROMPT     = " << int_prompt << ";"  << std::endl;
	outputFile << "ROLL           = " << roll << ";"  << std::endl;
	outputFile << "AMP            = " << amp << ";"  << std::endl;
	outputFile << "F_PROMPT       = " << f_prompt << ";"  << std::endl;
	outputFile << "AMP_LOW        = " << amp_low<< ";"  << std::endl;
	outputFile << "AMP_UP         = " << amp_up << ";"  << std::endl;
	outputFile << "DECO_SM        = " << deco_sm << ";"  << std::endl;
	outputFile << "N2_            = " << n2_ << ";"  << std::endl;
	outputFile << "FIT_L          = " << fit_l << ";"  << std::endl;
	outputFile << "FIT_U          = " << fit_u << ";"  << std::endl;
	outputFile << "A_FAST         = " << a_fast << ";"  << std::endl;
	outputFile << "TAU_FAST       = " << tau_fast << ";"  << std::endl;
	outputFile << "A_SLOW         = " << a_slow << ";"  << std::endl;
	outputFile << "TAU_SLOW       = " << tau_slow << ";"  << std::endl;
	outputFile << "SIGMA          = " << sigma << ";"  << std::endl;
	outputFile << "T_0            = " << t_0 << ";"  << std::endl;
  outputFile << "YERR           = " << yerr << ";"  << std::endl;
	outputFile << "//ProtoDUNE"  << std::endl;
  outputFile << "CHANNEL        = " << channel << ";" << std::endl; 

  
  outputFile.close();
}

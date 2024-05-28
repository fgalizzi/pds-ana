// Here
//
//*** MAIN ************************************
void cla::configDCR(){
//*********************************************
  gStyle->SetOptFit(1111); gStyle->SetOptTitle(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit"); 

  std::vector<double> x, avg, y2, y3, avg_wf, int_wf;
  avg.resize(MEMORYDEPTH);
  for (size_t i = 0; i < MEMORYDEPTH; i++) x.push_back( (double) i);

  std::vector< std::vector<double> > wf2, grad, wf;
  wf2.resize(N_WF);
  wf.resize(N_WF);
  grad.resize(N_WF);

  // Read and subtract the baseline
  read();

  SelCalib_WF(WFS, y2, MEMORYDEPTH, PREPULSE_TICKS, SAT_LOW, SAT_UP, BSL);
  TH1D* hI = BuildRawChargeHisto(y2, int_wf, MEMORYDEPTH, INT_LOW, INT_UP);
  Avg_Sel_WF (y2, y3, avg_wf, int_wf, SPE_LOW, SPE_UP);
 

  for(int i=0; i<y3.size()/MEMORYDEPTH; i++){
    for(int j=0; j<MEMORYDEPTH; j++) wf2[i].push_back(y3[i*MEMORYDEPTH+j]);
    //for(int j=0; j<MEMORYDEPTH; j++) wf[i].push_back(WFS[i*MEMORYDEPTH+j]);
   
//    for(int j=0; j<ITE; j++) MovingAverageWF(wf2[i], wf2[i], WIN);
    wf[i].resize(MEMORYDEPTH);
    TV1D_denoise(wf2[i], wf[i], wf2[i].size(), DEN);

    grad[i].resize(MEMORYDEPTH);
    vector_gradient(wf[i], grad[i]);
    for(int j=0; j<ITE; j++) MovingAverageWF(grad[i], grad[i], WIN);
  
  }

  for(int i=0; i<10; i++) DisplayWFs (wf[i], grad[i], MEMORYDEPTH, TICK_LEN, 1); 
  //for(int i=0; i<10; i++) DisplayWFs (grad[i], MEMORYDEPTH, TICK_LEN, 1); 
  //for(int i=0; i<N_WF; i++) for (int j=0; j<MEMORYDEPTH; j++ ) avg[j] +=
  //   grad[i][j];

//  DisplayWFs(avg, MEMORYDEPTH, TICK_LEN, 1);
}

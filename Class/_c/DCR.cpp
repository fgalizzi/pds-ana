// Here
//

void TV1D_denoise(vector<Double_t>& input, vector<Double_t>& output, const int width, const double lambda) {
  if (width>0) {                /*to avoid invalid memory access to input[0]*/
    int k=0, k0=0;            /*k: current sample location, k0: beginning of current segment*/
    double umin=lambda, umax=-lambda;    /*u is the dual variable*/
    double vmin=input[0]-lambda, vmax=input[0]+lambda;    /*bounds for the segment's value*/
    int kplus=0, kminus=0;  /*last positions where umax=-lambda, umin=lambda, respectively*/
    const double twolambda=2.0*lambda;    /*auxiliary variable*/
    const double minlambda=-lambda;        /*auxiliary variable*/
    for (;;) {                /*simple loop, the exit test is inside*/
      while (k==width-1) {    /*we use the right boundary condition*/
        if (umin<0.0) {            /*vmin is too high -> negative jump necessary*/
          do output[k0++]=vmin; while (k0<=kminus);
          umax=(vmin=input[kminus=k=k0])+(umin=lambda)-vmax;
        } else if (umax>0.0) {    /*vmax is too low -> positive jump necessary*/
          do output[k0++]=vmax; while (k0<=kplus);
          umin=(vmax=input[kplus=k=k0])+(umax=minlambda)-vmin;
        } else {
          vmin+=umin/(k-k0+1);
          do output[k0++]=vmin; while(k0<=k);
          return;
        }
      }
      if ((umin+=input[k+1]-vmin)<minlambda) {        /*negative jump necessary*/
        do output[k0++]=vmin; while (k0<=kminus);
        vmax=(vmin=input[kplus=kminus=k=k0])+twolambda;
        umin=lambda; umax=minlambda;
      } else if ((umax+=input[k+1]-vmax)>lambda) {    /*positive jump necessary*/
        do output[k0++]=vmax; while (k0<=kplus);
        vmin=(vmax=input[kplus=kminus=k=k0])-twolambda;
        umin=lambda; umax=minlambda;
      } else {        /*no jump necessary, we continue*/
        k++;
        if (umin>=lambda) {        /*update of vmin*/
          vmin+=(umin-lambda)/((kminus=k)-k0+1);
          umin=lambda;
        }
        if (umax<=minlambda) {    /*update of vmax*/
          vmax+=(umax+lambda)/((kplus=k)-k0+1);
          umax=minlambda;
        }
      }
    }
  }
}

void vector_gradient(std::vector<double>& wf, std::vector<double>& grad){
  double x1 = 0;
  double x2 = 0;
  int N = wf.size();

  for(int i=0; i<N-1; i++) grad[i] = (wf[i+1]-wf[i])*10.;
  grad[N] = 0;

}

//*** MAIN ************************************
void cla::DCR(){
//*********************************************
  /*gStyle->SetOptFit(1111); gStyle->SetOptTitle(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit"); 

  vector<double> x, avg_wf, int_wf;
  vector<vector<double>> y2, y3;
  avg_wf.resize(MEMORYDEPTH);
  x.resize(MEMORYDEPTH);

  for (size_t i = 0; i < MEMORYDEPTH; i++) x[i] = double(i);

  vector<vector<double>> wf2, grad, wf;
  wf2.resize(WFS.size(), vector<double>(MEMORYDEPTH));
  wf.resize(WFS.size(), vector<double>(MEMORYDEPTH));
  grad.resize(WFS.size(), vector<double>(MEMORYDEPTH)));

  // Read and subtract the baseline
  read();

  SelCalib_WF(WFS, y2, PREPULSE_TICKS, SAT_LOW, SAT_UP, BSL);
  TH1D* hI = BuildRawChargeHisto(y2, int_wf, INT_LOW, INT_UP, NBINS);
  Avg_Sel_WF (y2, y3, avg_wf, int_wf, SPE_LOW, SPE_UP);
 

  for(auto wf : WFS){
    for(int j=0; j<MEMORYDEPTH; j++) wf2[i].push_back(wf[j]);
   
//    for(int j=0; j<ITE; j++) MovingAverageWF(wf2[i], wf2[i], WIN);
    TV1D_denoise(wf2[i], wf[i], wf2[i].size(), DEN);
    vector_gradient(wf[i], grad[i]);
    for(int j=0; j<ITE; j++) MovingAverageWF(grad[i], grad[i], WIN);
  
  }

  y2.erase(y2.begin(), y2.end());

  // y2 contiene il gradiente di WFS
  for(auto v : grad) y2.insert(y2.end(), v.begin(), v.end());
*/
  printf("\n\n\nAL MOMENTO E' SOLO UN COMMENTO");
}

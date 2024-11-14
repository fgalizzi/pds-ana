
void double_expo(double* x, const double* p, size_t len){
  double a_f = p[0];
  double t_f = p[1];
  double a_s = p[2];
  double t_s = p[3];
  double sig = p[4];
  double t_0 = p[5];
  double con = p[6];

  for(size_t i=0; i<len; i++){
    double tt = double(i);
    double res = (a_f*exp(-(tt-t_0)/t_f)*exp(sig*sig/(2*t_f*t_f)))*TMath::Erfc(((t_0-tt)/sig+sig/t_f)/TMath::Power(2,0.5))/2. +
               (a_s*exp(-(tt-t_0)/t_s)*exp(sig*sig/(2*t_s*t_s)))*TMath::Erfc(((t_0-tt)/sig+sig/t_s)/TMath::Power(2,0.5))/2. +
                con;
    // if (isnan(res)==true || res<1.e-40) res = 0.;

    x[i] = res;
  }

}


void double_expo(double* x, const double* p, size_t len){
  double a_f = p[0];
  double t_f = p[1];
  double a_s = p[2];
  double t_s = p[3];
  double sig = p[4];
  double t_0 = p[5];
  double con = p[6];

  double tf_inv = 1./t_f;
  double ts_inv = 1./t_s;
  double sig_inv = 1./sig;
  double sig_tf = sig*tf_inv;
  double sig_ts = sig*ts_inv;
  double sig2_tf = sig*sig/(2*t_f*t_f);
  double sig2_ts = sig*sig/(2*t_s*t_s);
  double sqrt2_inv = 1.0 / TMath::Sqrt2();

  for(size_t i=0; i<len; i++){
    double td = t_0-double(i);
    double res = a_f*exp(td*tf_inv + sig2_tf)*TMath::Erfc((td*sig_inv+sig_tf)*sqrt2_inv)*0.5 +
                 a_s*exp(td*ts_inv + sig2_ts)*TMath::Erfc((td*sig_inv+sig_ts)*sqrt2_inv)*0.5 +
                 con;
    // if (isnan(res)==true || res<1.e-40) res = 0.;
    
    x[i] = res;
  }

}

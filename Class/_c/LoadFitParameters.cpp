
void cla::LoadFitParameters(TF1* f) {
  double mu0  = f->GetParameter(0);
  double gain = f->GetParameter(1);
  double s0   = f->GetParameter(2);
  mu0_low     = mu0-s0*2;
  mu0_up      = mu0+s0;
  spe_low     = mu0+gain-s0;
  spe_up      = mu0+gain+s0;

  return;
};


void cla::LoadFitParameters(TF1* fgaus) {
  double mu0 = fgaus->GetParameter(0);
  double gain= fgaus->GetParameter(1);
  double s0  = fgaus->GetParameter(2);
  MU0_LOW     = mu0-s0*2;
  MU0_UP      = mu0+s0;
  SPE_LOW     = mu0+gain-s0;
  SPE_UP      = mu0+gain+s0;

  return;
};

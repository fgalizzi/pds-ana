
void cla::TooAnnoying(){
  manual = 1;
  fit_low = -1.8;
  fit_up = 6.8;
  SelfTrigger();
  fit_low = 1444;
  fit_up = 1453;
  Jitter();
}

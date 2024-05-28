
void cla::TooAnnoying(){
  LoadFitParameters(fgaus);  
  AverageWF();
  SPE();
  mov_win = 1;
  win = 20;
  AverageWF();
  SPE();
}

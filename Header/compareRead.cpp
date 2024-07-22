// compare speed 

void compareRead(){
  std::string fileName = "Ch1.dat";
  std::vector<double> y;

  y.resize(1e6, 0.);
  size_t len = y.size();
  double t=0;

  auto start = std::chrono::high_resolution_clock::now();
  for(int i=0; i<10; i++) for(size_t j=0; j<len; j++) t += y[j]; 
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end - start;
  std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;

  start = std::chrono::high_resolution_clock::now();
  for(int i=0; i<10; i++) t = std::reduce(y.begin(), y.end());
  end = std::chrono::high_resolution_clock::now();
  duration = end - start;
  std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;



}

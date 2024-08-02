// compare speed
//

#include "G_Utility.hpp"

void compareRead(){

  std::vector<double> vec;
  for(size_t i=0; i<50; i++)
    for(size_t j=0; j<i; j++) vec.push_back(double(i));

  double t;

  auto start = std::chrono::high_resolution_clock::now();
  for(int i=0; i<100; i++) t = Vector_MPV(vec); 
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end - start;
  std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;
  std::cout << t << std::endl;
  std::cout << *std::max_element(vec.begin(), vec.end()) << std::endl;
  // start = std::chrono::high_resolution_clock::now();
  // for(int i=0; i<10000; i++) t = Vector_MPV2(vec); 
  // end = std::chrono::high_resolution_clock::now();
  // duration = end - start;
  // std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;
  //
}

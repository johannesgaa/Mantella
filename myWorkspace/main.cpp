#include "orbitalMechanics.hpp"

int main(){
  
  std::cout << " --- stumpffC --- " << std::endl;
  std::cout << "stumpffC(100) " << mant::itd::stumpffC(100) << std::endl;
  std::cout << "stumpffC(-100) " << mant::itd::stumpffC(-100) << std::endl;
  std::cout << "stumpffC(0) " << mant::itd::stumpffC(0) << std::endl << std::endl;

  std::cout << " --- stumpffS --- " << std::endl;  
  std::cout << "stumpffS(100) " << mant::itd::stumpffS(100) << std::endl;
  std::cout << "stumpffS(-100) " << mant::itd::stumpffS(-100) << std::endl;
  std::cout << "stumpffS(0) " << mant::itd::stumpffS(0) << std::endl << std::endl;
  
  std::cout << " --- Lambert --- " << std::endl;
  std::pair<arma::Col<double>, arma::Col<double>> velocities = mant::itd::lambert(
    {5000.0, 10000.0, 2100.0}, 
    {-14600.0, 2500.0, 7000.0}, 
    3600.0, 
    true); 
  std::cout << "v1 " << std::endl << velocities.first << std::endl;
  std::cout << "v2 " << std::endl << velocities.second << std::endl << std::endl;
  
  std::cout << " --- Planet --- " << std::endl;
  
  return 0;
}

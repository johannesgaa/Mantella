#include "orbitalMechanics.hpp"

int main(){
  
  std::cout << "stumpffC(100) " << mant::itd::stumpffC(100) << std::endl;
  std::cout << "stumpffC(-100) " << mant::itd::stumpffC(-100) << std::endl;
  std::cout << "stumpffC(0) " << mant::itd::stumpffC(0) << std::endl << std::endl;
  
  std::cout << "stumpffS(100) " << mant::itd::stumpffS(100) << std::endl;
  std::cout << "stumpffS(-100) " << mant::itd::stumpffS(-100) << std::endl;
  std::cout << "stumpffS(0) " << mant::itd::stumpffS(0) << std::endl;
  
  return 0;
}

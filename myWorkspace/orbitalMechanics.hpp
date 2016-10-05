#pragma once

// C++ standard library
#include <cmath>
#include <iostream>

// Armadillo
#include <armadillo>

// Mantella
#include "../include/mantella" //apt-get install liblapacke-dev

namespace mant{
  namespace itd{
    
    const double heliocentricGravitationalConstant = 1.32712440018e11; //km^3 / s^2
    
    double stumpffC(
        const double parameter){
      
      if(parameter > 0){
        
        return (1.0 - std::cos(std::sqrt(parameter))) / parameter;
        
      }else if (parameter < 0){
        
        return (std::cosh(std::sqrt(-parameter)) - 1.0) / -parameter;
        
      }else{
        
        return 0.5;
        
      }
    }
    
    double stumpffS(
        const double parameter){
      
      if(parameter > 0){
        
        return (std::sqrt(parameter) - std::sin(std::sqrt(parameter))) / std::pow(std::sqrt(parameter), 3.0);
        
      }else if (parameter < 0){
        
        return (std::sinh(std::sqrt(-parameter)) - std::sqrt(-parameter)) / std::pow(std::sqrt(-parameter), 3.0);
        
      }else{
        
        return 1.0 / 6.0;
        
      }
    }
    
    
  }
}







#pragma once

// C++ standard library
#include <cmath>
#include <iostream>
#include <functional>

// Armadillo
#include <armadillo>

// Mantella
#include <mantella>

namespace mant{
  namespace itd{
    
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
    
    std::pair<arma::Col<double>::fixed<3>, arma::Col<double>::fixed<3>> lambert(
        const arma::Col<double>::fixed<3> &startPosition,
        const arma::Col<double>::fixed<3> &endPosition,
        const double flightTime,
        const bool isPrograde){
      
      double startPositionNorm = arma::norm(startPosition);
      double endPositionNorm = arma::norm(endPosition);
      
      arma::Col<double> positionsCrossProduct = arma::cross(startPosition, endPosition);
      double theta = std::acos(arma::norm_dot(startPosition, endPosition));
      
      if(isPrograde){
        if(positionsCrossProduct(2) < 0){
          theta = 2.0 * arma::datum::pi - theta;
        }
      }else{
        if(positionsCrossProduct(2) >= 0){
          theta = 2.0 * arma::datum::pi - theta;
        }
      }

      double a = std::sin(theta) * std::sqrt((startPositionNorm * endPositionNorm) / (1.0 - std::cos(theta)));
      
      std::function<double(const double)> yFunction = [&](
          const double parameter){ 
        return startPositionNorm + endPositionNorm + a * ((parameter * stumpffS(parameter) - 1.0) / std::sqrt(stumpffC(parameter)));
      };
      
      double z = mant::brent(
        [&](
            const double parameter) { 

          return stumpffS(parameter) * std::pow(yFunction(parameter) / stumpffC(parameter), 3.0/2.0) + a * std::sqrt(yFunction(parameter)) - std::sqrt(3.986e+5) * flightTime; //heliocentric should be 1.32712440018e11
          
        }, -2.0, 2.0, 1000); //TODO: remove magic numbers
      
      double y = yFunction(z);

      double f = 1.0 - y / startPositionNorm;
      double g = a * std::sqrt(y / 3.986e+5); //heliocentric should be 1.32712440018e11
      double gDot = 1.0 - y / endPositionNorm;
      
      return {
        {(1.0 / g) * (endPosition - f * startPosition)}, 
        {(1.0 / g) * (gDot * endPosition - startPosition)}
      };
      
    }
    
  }
}







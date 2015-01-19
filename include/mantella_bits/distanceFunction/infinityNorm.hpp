#pragma once

// Mantella
#include <mantella_bits/distanceFunction.hpp>

namespace mant {
  class InfinityNorm : public DistanceFunction<double> {
    protected:
      double getDistanceImplementation(
          const arma::Col<double>& parameter) const noexcept override;

      arma::Col<double> getNeighbourImplementation(
          const arma::Col<double>& parameter,
          const arma::Col<double>& minimalDistance,
          const arma::Col<double>& maximalDistance) const noexcept override;
  };
}
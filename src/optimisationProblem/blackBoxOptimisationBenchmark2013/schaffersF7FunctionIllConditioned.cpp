#include <hop_bits/optimisationProblem/blackBoxOptimisationBenchmark2013/schaffersF7FunctionIllConditioned.hpp>

// Cereal
#include <cereal/archives/json.hpp>
#include <cereal/types/polymorphic.hpp>

namespace hop {
  namespace bbob2013 {
    double SchaffersF7FunctionIllConditioned::getObjectiveValueImplementation(
        const arma::Col<double>& parameter) const noexcept {
      const arma::Col<double>& z = arma::square(delta_ % (rotationQ_ * getAsymmetricTransformation(0.5, rotationR_ * (parameter - translation_))));
      const arma::Col<double>& s = arma::pow(z.subvec(0, z.n_elem - 2) + z.subvec(1, z.n_elem - 1), 0.25);

      return std::pow(arma::mean(s % (1.0 + arma::square(arma::sin(50.0 * arma::pow(s, 0.4))))), 2.0) + 10.0 * getPenality(parameter);
    }

    std::string SchaffersF7FunctionIllConditioned::to_string() const noexcept {
      return "SchaffersF7FunctionIllConditioned";
    }
  }
}

CEREAL_REGISTER_TYPE(hop::bbob2013::SchaffersF7FunctionIllConditioned)
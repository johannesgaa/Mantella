#include <hop_bits/optimisationProblem/blackBoxOptimisationBenchmark2015/schwefelFunction.hpp>

// Cereal
#include <cereal/archives/json.hpp>
#include <cereal/types/polymorphic.hpp>

namespace hop {
  namespace bbob2015 {
    double SchwefelFunction::getObjectiveValueImplementation(
        const arma::Col<double>& parameter) const noexcept {
      const arma::Col<double>& xOpt = arma::abs(4.2096874633 * one_);
      const arma::Col<double>& xHat = 2.0 * one_ % parameter;

      arma::Col<double> zHat(xHat.n_elem);
      zHat.at(0) = xHat.at(0);
      zHat.tail(zHat.n_elem - 1) = xHat.tail(zHat.n_elem - 1) + 0.25 * (xHat.head(xHat.n_elem - 1) - xOpt.head(xOpt.n_elem - 1));

      const arma::Col<double>& z = 100.0 * (delta_ % (zHat - xOpt) + xOpt);

      return 0.01 * (418.9828872724339 - arma::mean(z % arma::sin(arma::sqrt(arma::abs(z))))) + 100.0 * getPenality(z / 100.0);
    }

    std::string SchwefelFunction::to_string() const noexcept {
      return "SchwefelFunction";
    }
  }
}

CEREAL_REGISTER_TYPE(hop::bbob2015::SchwefelFunction)
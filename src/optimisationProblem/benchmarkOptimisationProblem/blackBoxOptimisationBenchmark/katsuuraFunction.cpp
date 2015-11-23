#include "mantella_bits/optimisationProblem/benchmarkOptimisationProblem/blackBoxOptimisationBenchmark/katsuuraFunction.hpp"

// C++ standard library
#include <cassert>
#include <cmath>

// Mantella
#include "mantella_bits/assert.hpp"
#include "mantella_bits/probability.hpp"

// This implementation contains a lot of *magic numbers* and behaviour, introduced by the black-box optimisation benchmark, but only partially explained in the paper.
// So don't expect use to explain the unexplained.
// @see N. Hansen et al., Real-Parameter Black-Box Optimization Benchmarking 2010: Experimental Setup. Research Report RR-7215, INRIA, 2010.
namespace mant {
  namespace bbob {
    KatsuuraFunction::KatsuuraFunction(
        const arma::uword numberOfDimensions)
        : BlackBoxOptimisationBenchmark(numberOfDimensions),
          parameterConditioning_(getParameterConditioning(10.0)),
          rotationQ_(randomRotationMatrix(numberOfDimensions_)) {
      setParameterTranslation(getRandomParameterTranslation());
      setParameterRotation(randomRotationMatrix(numberOfDimensions_));

#if defined(SUPPORT_MPI)
      MPI_Bcast(rotationQ_.memptr(), static_cast<int>(rotationQ_.n_elem), MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

      setObjectiveFunction([this](
                               const arma::Col<double>& parameter) {
          assert(parameter.n_elem == numberOfDimensions_);
            
          arma::Col<double> z = rotationQ_ * (parameterConditioning_ % parameter);

          double product = 1.0;
          for (arma::uword n = 0; n < z.n_elem; ++n) {
            const double value = z(n);

            double sum = 0.0;
            for (arma::uword k = 1; k < 33; ++k) {
              const double power = std::pow(2.0, k);
              sum += std::abs(power * value - std::round(power * value)) / power;
            }

            product *= std::pow(1.0 + (static_cast<double>(n) + 1.0) * sum, 10.0 / std::pow(numberOfDimensions_, 1.2));
          }

          return 10.0 / std::pow(numberOfDimensions_, 2.0) * (product - 1.0);
      },
          "BBOB Katsuura Function");
    }
  }
}

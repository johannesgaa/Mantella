// Catch
#include <catch.hpp>

// C++ Standard Library
#include <cstdlib>
#include <string>

// Armadillo
#include <armadillo>

// Boost
#include <boost/filesystem.hpp>

// HOP
#include <hop>

extern boost::filesystem::path testDirectory;

TEST_CASE("CompositeGriewankRosenbrockFunctionF8F2", "") {
  for (const auto& numberOfDimensions : {2, 40}) {
    hop::bbob2013::CompositeGriewankRosenbrockFunctionF8F2 compositeGriewankRosenbrockFunctionF8F2(numberOfDimensions);

    arma::Mat<double> parameters;
    parameters.load(testDirectory.string() + "/data/optimisationProblem/benchmark/blackBoxOptimisationBenchmark2013/parameters,dim" + std::to_string(numberOfDimensions) +".mat");

    arma::Mat<double> rotationR;
    rotationR.load(testDirectory.string() + "/data/optimisationProblem/benchmark/blackBoxOptimisationBenchmark2013/rotationR,dim" + std::to_string(numberOfDimensions) +".mat");

    arma::Col<double> expected;
    expected.load(testDirectory.string() + "/data/optimisationProblem/benchmark/blackBoxOptimisationBenchmark2013/expectedCompositeGriewankRosenbrockFunctionF8F2,dim" + std::to_string(numberOfDimensions) +".mat");

    // TODO Adjust objectiveValueTranslation
    compositeGriewankRosenbrockFunctionF8F2.setObjectiveValueTranslation(0);
    compositeGriewankRosenbrockFunctionF8F2.setRotationR(rotationR);

    for (std::size_t n = 0; n < parameters.n_cols; ++n) {
      CHECK(compositeGriewankRosenbrockFunctionF8F2.getObjectiveValue(parameters.col(n)) == Approx(expected.at(n)));
    }
  }
}

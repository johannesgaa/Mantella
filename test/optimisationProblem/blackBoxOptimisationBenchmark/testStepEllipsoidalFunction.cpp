// Catch
#include <catch.hpp>

// C++ Standard Library
#include <cstdlib>
#include <string>

// Armadillo
#include <armadillo>

// Mantella
#include <mantella>

extern std::string testDirectory;

TEST_CASE("bbob::StepEllipsoidalFunction", "") {
  for (const auto& numberOfDimensions : {2, 40}) {
    mant::bbob::StepEllipsoidalFunction stepEllipsoidalFunction(numberOfDimensions);

    arma::Mat<double> parameters;
    parameters.load(testDirectory + "/data/optimisationProblem/blackBoxOptimisationBenchmark/parameters,dim" + std::to_string(numberOfDimensions) +".mat");

    arma::Col<double> translation;
    translation.load(testDirectory + "/data/optimisationProblem/blackBoxOptimisationBenchmark/translation,dim" + std::to_string(numberOfDimensions) +".mat");

    arma::Mat<double> rotationR;
    rotationR.load(testDirectory + "/data/optimisationProblem/blackBoxOptimisationBenchmark/rotationR,dim" + std::to_string(numberOfDimensions) +".mat");

    arma::Mat<double> rotationQ;
    rotationQ.load(testDirectory + "/data/optimisationProblem/blackBoxOptimisationBenchmark/rotationQ,dim" + std::to_string(numberOfDimensions) +".mat");

    arma::Col<double> expected;
    expected.load(testDirectory + "/data/optimisationProblem/blackBoxOptimisationBenchmark/expectedStepEllipsoidalFunction,dim" + std::to_string(numberOfDimensions) +".mat");

    stepEllipsoidalFunction.setObjectiveValueTranslation(0);
    stepEllipsoidalFunction.setParameterTranslation(translation);
    stepEllipsoidalFunction.setParameterRotationR(rotationR);
    stepEllipsoidalFunction.setParameterRotationQ(rotationQ);

    for (std::size_t n = 0; n < parameters.n_cols; ++n) {
      CHECK(stepEllipsoidalFunction.getObjectiveValue(parameters.col(n)) == Approx(expected.at(n)));
    }
  }

  SECTION("Returns the specified class name.") {
    CHECK(mant::bbob::StepEllipsoidalFunction(5).toString() == "bbob_step_ellipsoidal_function");
  }
}

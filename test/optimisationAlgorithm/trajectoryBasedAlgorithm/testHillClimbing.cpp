// Catch
#include <catch.hpp>

// C++ Standard Library
#include <memory>
#include <iostream>

// Armadillo
#include <armadillo>

// HOP
#include <hop>

class TestHillClimbing : public hop::HillClimbing {
  public:
    TestHillClimbing(
        const std::shared_ptr<hop::OptimisationProblem> optimisationProblem)
      : HillClimbing(optimisationProblem),
        velocityIndex_(0){
      velocities_.load("/Users/SRA/Documents/workspace/OnlineOptimisation/test/data/testVel.mat");
    }

    arma::Col<double> getVelocity() {
      return velocities_.col(velocityIndex_++);
    }

  protected:
    unsigned int velocityIndex_;
    arma::mat velocities_;
};

class TestHillClimbingProblem : public hop::OptimisationProblem {
  public:
    TestHillClimbingProblem(
        const unsigned int numberOfDimensions)
      : OptimisationProblem(numberOfDimensions),
        objectiveValueIndex_(0) {
      objectiveValues_.load("/Users/SRA/Documents/workspace/OnlineOptimisation/test/data/testObj.mat");
    }

    arma::Mat<double> getParameterHistory() const noexcept {
      return parameterHistory_;
    }

  protected:
    unsigned int objectiveValueIndex_;
    arma::Col<double> objectiveValues_;
    int n;

    arma::Mat<double> parameterHistory_;

    double getObjectiveValueImplementation(
       const arma::Col<double>& parameter) const override {
       //parameterHistory_.insert_cols(n, parameter);

      return objectiveValues_.at(objectiveValueIndex_);
    }

    std::string to_string() const noexcept {
      return "HillClimbing";
    }
};

TEST_CASE("Hill climbing", "") {
  std::shared_ptr<TestHillClimbingProblem> testHillClimbingProblem(new TestHillClimbingProblem(4));

  TestHillClimbing testHillClimbing(testHillClimbingProblem);
  testHillClimbing.setInitialParameter(arma::zeros<arma::Col<double>>(testHillClimbingProblem->getNumberOfDimensions()));
  testHillClimbing.setMaximalNumberOfIterations(4);

  testHillClimbing.optimise();
  arma::Mat<double> actualParameterHistory = testHillClimbingProblem->getParameterHistory();
  actualParameterHistory.load("/Users/SRA/Documents/workspace/OnlineOptimisation/test/data/testExp.mat");
  arma::Mat<double> expectedParameterHistory;
  expectedParameterHistory.load("/Users/SRA/Documents/workspace/OnlineOptimisation/test/data/testExp.mat");

  for(std::size_t n = 0; n < expectedParameterHistory.n_cols; ++n) {
    arma::Col<double> expectedParameter = expectedParameterHistory.col(n);
    arma::Col<double> actualParameter = actualParameterHistory.col(n);

    for (std::size_t k = 0; k < expectedParameter.n_elem; ++k) {
      CHECK(actualParameter.at(k) == Approx(expectedParameter.at(k)));
    }
  }
}


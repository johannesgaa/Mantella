// Catch
#include <catch.hpp>
#include <helper.hpp>

// C++ Standard Library
#include <memory>
#include <iostream>
#include <cmath>

// Armadillo
#include <armadillo>

// Mantella
#include <mantella>

extern boost::filesystem::path testDirectory;
static std::string dataPath_ = "/data/optimisationAlgorithm/trajectoryBasedAlgrorithm/hillClimbing/";


class TestHillClimbing : public mant::HillClimbing<double, mant::EuclideanDistance> {
  public:
    TestHillClimbing(
        const std::shared_ptr<mant::OptimisationProblem<double>> optimisationProblem)
      : mant::HillClimbing<double, mant::EuclideanDistance>(optimisationProblem),
        velocityIndex_(0){
        
    }

    arma::Col<double> getVelocity() noexcept {
      return velocities_.col(velocityIndex_++);
    }

    void setVelocitys(arma::mat velocities){
      velocities_ = velocities;
    }

  protected:
    unsigned int velocityIndex_;
    arma::mat velocities_;
};

class TestHillClimbingProblem : public mant::OptimisationProblem<double> {
  public:
    TestHillClimbingProblem(
        const unsigned int numberOfDimensions)
      : mant::OptimisationProblem<double>(numberOfDimensions) {
      softConstraintsValuesIndex_ = 0;
      objectiveValuesIndex_ = 0;
      parameterHistory_.clear();
    }


    std::vector<arma::Col<double>> getParameterHistory() const noexcept {
      return parameterHistory_;
    }

    void setObjectiveValues(arma::Col<double> objectiveValues) {
      objectiveValues_ = objectiveValues;
    }

    void setSoftConstraintsValues(arma::Col<double> softConstraintsValues) {
      softConstraintsValues_ = softConstraintsValues;
    }

  protected:
    static unsigned int softConstraintsValuesIndex_;
    static unsigned int objectiveValuesIndex_;
    arma::Col<double> objectiveValues_;
    arma::Col<double> softConstraintsValues_;

    static std::vector<arma::Col<double>> parameterHistory_;

    double getSoftConstraintsValueImplementation(
       const arma::Col<double>& parameter) const noexcept override {

      return softConstraintsValues_.at(softConstraintsValuesIndex_++);
    }

    double getObjectiveValueImplementation(
       const arma::Col<double>& parameter) const noexcept override {
       parameterHistory_.push_back(parameter);

      return objectiveValues_.at(objectiveValuesIndex_++);
    }

    std::string to_string() const noexcept {
      return "TestHillClimbing";
    }
};

decltype(TestHillClimbingProblem::parameterHistory_) TestHillClimbingProblem::parameterHistory_;
decltype(TestHillClimbingProblem::softConstraintsValuesIndex_) TestHillClimbingProblem::softConstraintsValuesIndex_;
decltype(TestHillClimbingProblem::objectiveValuesIndex_) TestHillClimbingProblem::objectiveValuesIndex_;



TEST_CASE("Hill climbing", "") {
   // Set OptimisationProblem values
   arma::Col<double> upperBounds;
   upperBounds.load(testDirectory.string() + dataPath_ + "upperBounds[2.1].mat"); // The parameter is optional

   arma::Col<double> lowerBounds;
   lowerBounds.load(testDirectory.string() + dataPath_ + "lowerBounds[2.1].mat"); // The parameter is optional

   arma::Col<double> objectiveValues;
   objectiveValues.load(testDirectory.string() + dataPath_ + "objectiveValues[2.1].mat");

   arma::Col<double> softConstraintsValues;
   softConstraintsValues.load(testDirectory.string() + dataPath_ + "softConstraintsValues[2.1].mat");

   // Init OptimisationProblem
   std::shared_ptr<TestHillClimbingProblem> testHillClimbingProblem(new TestHillClimbingProblem(upperBounds.n_elem));
   testHillClimbingProblem->setUpperBounds(upperBounds);
   testHillClimbingProblem->setLowerBounds(lowerBounds);
   testHillClimbingProblem->setAcceptableObjectiveValue(-500.0);
   testHillClimbingProblem->setObjectiveValues(objectiveValues);
   testHillClimbingProblem->setSoftConstraintsValues(softConstraintsValues);

   // Set OptimisationAlgorithm values
   arma::mat velocities;
   velocities.load(testDirectory.string() + dataPath_ + "velocities[2.1].mat");

   arma::Col<double> initialParameter;
   initialParameter.load(testDirectory.string() + dataPath_ + "initialParameter[2.1].mat"); // The parameter is optional

   // Init OptimisationAlgorithm
   TestHillClimbing testHillClimbing(testHillClimbingProblem);
   testHillClimbing.setVelocitys(velocities);
   testHillClimbing.setInitialParameter(initialParameter);
   testHillClimbing.setMaximalStepSize(2.0);

   // Set Expected values
   arma::Mat<double> expectedParameterHistory;
   expectedParameterHistory.load(testDirectory.string() + dataPath_ + "expected[2.1].mat");

   // Run hill climbing algorithm
  testHillClimbing.optimise();
  // Comparing of candidateParameters
  std::vector<arma::Col<double>> actualParameterHistory = testHillClimbingProblem->getParameterHistory();
  compare(actualParameterHistory,expectedParameterHistory);
}

TEST_CASE("HillClimbing Test Maximalstepsize", "") {
  // name for the expected data
  std::string expectedName = "1.1";

  //Set OptimisationProblem values
  arma::Col<double> upperBounds;
  upperBounds.load(testDirectory.string() + dataPath_ + "upperBounds["+ expectedName +"].mat");

  arma::Col<double> lowerBounds;
  lowerBounds.load(testDirectory.string() + dataPath_ + "lowerBounds["+ expectedName +"].mat");

  arma::Col<double> objectiveValues;
  objectiveValues.load(testDirectory.string() + dataPath_ + "objectiveValues["+ expectedName +"].mat");

  arma::Col<double> softConstraintsValues;
  softConstraintsValues.load(testDirectory.string() + dataPath_ + "softConstraintsValues["+ expectedName +"].mat");

  // Init OptimisationProblem
  std::shared_ptr<TestHillClimbingProblem> testHillClimbingProblem(new TestHillClimbingProblem(upperBounds.n_elem));
  testHillClimbingProblem->setUpperBounds(upperBounds);
  testHillClimbingProblem->setLowerBounds(lowerBounds);
  testHillClimbingProblem->setObjectiveValues(objectiveValues);
  testHillClimbingProblem->setSoftConstraintsValues(softConstraintsValues);
  testHillClimbingProblem->setAcceptableObjectiveValue(10.0);

  // Set OptimisationAlgorithm values
  arma::mat velocities;
  velocities.load(testDirectory.string() + dataPath_ + "velocities["+ expectedName +"].mat");

  arma::Col<double> initialParameter;
  initialParameter.load(testDirectory.string() + dataPath_ + "initialParameter["+ expectedName +"].mat");

  // Init OptimisationAlgorithm
  TestHillClimbing testHillClimbing(testHillClimbingProblem);
  testHillClimbing.setVelocitys(velocities);
  testHillClimbing.setInitialParameter(initialParameter);

  // Set Expected values
  arma::Mat<double> expectedParameterHistory;
  expectedParameterHistory.load(testDirectory.string() + dataPath_ + "expected["+ expectedName +"].mat");

  // Run hill climbing algorithm
  testHillClimbing.optimise();

  // Comparing of candidateParameters
  std::vector<arma::Col<double>> actualParameterHistory = testHillClimbingProblem->getParameterHistory();
  compare(actualParameterHistory,expectedParameterHistory);
}

TEST_CASE("HillClimbing Check initialParameter at each limit", "") {
   // name for the expected data
   std::string expectedName = "1.2";

   //Set OptimisationProblem values
   arma::Col<double> upperBounds;
   upperBounds.load(testDirectory.string() + dataPath_ + "upperBounds["+ expectedName +"].mat");

   arma::Col<double> lowerBounds;
   lowerBounds.load(testDirectory.string() + dataPath_ + "lowerBounds["+ expectedName +"].mat");

   arma::Col<double> objectiveValues;
   objectiveValues.load(testDirectory.string() + dataPath_ + "objectiveValues["+ expectedName +"].mat");

   arma::Col<double> softConstraintsValues;
   softConstraintsValues.load(testDirectory.string() + dataPath_ + "softConstraintsValues["+ expectedName +"].mat");

   // Init OptimisationProblem
   std::shared_ptr<TestHillClimbingProblem> testHillClimbingProblem(new TestHillClimbingProblem(upperBounds.n_elem));
   testHillClimbingProblem->setUpperBounds(upperBounds);
   testHillClimbingProblem->setLowerBounds(lowerBounds);
   testHillClimbingProblem->setObjectiveValues(objectiveValues);
   testHillClimbingProblem->setSoftConstraintsValues(softConstraintsValues);
   testHillClimbingProblem->setAcceptableObjectiveValue(10.0);

   // Set OptimisationAlgorithm values
   arma::mat velocities;
   velocities.load(testDirectory.string() + dataPath_ + "velocities["+ expectedName +"].mat");

   arma::Col<double> initialParameter;
   initialParameter.load(testDirectory.string() + dataPath_ + "initialParameter["+ expectedName +"].mat");

   // Init OptimisationAlgorithm
   TestHillClimbing testHillClimbing(testHillClimbingProblem);
   testHillClimbing.setVelocitys(velocities);
   testHillClimbing.setInitialParameter(initialParameter);

   // Set Expected values
   arma::Mat<double> expectedParameterHistory;
   expectedParameterHistory.load(testDirectory.string() + dataPath_ + "expected["+ expectedName +"].mat");

   // Run hill climbing algorithm
   testHillClimbing.optimise();

   // Comparing of candidateParameters
   std::vector<arma::Col<double>> actualParameterHistory = testHillClimbingProblem->getParameterHistory();
   compare(actualParameterHistory,expectedParameterHistory);
}

TEST_CASE("HillClimbing Check initialParameter at one limit", "") {
  // name for the expected data
  std::string expectedName = "1.3";

  //Set OptimisationProblem values
  arma::Col<double> upperBounds;
  upperBounds.load(testDirectory.string() + dataPath_ + "upperBounds["+ expectedName +"].mat");

  arma::Col<double> lowerBounds;
  lowerBounds.load(testDirectory.string() + dataPath_ + "lowerBounds["+ expectedName +"].mat");

  arma::Col<double> objectiveValues;
  objectiveValues.load(testDirectory.string() + dataPath_ + "objectiveValues["+ expectedName +"].mat");

  arma::Col<double> softConstraintsValues;
  softConstraintsValues.load(testDirectory.string() + dataPath_ + "softConstraintsValues["+ expectedName +"].mat");

  // Init OptimisationProblem
  std::shared_ptr<TestHillClimbingProblem> testHillClimbingProblem(new TestHillClimbingProblem(upperBounds.n_elem));
  testHillClimbingProblem->setUpperBounds(upperBounds);
  testHillClimbingProblem->setLowerBounds(lowerBounds);
  testHillClimbingProblem->setObjectiveValues(objectiveValues);
  testHillClimbingProblem->setSoftConstraintsValues(softConstraintsValues);
  testHillClimbingProblem->setAcceptableObjectiveValue(10.0);

  // Set OptimisationAlgorithm values
  arma::mat velocities;
  velocities.load(testDirectory.string() + dataPath_ + "velocities["+ expectedName +"].mat");

  arma::Col<double> initialParameter;
  initialParameter.load(testDirectory.string() + dataPath_ + "initialParameter["+ expectedName +"].mat");

  // Init OptimisationAlgorithm
  TestHillClimbing testHillClimbing(testHillClimbingProblem);
  testHillClimbing.setVelocitys(velocities);
  testHillClimbing.setInitialParameter(initialParameter);

  // Set Expected values
  arma::Mat<double> expectedParameterHistory;
  expectedParameterHistory.load(testDirectory.string() + dataPath_ + "expected["+ expectedName +"].mat");

  // Run hill climbing algorithm
  testHillClimbing.optimise();

  // Comparing of candidateParameters
  std::vector<arma::Col<double>> actualParameterHistory = testHillClimbingProblem->getParameterHistory();
  compare(actualParameterHistory,expectedParameterHistory);
}

TEST_CASE("HillClimbing Check initialParameter in-range", "") {
  // name for the expected data
  std::string expectedName = "1.4";

  //Set OptimisationProblem values
  arma::Col<double> upperBounds;
  upperBounds.load(testDirectory.string() + dataPath_ + "upperBounds["+ expectedName +"].mat");

  arma::Col<double> lowerBounds;
  lowerBounds.load(testDirectory.string() + dataPath_ + "lowerBounds["+ expectedName +"].mat");

  arma::Col<double> objectiveValues;
  objectiveValues.load(testDirectory.string() + dataPath_ + "objectiveValues["+ expectedName +"].mat");

  arma::Col<double> softConstraintsValues;
  softConstraintsValues.load(testDirectory.string() + dataPath_ + "softConstraintsValues["+ expectedName +"].mat");

  // Init OptimisationProblem
  std::shared_ptr<TestHillClimbingProblem> testHillClimbingProblem(new TestHillClimbingProblem(upperBounds.n_elem));
  testHillClimbingProblem->setUpperBounds(upperBounds);
  testHillClimbingProblem->setLowerBounds(lowerBounds);
  testHillClimbingProblem->setObjectiveValues(objectiveValues);
  testHillClimbingProblem->setSoftConstraintsValues(softConstraintsValues);
  testHillClimbingProblem->setAcceptableObjectiveValue(10.0);

  // Set OptimisationAlgorithm values
  arma::mat velocities;
  velocities.load(testDirectory.string() + dataPath_ + "velocities["+ expectedName +"].mat");

  arma::Col<double> initialParameter;
  initialParameter.load(testDirectory.string() + dataPath_ + "initialParameter["+ expectedName +"].mat");

  // Init OptimisationAlgorithm
  TestHillClimbing testHillClimbing(testHillClimbingProblem);
  testHillClimbing.setVelocitys(velocities);
  testHillClimbing.setInitialParameter(initialParameter);

  // Set Expected values
  arma::Mat<double> expectedParameterHistory;
  expectedParameterHistory.load(testDirectory.string() + dataPath_ + "expected["+ expectedName +"].mat");

  // Run hill climbing algorithm
  testHillClimbing.optimise();

  // Comparing of candidateParameters
  std::vector<arma::Col<double>> actualParameterHistory = testHillClimbingProblem->getParameterHistory();
  compare(actualParameterHistory,expectedParameterHistory);
}

TEST_CASE("HillClimbing Check initialParameter out of range maximal limit", "") {
  // name for the expected data
  std::string expectedName = "1.5";

  //Set OptimisationProblem values
  arma::Col<double> upperBounds;
  upperBounds.load(testDirectory.string() + dataPath_ + "upperBounds["+ expectedName +"].mat");

  arma::Col<double> lowerBounds;
  lowerBounds.load(testDirectory.string() + dataPath_ + "lowerBounds["+ expectedName +"].mat");

  arma::Col<double> objectiveValues;
  objectiveValues.load(testDirectory.string() + dataPath_ + "objectiveValues["+ expectedName +"].mat");

  arma::Col<double> softConstraintsValues;
  softConstraintsValues.load(testDirectory.string() + dataPath_ + "softConstraintsValues["+ expectedName +"].mat");

  // Init OptimisationProblem
  std::shared_ptr<TestHillClimbingProblem> testHillClimbingProblem(new TestHillClimbingProblem(upperBounds.n_elem));
  testHillClimbingProblem->setUpperBounds(upperBounds);
  testHillClimbingProblem->setLowerBounds(lowerBounds);
  testHillClimbingProblem->setObjectiveValues(objectiveValues);
  testHillClimbingProblem->setSoftConstraintsValues(softConstraintsValues);
  testHillClimbingProblem->setAcceptableObjectiveValue(10.0);

  // Set OptimisationAlgorithm values
  arma::mat velocities;
  velocities.load(testDirectory.string() + dataPath_ + "velocities["+ expectedName +"].mat");

  arma::Col<double> initialParameter;
  initialParameter.load(testDirectory.string() + dataPath_ + "initialParameter["+ expectedName +"].mat");

  // Init OptimisationAlgorithm
  TestHillClimbing testHillClimbing(testHillClimbingProblem);
  testHillClimbing.setVelocitys(velocities);
  testHillClimbing.setInitialParameter(initialParameter);

  // Set Expected values
  arma::Mat<double> expectedParameterHistory;
  expectedParameterHistory.load(testDirectory.string() + dataPath_ + "expected["+ expectedName +"].mat");

  // Run hill climbing algorithm
  testHillClimbing.optimise();

  // Comparing of candidateParameters
  std::vector<arma::Col<double>> actualParameterHistory = testHillClimbingProblem->getParameterHistory();
  compare(actualParameterHistory,expectedParameterHistory);
}

TEST_CASE("HillClimbing Check initialParameter out of range arbitrary values", "") {
  // name for the expected data
  std::string expectedName = "1.6";

  //Set OptimisationProblem values
  arma::Col<double> upperBounds;
  upperBounds.load(testDirectory.string() + dataPath_ + "upperBounds["+ expectedName +"].mat");

  arma::Col<double> lowerBounds;
  lowerBounds.load(testDirectory.string() + dataPath_ + "lowerBounds["+ expectedName +"].mat");

  arma::Col<double> objectiveValues;
  objectiveValues.load(testDirectory.string() + dataPath_ + "objectiveValues["+ expectedName +"].mat");

  arma::Col<double> softConstraintsValues;
  softConstraintsValues.load(testDirectory.string() + dataPath_ + "softConstraintsValues["+ expectedName +"].mat");

  // Init OptimisationProblem
  std::shared_ptr<TestHillClimbingProblem> testHillClimbingProblem(new TestHillClimbingProblem(upperBounds.n_elem));
  testHillClimbingProblem->setUpperBounds(upperBounds);
  testHillClimbingProblem->setLowerBounds(lowerBounds);
  testHillClimbingProblem->setObjectiveValues(objectiveValues);
  testHillClimbingProblem->setSoftConstraintsValues(softConstraintsValues);
  testHillClimbingProblem->setAcceptableObjectiveValue(10.0);

  // Set OptimisationAlgorithm values
  arma::mat velocities;
  velocities.load(testDirectory.string() + dataPath_ + "velocities["+ expectedName +"].mat");

  arma::Col<double> initialParameter;
  initialParameter.load(testDirectory.string() + dataPath_ + "initialParameter["+ expectedName +"].mat");

  // Init OptimisationAlgorithm
  TestHillClimbing testHillClimbing(testHillClimbingProblem);
  testHillClimbing.setVelocitys(velocities);
  testHillClimbing.setInitialParameter(initialParameter);

  // Set Expected values
  arma::Mat<double> expectedParameterHistory;
  expectedParameterHistory.load(testDirectory.string() + dataPath_ + "expected["+ expectedName +"].mat");

  // Run hill climbing algorithm
  testHillClimbing.optimise();

  // Comparing of candidateParameters
  std::vector<arma::Col<double>> actualParameterHistory = testHillClimbingProblem->getParameterHistory();
  compare(actualParameterHistory,expectedParameterHistory);
}
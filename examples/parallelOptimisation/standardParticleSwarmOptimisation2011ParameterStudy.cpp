// C++ Standard Library
#include <cstdlib>
#include <array>
#include <stdexcept>

// Armadillo
#include <armadillo>

// Boost
#include <boost/filesystem.hpp>

// Cereal
#include <cereal/archives/json.hpp>
#include <cereal/types/polymorphic.hpp>

// HOP
#include <hop>

int main (const int argc, const char* argv[]) {
  try {
    arma::arma_rng::set_seed_random();

    unsigned int optimisationProblemIndex = std::stoi(argv[1]);
    unsigned int numberOfDimensions = std::stoi(argv[2]);
    unsigned int populationSize = std::stoi(argv[3]);

    std::array<std::shared_ptr<hop::OptimisationProblem>, 24> optimisationProblems({
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::SphereFunction(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::EllipsoidalFunction(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::RastriginFunction(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::BuecheRastriginFunction(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::LinearSlope(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::AttractiveSectorFunction(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::StepEllipsoidalFunction(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::RosenbrockFunction(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::RosenbrockFunctionRotated(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::EllipsoidalFunctionRotated(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::DiscusFunction(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::BentCigarFunction(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::SharpRidgeFunction(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::DifferentPowersFunction(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::RastriginFunctionRotated(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::WeierstrassFunction(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::SchaffersF7Function(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::SchaffersF7FunctionIllConditioned(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::CompositeGriewankRosenbrockFunctionF8F2(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::SchwefelFunction(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::GallaghersGaussian101mePeaksFunction(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::GallaghersGaussian21hiPeaksFunction(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::KatsuuraFunction(numberOfDimensions)),
      std::shared_ptr<hop::OptimisationProblem>(new hop::bbob2013::LunacekBiRastriginFunction(numberOfDimensions)),
    });

    std::shared_ptr<hop::OptimisationProblem> optimisationProblem = optimisationProblems.at(optimisationProblemIndex);

    arma::Col<double> localAttractionParameters = arma::linspace(0, 4, 10);
    arma::Col<double> globalAttractionParameters = arma::linspace(0, 4, 10);
    arma::Col<double> accelerationParameters = arma::linspace(0, 4, 10);
    arma::Col<double> neighbourhoodParameters = arma::linspace(0.1, 1, 5);

    unsigned int iterations = 30;

    hop::StandardParticleSwarmOptimisation2011 optimisationAlgorithm(optimisationProblem, populationSize);
    optimisationAlgorithm.setMaximalNumberOfIterations(10000);

    for(auto localAttractionParameter : localAttractionParameters) {
      optimisationAlgorithm.setLocalAttraction(localAttractionParameter);

      for(auto globalAttractionParameter : globalAttractionParameters) {
        optimisationAlgorithm.setGlobalAttraction(globalAttractionParameter);

        for(auto accelerationParameter: accelerationParameters) {
          optimisationAlgorithm.setAcceleration(accelerationParameter);

          for(auto neighbourhoodParameter : neighbourhoodParameters) {
            optimisationAlgorithm.setNeighbourhoodProbability(neighbourhoodParameter);

                      boost::filesystem::path filepath = boost::filesystem::path("./parameterisation_prob" + hop::to_string(optimisationProblem) + "_dim" + std::to_string(numberOfDimensions)+ "_pop" + std::to_string(populationSize) + "_local" + std::to_string(localAttractionParameter) + "_global" + std::to_string(globalAttractionParameter) + "_acc" + std::to_string(accelerationParameter) + "_neigh" + std::to_string(neighbourhoodParameter) + ".mat");

             if(!boost::filesystem::exists(filepath)) {
              arma::Mat<double> results(iterations, 6 + numberOfDimensions);
              for(std::size_t n = 0; n < iterations; ++n) {
                optimisationAlgorithm.optimise();

                results.at(n, 0) = optimisationProblem->getNumberOfEvaluations();
                results.at(n, 1) = optimisationProblem->getNumberOfDistinctEvaluations();
                results.at(n, 2) = optimisationAlgorithm.getNumberOfIterations();
                results.at(n, 3) = optimisationAlgorithm.getBestObjectiveValue() - optimisationProblem->getAcceptableObjectiveValue();
                results.at(n, 4) = optimisationAlgorithm.isFinished();
                results.at(n, 5) = optimisationAlgorithm.isTerminated();
                results(n, arma::span(6, 6 + numberOfDimensions - 1)) = optimisationAlgorithm.getBestParameter().t();
              }

              results.save(filepath.string(), arma::raw_ascii);
            }
          }
        }
      }
    }
    return EXIT_SUCCESS;
  } catch(const std::exception& exception) {
    std::cout << exception.what();
  }
}
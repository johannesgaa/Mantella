// Catch// Catch
#include <catch.hpp>

// C++ standard library
#include <random>
#include <memory>

// Armadillo
#include <armadillo>

// Mantella
#include <mantella>

TEST_CASE("BestSamplesSelection") {
  SECTION(".select") {
  
  }

  SECTION(".toString") {
    SECTION("Returns the expected class name.") {
      mant::cacheSamples = true;
      std::shared_ptr<mant::OptimisationProblem> optimisationProblem(new mant::bbob::SphereFunction(std::uniform_int_distribution<arma::uword>(1, 10)(mant::Rng::getGenerator())));
      mant::RandomSearch randomSearch(optimisationProblem);
      randomSearch.optimise();
      CHECK(mant::BestSamplesSelection(optimisationProblem->getCachedSamples(), std::uniform_int_distribution<arma::uword>(1, optimisationProblem->getCachedSamples().size())(mant::Rng::getGenerator())).toString() == "best_samples_selection");
    }
  }
}
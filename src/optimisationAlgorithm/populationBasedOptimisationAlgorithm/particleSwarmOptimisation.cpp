#include "mantella_bits/optimisationAlgorithm/populationBasedOptimisationAlgorithm/particleSwarmOptimisation.hpp"

// C++ standard library
#include <cmath>
#include <functional>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <vector>

// Mantella
#include "mantella_bits/optimisationProblem.hpp"
#include "mantella_bits/probability.hpp"

namespace mant {
  ParticleSwarmOptimisation::ParticleSwarmOptimisation()
      : PopulationBasedOptimisationAlgorithm() {
#if defined(SUPPORT_MPI)
    MPI_Op_create(&mpiOpBestSample, 1, &MPI_Op_best_sample);

    setCommunicationFunctions(
        {{[this](
              const arma::uword numberOfDimensions_,
              const arma::mat& parameters_,
              const arma::rowvec& objectiveValues_,
              const arma::rowvec& differences_) {
            MPI_Datatype MANT_MPI_SAMPLE;
            MPI_Type_contiguous(static_cast<int>(2 + numberOfDimensions_), MPI_DOUBLE, &MANT_MPI_SAMPLE);
            MPI_Type_commit(&MANT_MPI_SAMPLE);
            arma::vec mpiSendSample(2 + numberOfDimensions_);
            mpiSendSample(0) = static_cast<double>(numberOfDimensions_);
            mpiSendSample(1) = bestFoundObjectiveValue_;
            mpiSendSample.tail(numberOfDimensions_) = bestFoundParameter_;

            arma::vec mpiReceiveSample(arma::size(mpiSendSample));
            MPI_Allreduce(mpiSendSample.memptr(), mpiReceiveSample.memptr(), 1, MANT_MPI_SAMPLE, MPI_Op_best_sample, MPI_COMM_WORLD);

            bestFoundObjectiveValue_ = mpiReceiveSample(1);
            bestFoundParameter_ = mpiReceiveSample.tail(numberOfDimensions_);
            return parameters_;
        }, "Find the cluster-wide best found sample"}});
#endif
        
    setInitialisingFunctions(
        {{[this](
              const arma::uword numberOfDimensions_,
              const arma::mat& initialParameters_) {
            velocities_ = uniformRandomNumbers(numberOfDimensions_, populationSize_, std::uniform_real_distribution<double>(-1.0, 1.0)) - initialParameters_;

            return initialParameters_;
          },
          "Velocity initialisation"},
         {[this](
              const arma::uword numberOfDimensions_,
              const arma::mat& initialParameters_) {
            localBestSolutions_ = initialParameters_;

            return initialParameters_;
          },
          "Local best solution initialisation"},
         {[this](
              const arma::uword numberOfDimensions_,
              const arma::mat& initialParameters_) {
            localBestObjectiveValues_.set_size(populationSize_);
            localBestObjectiveValues_.fill(std::numeric_limits<double>::infinity());

            return initialParameters_;
          },
          "Local best objective values initialisation"}});

    setNextParametersFunctions(
        {{[this](
              const arma::uword numberOfDimensions_,
              const arma::mat& parameters_,
              const arma::rowvec& objectiveValues_,
              const arma::rowvec& differences_) {
            for (arma::uword n = 0; n < populationSize_; ++n) {
              if (objectiveValues_(n) < localBestObjectiveValues_(n)) {
                localBestSolutions_.col(n) = parameters_.col(n);
                localBestObjectiveValues_(n) = objectiveValues_(n);
              }
            }
            
            arma::mat particles(arma::size(parameters_));

#pragma omp parallel
            {
#pragma omp for schedule(static)
            for (arma::uword n = 0; n < populationSize_; ++n) {
              const arma::vec& particle = parameters_.col(n);

              arma::vec attractionCenter = (getMaximalLocalAttraction() * (localBestSolutions_.col(n) - particle) + getMaximalGlobalAttraction() * (getBestFoundParameter() - particle)) / 3.0;
              const arma::vec& velocity = getMaximalAcceleration() * uniformRandomNumbers(numberOfDimensions_) % velocities_.col(n) + randomNeighbour(attractionCenter, 0, arma::norm(attractionCenter));

              particles.col(n) = particle + velocity;
              velocities_.col(n) = velocity;
            }
            }

            return particles;
          },
          "Particle swarm optimisation"}});

    auto boundariesHandlingFunctions = getBoundariesHandlingFunctions();
    boundariesHandlingFunctions.push_back(
        {[this](
             const arma::mat& parameters_,
             const arma::umat& isBelowLowerBound_,
             const arma::umat& isAboveUpperBound_) {
           velocities_.elem(arma::find(isBelowLowerBound_)) *= -0.5;
           velocities_.elem(arma::find(isAboveUpperBound_)) *= -0.5;

           return parameters_;
         },
         "Halve velocity and revert direction"});
    setBoundariesHandlingFunctions(boundariesHandlingFunctions);

    setMaximalAcceleration(1.0 / (2.0 * std::log(2.0)));
    setMaximalLocalAttraction(0.5 + std::log(2.0));
    setMaximalGlobalAttraction(maximalLocalAttraction_);
  }

  void ParticleSwarmOptimisation::optimise(
      OptimisationProblem& optimisationProblem) {
    optimise(optimisationProblem, uniformRandomNumbers(optimisationProblem.numberOfDimensions_, populationSize_));
  }

  void ParticleSwarmOptimisation::setMaximalAcceleration(
      const double maximalAcceleration) {
    maximalAcceleration_ = maximalAcceleration;
  }

  double ParticleSwarmOptimisation::getMaximalAcceleration() const {
    return maximalAcceleration_;
  }

  void ParticleSwarmOptimisation::setMaximalLocalAttraction(
      const double maximalLocalAttraction) {
    maximalLocalAttraction_ = maximalLocalAttraction;
  }

  double ParticleSwarmOptimisation::getMaximalLocalAttraction() const {
    return maximalLocalAttraction_;
  }

  void ParticleSwarmOptimisation::setMaximalGlobalAttraction(
      const double maximalGlobalAttraction) {
    maximalGlobalAttraction_ = maximalGlobalAttraction;
  }

  double ParticleSwarmOptimisation::getMaximalGlobalAttraction() const {
    return maximalGlobalAttraction_;
  }
}

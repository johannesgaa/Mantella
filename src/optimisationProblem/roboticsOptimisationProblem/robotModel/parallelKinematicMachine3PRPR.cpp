#include <mantella_bits/optimisationProblem/roboticsOptimisationProblem/robotModel/parallelKinematicMachine3PRPR.hpp>

// C++ standard library
#include <cassert>

// Mantella
#include <mantella_bits/helper/assert.hpp>
#include <mantella_bits/helper/geometry.hpp>

namespace mant {
  namespace robotics {
    ParallelKinematicMachine3PRPR::ParallelKinematicMachine3PRPR(
        const arma::uword numberOfRedundantJoints) 
      : RobotModel(3, numberOfRedundantJoints) {
      setMinimalActiveJointsActuation({0.1, 0.1, 0.1});
      setMaximalActiveJointsActuation({1.2, 1.2, 1.2});

      setEndEffectorJointPositions({
        -0.000066580445834, 0.106954081945581,
        -0.092751709777083, -0.053477040972790,
        0.092818290222917, -0.053477040972790});

      setRedundantJointStartPositions({
        0.1, 1.0392,
        0.0, 0.8,
        1.2, 0.8});

      setRedundantJointEndPositions({
        1.1, 1.0392,
        0.0, -0.2,
        1.2, -0.2});

      redundantJointStartToEndPositions_ = redundantJointEndPositions_ - redundantJointStartPositions_;
      redundantJointIndicies_ = arma::find(arma::any(redundantJointStartToEndPositions_));

      redundantJointAngleSines_.set_size(redundantJointIndicies_.n_elem);
      redundantJointAngleCosines_.set_size(redundantJointIndicies_.n_elem);

      for (arma::uword n = 0; n < redundantJointIndicies_.n_elem; ++n) {
        const double redundantJointAngle = std::atan2(redundantJointStartToEndPositions_(1, n), redundantJointStartToEndPositions_(0, n));
        redundantJointAngleSines_(n) = std::sin(redundantJointAngle);
        redundantJointAngleCosines_(n) = std::cos(redundantJointAngle);
      }
    }
            
    void ParallelKinematicMachine3PRPR::setEndEffectorJointPositions(
        arma::Mat<double>::fixed<2, 3> endEffectorJointPositions) {
      endEffectorJointPositions_ = endEffectorJointPositions;
    }
    
    void ParallelKinematicMachine3PRPR::setRedundantJointStartPositions(
        arma::Mat<double>::fixed<2, 3> redundantJointStartPositions) {
      redundantJointStartPositions_ = redundantJointStartPositions;
    }
    
    void ParallelKinematicMachine3PRPR::setRedundantJointEndPositions(
        arma::Mat<double>::fixed<2, 3> redundantJointEndPositions) {
      redundantJointEndPositions_ = redundantJointEndPositions;
    }
        
    arma::Mat<double>::fixed<2, 3> ParallelKinematicMachine3PRPR::getEndEffectorJointPositions() const {
      return endEffectorJointPositions_;
    }
    
    arma::Mat<double>::fixed<2, 3> ParallelKinematicMachine3PRPR::getRedundantJointStartPositions() const {
      return redundantJointStartPositions_;
    }
    
    arma::Mat<double>::fixed<2, 3> ParallelKinematicMachine3PRPR::getRedundantJointEndPositions() const {
      return redundantJointEndPositions_;
    }

    arma::Cube<double> ParallelKinematicMachine3PRPR::getModelImplementation(
        const arma::Col<double>& endEffectorPose,
        const arma::Row<double>& redundantJointsActuation) const {
      assert(redundantJointsActuation.n_elem == numberOfRedundantJoints_);
      assert(!arma::any(redundantJointsActuation < minimalRedundantJointsActuation_) && !arma::any(redundantJointsActuation > maximalRedundantJointsActuation_));

      arma::Cube<double> model;

      const arma::Col<double>& endEffectorPosition = endEffectorPose.subvec(0, 1);
      const double& endEffectorAngle = endEffectorPose(2);

      model.slice(0) = redundantJointStartPositions_;
      for (arma::uword n = 0; n < redundantJointIndicies_.n_elem; ++n) {
        const arma::uword& redundantJointIndex = redundantJointIndicies_(n);
        model.slice(0).col(redundantJointIndex) += redundantJointsActuation(redundantJointIndex) * redundantJointStartToEndPositions_.col(redundantJointIndex);
      }

      model.slice(1) = get2DRotation(endEffectorAngle) * endEffectorJointPositions_;
      model.slice(1).each_col() += endEffectorPosition;

      return model;
    }

    arma::Row<double> ParallelKinematicMachine3PRPR::getActuationImplementation(
        const arma::Col<double>& endEffectorPose,
        const arma::Row<double>& redundantJointsActuation) const {
      assert(redundantJointsActuation.n_elem == numberOfRedundantJoints_);
      assert(!arma::any(redundantJointsActuation < minimalRedundantJointsActuation_) && !arma::any(redundantJointsActuation > maximalRedundantJointsActuation_));
      
      const arma::Cube<double>& model = getModel(endEffectorPose, redundantJointsActuation);

      const arma::Mat<double>& baseJoints = model.slice(0);
      const arma::Mat<double>& endEffectorJoints = model.slice(1);

      return arma::sqrt(arma::sum(arma::square(endEffectorJoints - baseJoints)));
    }

    double ParallelKinematicMachine3PRPR::getEndEffectorPoseErrorImplementation(
        const arma::Col<double>& endEffectorPose,
        const arma::Row<double>& redundantJointsActuation) const {
      assert(redundantJointsActuation.n_elem == numberOfRedundantJoints_);
      assert(!arma::any(redundantJointsActuation < minimalRedundantJointsActuation_) && !arma::any(redundantJointsActuation > maximalRedundantJointsActuation_));
      
      const arma::Cube<double>& model = getModel(endEffectorPose, redundantJointsActuation);

      const arma::Mat<double>& baseJointPositions = model.slice(1);

      const arma::Mat<double>& endEffectorJointPositions = model.slice(1);
      arma::Mat<double> endEffectorJointPositionsRotated = endEffectorJointPositions;
      endEffectorJointPositionsRotated.each_col() -= endEffectorPose.subvec(0, 1);

      const arma::Mat<double>& baseToEndEffectorJointPositions = endEffectorJointPositions - baseJointPositions;
      const arma::Row<double>& baseToEndEffectorJointsActuation = arma::sqrt(arma::sum(arma::square(baseToEndEffectorJointPositions)));

      if (arma::any(baseToEndEffectorJointsActuation < minimalActiveJointsActuation_) || arma::any(baseToEndEffectorJointsActuation > maximalActiveJointsActuation_)) {
        return 0.0;
      }

      arma::Mat<double> forwardKinematic;
      forwardKinematic.head_rows(2) = baseToEndEffectorJointPositions;
      forwardKinematic.row(2) = -forwardKinematic.row(0) % endEffectorJointPositionsRotated.row(1) + forwardKinematic.row(1) % endEffectorJointPositionsRotated.row(0);

      arma::Mat<double> inverseKinematic(3, 3 + redundantJointIndicies_.n_elem, arma::fill::zeros);
      inverseKinematic.diag() = -baseToEndEffectorJointsActuation;
      for (arma::uword n = 0; n < redundantJointIndicies_.n_elem; ++n) {
        const arma::uword& redundantJointIndex = redundantJointIndicies_(n);
        inverseKinematic(n, 3 + n) = forwardKinematic(redundantJointIndex, 0) * redundantJointAngleCosines_(n) + forwardKinematic(redundantJointIndex, 1) * redundantJointAngleSines_(n);
      }

      return -1.0 / arma::cond(arma::solve(forwardKinematic.t(), inverseKinematic));
    }
    
    std::string ParallelKinematicMachine3PRPR::toString() const {
      return "robotics_parallel_kinematic_machine_3prpr";
    }
  }
}
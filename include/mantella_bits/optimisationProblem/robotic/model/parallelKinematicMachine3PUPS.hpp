namespace mant {
  namespace robotic {
    class ParallelKinematicMachine3PUPS {
      public:
        inline explicit ParallelKinematicMachine3PUPS() noexcept;

        inline arma::Row<double>::fixed<3> getMinimalActiveJointActuations() const noexcept;

        inline void setMinimalActiveJointActuations(
            const arma::Row<double>::fixed<3>& minimalActiveJointActuations) noexcept;

        inline arma::Row<double>::fixed<3> getMaximalActiveJointActuations() const noexcept;

        inline void setMaximalActiveJointActuations(
            const arma::Row<double>::fixed<3>& maximalActiveJointActuations) noexcept;

        inline arma::Mat<double>::fixed<3, 3> getEndEffectorJointPositions() const noexcept;

        inline void setEndEffectorJointPositions(
            const arma::Mat<double>::fixed<3, 3>& endEffectorJointPositions) noexcept;

        inline arma::Mat<double>::fixed<3, 3> getRedundantJointStartPositions() const noexcept;

        arma::Mat<double>::fixed<3, 3> endEffectorJointsRelative_;
        inline void setRedundantJointStartPositions(
            const arma::Mat<double>::fixed<3, 3>& redundantJointStartPositions) noexcept;

        inline arma::Mat<double>::fixed<3, 3> getRedundantJointEndPositions() const noexcept;

        arma::Mat<double>::fixed<3, 3> redundantJointStarts_;
        arma::Mat<double>::fixed<3, 3> redundantJointEnds_;
        inline void setRedundantJointEndPositions(
            const arma::Mat<double>::fixed<3, 3>& redundantJointEndPositions) noexcept;

        inline std::vector<arma::Mat<double>::fixed<3, 3>> getModel(
            const arma::Col<double>::fixed<6>& endEffectorPose,
            const arma::Row<double>& redundantJointActuations) const;

        inline arma::Row<double>::fixed<3> getActuation(
            const arma::Col<double>::fixed<6>& endEffectorPose,
            const arma::Row<double>& redundantJointActuations) const;

        inline arma::Col<double>::fixed<6> getEndEffectorPose(
            const arma::Row<double>::fixed<3>& actuations,
            const arma::Row<double>& redundantJointActuations) const;

        inline double getEndEffectorPoseAccuracy(
            const arma::Col<double>::fixed<6>& endEffectorPose,
            const arma::Row<double>& redundantJointActuations) const;

      protected:
        arma::Row<double>::fixed<3> minimalActiveJointActuations_;
        arma::Row<double>::fixed<3> maximalActiveJointActuations_;

        arma::Mat<double>::fixed<3, 3> redundantJointsStartToEnd_;
        arma::Col<unsigned int> redundantJointIndicies_;
        arma::Mat<double> redundantJointAngles_;
    };

    //
    // Implementation
    //

    inline ParallelKinematicMachine3PUPS::ParallelKinematicMachine3PUPS() noexcept {
      setMinimalActiveJointActuations({0.39, 0.39, 0.39});
      setMaximalActiveJointActuations({0.95, 0.95, 0.95});

      setEndEffectorJointPositions({
        -0.025561381023353, 0.086293776138137, 0.12,
        0.087513292835791, -0.021010082747031, 0.12,
        -0.061951911812438, -0.065283693391106, 0.12});

      setRedundantJointStartPositions({
        -0.463708870031622, 0.417029254828353, -0.346410161513775,
        0.593012363818459, 0.193069033993384, -0.346410161513775,
        -0.129303493786837, -0.610098288821738, -0.346410161513775});

      setRedundantJointEndPositions({
        -0.247202519085512, 0.292029254828353, 0.086602540378444,
        0.376506012872349, 0.068069033993384, 0.086602540378444,
        -0.129303493786837, -0.360098288821738, 0.086602540378444});

      redundantJointStartToEndPositions_ = redundantJointEndPositions_ - redundantJointStartPositions_;
      redundantJointIndicies_ = arma::find(arma::any(redundantJointStartToEndPositions_));

      redundantJointAngles_.set_size(3, redundantJointIndicies_.n_elem);

      for (std::size_t n = 0; n < redundantJointIndicies_.n_elem; ++n) {
        const double& redundantJointXAngle = std::atan2(redundantJointsStartToEnd_.at(1, n), redundantJointsStartToEnd_.at(0, n));
        const double& redundantJointYAngle = std::atan2(redundantJointsStartToEnd_.at(2, n), redundantJointsStartToEnd_.at(1, n));
        redundantJointAngles_.col(n) = arma::Col<double>::fixed<3>({std::cos(redundantJointXAngle) * std::cos(redundantJointYAngle), std::sin(redundantJointXAngle) * std::cos(redundantJointYAngle), std::sin(redundantJointYAngle)});
      }
    }

    inline arma::Row<double>::fixed<3> ParallelKinematicMachine3PUPS::getMinimalActiveJointActuations() const noexcept {
      return minimalActiveJointActuations_;
    }

    inline void ParallelKinematicMachine3PUPS::setMinimalActiveJointActuations(
        const arma::Row<double>::fixed<3>& minimalActiveJointActuations) noexcept {
      minimalActiveJointActuations_ = minimalActiveJointActuations;
    }

    inline arma::Row<double>::fixed<3> ParallelKinematicMachine3PUPS::getMaximalActiveJointActuations() const noexcept {
      return maximalActiveJointActuations_;
    }

    inline void ParallelKinematicMachine3PUPS::setMaximalActiveJointActuations(
        const arma::Row<double>::fixed<3>& maximalActiveJointActuations) noexcept {
      maximalActiveJointActuations_ = maximalActiveJointActuations;
    }

    inline arma::Mat<double>::fixed<3, 3> ParallelKinematicMachine3PUPS::getEndEffectorJointPositions() const noexcept {
      return endEffectorJointPositions_;
    }

    inline void ParallelKinematicMachine3PUPS::setEndEffectorJointPositions(
        const arma::Mat<double>::fixed<3, 3>& endEffectorJointPositions) noexcept {
      endEffectorJointPositions_ = endEffectorJointPositions;
    }

    inline arma::Mat<double>::fixed<3, 3> ParallelKinematicMachine3PUPS::getRedundantJointStartPositions() const noexcept {
      return redundantJointStartPositions_;
    }

    inline void ParallelKinematicMachine3PUPS::setRedundantJointStartPositions(
        const arma::Mat<double>::fixed<3, 3>& redundantJointStartPositions) noexcept {
      redundantJointStartPositions_ = redundantJointStartPositions;
    }

    inline arma::Mat<double>::fixed<3, 3> ParallelKinematicMachine3PUPS::getRedundantJointEndPositions() const noexcept {
      return redundantJointEndPositions_;
    }

    inline void ParallelKinematicMachine3PUPS::setRedundantJointEndPositions(
        const arma::Mat<double>::fixed<3, 3>& redundantJointEndPositions) noexcept {
      redundantJointEndPositions_ = redundantJointEndPositions;
    }

    inline std::vector<arma::Mat<double>::fixed<3, 3>> ParallelKinematicMachine3PUPS::getModel(
        const arma::Col<double>::fixed<6>& endEffectorPose,
        const arma::Row<double>& redundantJointActuations) const {
      if (arma::any(arma::vectorise(redundantJointActuations < 0)) || arma::any(arma::vectorise(redundantJointActuations > 1))) {
        throw std::logic_error("All values for the actuation of redundantion joints must be between [0, 1].");
      }
      // TODO Check number of redundantJointActuations vs redudantant elements

      const arma::Col<double>::fixed<3>& endEffectorPosition = endEffectorPose.subvec(0, 2);
      const double& endEffectorRollAngle = endEffectorPose.at(3);
      const double& endEffectorPitchAngle = endEffectorPose.at(4);
      const double& endEffectorYawAngle = endEffectorPose.at(5);

      arma::Mat<double>::fixed<3, 3> baseJoints = redundantJointStarts_;
      for (std::size_t n = 0; n < redundantJointIndicies_.n_elem; n++) {
        const unsigned int& redundantJointIndex = redundantJointIndicies_.at(n);
        baseJoints.col(redundantJointIndex) += redundantJointActuations.at(redundantJointIndex) * redundantJointsStartToEnd_.col(redundantJointIndex);
      }

      arma::Mat<double>::fixed<3, 3> endEffectorJoints = get3DRotationMatrix(endEffectorRollAngle, endEffectorPitchAngle, endEffectorYawAngle) * endEffectorJointsRelative_;
      endEffectorJoints.each_col() += endEffectorPosition;

      std::vector<arma::Mat<double>::fixed<3, 3>> model;

      model.push_back(baseJoints);
      model.push_back(endEffectorJoints);

      return model;
    }

    inline arma::Col<double>::fixed<6> ParallelKinematicMachine3PUPS::getEndEffectorPose(
        const arma::Row<double>::fixed<3>& actuations,
        const arma::Row<double>& redundantJointActuations) const {
      // TODO Direct kinematic (estimate position, using a simple HillCLimber algorithm)
      return {0, 0, 0, 0, 0, 0};
    }

    inline arma::Row<double>::fixed<3> ParallelKinematicMachine3PUPS::getActuation(
        const arma::Col<double>::fixed<6>& endEffectorPose,
        const arma::Row<double>& redundantJointActuations) const {
      const std::vector<arma::Mat<double>::fixed<3, 3>>& model = getModel(endEffectorPose, redundantJointActuations);

      const arma::Mat<double>::fixed<3, 3>& baseJoints = model.at(0);
      const arma::Mat<double>::fixed<3, 3>& endEffectorJoints = model.at(1);

      return arma::sqrt(arma::sum(arma::square(endEffectorJoints - baseJoints)));
    }

    inline double ParallelKinematicMachine3PUPS::getEndEffectorPoseAccuracy(
        const arma::Col<double>::fixed<6>& endEffectorPose,
        const arma::Row<double>& redundantJointActuations) const {
      const std::vector<arma::Mat<double>::fixed<3, 3>>& model = getModel(endEffectorPose, redundantJointActuations);

      const arma::Mat<double>::fixed<3, 3>& baseJoints = model.at(0);

      const arma::Mat<double>::fixed<3, 3>& endEffectorJoints = model.at(1);
      arma::Mat<double>::fixed<3, 3> endEffectorJointsRotated = endEffectorJoints;
      endEffectorJointsRotated.each_col() -= endEffectorPose.subvec(0, 1);

      const arma::Mat<double>::fixed<3, 3>& baseToEndEffectorJointPositions = endEffectorJoints - baseJoints;
      const arma::Row<double>::fixed<3>& baseToEndEffectorJointActuations = arma::sqrt(arma::sum(arma::square(baseToEndEffectorJointPositions)));

      if (arma::any(baseToEndEffectorJointActuations < minimalActiveJointActuations_) || arma::any(baseToEndEffectorJointActuations > maximalActiveJointActuations_)) {
        return 0.0;
      }

      arma::Mat<double>::fixed<6, 3> forwardKinematic;
      forwardKinematic.rows(0, 2) = baseToEndEffectorJointPositions;
      for (std::size_t n = 0; n < baseToEndEffectorJointPositions.n_cols; ++n) {
        forwardKinematic.submat(3, n, 5, n) = arma::cross(endEffectorJointsRotated.col(n), baseToEndEffectorJointPositions.col(n));
      }

      arma::Mat<double> inverseKinematic(6, 3 + redundantJointIndicies_.n_elem, arma::fill::zeros);
      inverseKinematic.diag() = -arma::sqrt(arma::sum(arma::square(baseToEndEffectorJointPositions)));
      for (std::size_t n = 0; n < redundantJointIndicies_.n_elem; ++n) {
        const unsigned int& redundantJointIndex = redundantJointIndicies_.at(n);
        inverseKinematic.at(n, 3 + n) = arma::dot(baseToEndEffectorJointPositions.col(redundantJointIndex), redundantJointAngles_.col(redundantJointIndex));
      }

      return -1.0 / arma::cond(arma::solve(forwardKinematic.t(), inverseKinematic));
    }
  }
}
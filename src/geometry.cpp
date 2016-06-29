#include "mantella_bits/geometry.hpp"
#include "mantella_bits/config.hpp"

// C++ standard library
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace mant {
  arma::mat::fixed<2, 2> rotationMatrix2d(
      const double angle) {
    // clang-format off
    return arma::mat::fixed<2, 2>({
      std::cos(angle), -std::sin(angle),
      std::sin(angle),  std::cos(angle)});
    // clang-format on
  }

  arma::mat::fixed<3, 3> rotationMatrix3d(
      const double rollAngle,
      const double pitchAngle,
      const double yawAngle) {
    // clang-format off
    return arma::mat::fixed<3, 3>({
      // First column
      std::cos(yawAngle) * std::cos(pitchAngle),
      std::cos(yawAngle) * std::sin(pitchAngle) * std::sin(rollAngle) - std::sin(yawAngle) * std::cos(rollAngle),
      std::cos(yawAngle) * std::sin(pitchAngle) * std::cos(rollAngle) + std::sin(yawAngle) * std::sin(rollAngle),
      // Second column
      std::sin(yawAngle) * std::cos(pitchAngle),
      std::sin(yawAngle) * std::sin(pitchAngle) * std::sin(rollAngle) + std::cos(yawAngle) * std::cos(rollAngle),
      std::sin(yawAngle) * std::sin(pitchAngle) * std::cos(rollAngle) - std::cos(yawAngle) * std::sin(rollAngle),
      // Last column
      -std::sin(pitchAngle),
      std::cos(pitchAngle) * std::sin(rollAngle),
      std::cos(pitchAngle) * std::cos(rollAngle)});
    // clang-format on
  }

  std::vector<arma::vec::fixed<2>> circleCircleIntersections(
      const arma::vec::fixed<2>& firstCentre,
      const double firstRadius,
      const arma::vec::fixed<2>& secondCentre,
      const double secondRadius) {
    if (firstRadius < 0) {
      throw std::domain_error("circleCircleIntersections: The first radius must be positive (including 0).");
    } else if (secondRadius < 0) {
      throw std::domain_error("circleCircleIntersections: The second radius must be positive (including 0).");
    }

    /* The circle circle intersection points are calculated as follows:
     * 0. We assume that both centers are on the x-axis and `firstCentre` is at (0, 0). 
     *    **Note:** This assumptions are lifted later on.
     * 1. Calculate the distance between both centers.
     * 2. Calculate the x-part of the intersection point.
     *
     *   y-axis
     *      ^
     *      |        
     *      |
     *      |----- distance ----|
     *      |                   |   
     *      +-----------------------------> x-axis
     *      ^         ^         ^
     *      |         |         |
     * `firstCenter`  |   `secondCentrer`
     *                |
     *        intersection point (only the x-part)
     *
     * After the first to steps we should get a picture like above.
     * At this point, we either got 2 intersection points, when the distance between the x-part of the intersection point and `firstCenter` is less than `firstRadius`, 1 intersection point if both are equal and no intersection if its greater.
     *
     * A: We got 1 intersection
     *   3. Generate a unit vector pointing from `firstCentre` to `secondCentre`.
     *   4. Set the y-part to 0.
     *   5. Scale, rotate and translate the (x, y)-coordinate to be within the actual coordinate system (removing the assumption from 1.)
     *   
     * B: We got 2 intersections
     *   3. Generate a unit vector pointing from `firstCentre` to `secondCentre`.
     *   4. Calculate the y-part.
     *   5. Scale, rotate and translate the (x, y)-coordinate to be within the actual coordinate system (removing the assumption from 1.)
     */

    const double distance = arma::norm(secondCentre - firstCentre);
    if (distance < ::mant::machinePrecision && std::abs(firstRadius - secondRadius) < ::mant::machinePrecision) {
      // Both circles are identical ...
      if (firstRadius > ::mant::machinePrecision) {
        // ... but have a non-zero radius.
        throw std::invalid_argument("circleCircleIntersections: Both centers and radii (> 0) are identical, resulting in infinite intersections.");
      }
      // ... and dots.
      return {firstCentre};
    }

    if (distance > firstRadius + secondRadius || distance < std::abs(firstRadius - secondRadius)) {
      // Both circles are either to far away or to close.
      return {};
    }

    const double x = (std::pow(firstRadius, 2.0) - std::pow(secondRadius, 2.0) + std::pow(distance, 2.0)) / (2.0 * distance);
    const arma::vec::fixed<2>& unitVector = (secondCentre - firstCentre) / distance;

    if (std::abs(firstRadius - std::abs(x)) < ::mant::machinePrecision) {
      // One intersection
      return std::vector<arma::vec::fixed<2>>({{firstCentre + unitVector * x}});
    } else {
      // Two intersections
      const double y = std::sqrt(std::pow(firstRadius, 2.0) - std::pow(x, 2.0));
      return std::vector<arma::vec::fixed<2>>({{firstCentre(0) + unitVector(0) * x - unitVector(1) * y, firstCentre(1) + unitVector(1) * x + unitVector(0) * y}, {firstCentre(0) + unitVector(0) * x + unitVector(1) * y, firstCentre(1) + unitVector(1) * x - unitVector(0) * y}});
    }
  }

  std::vector<arma::vec::fixed<3>> circleSphereIntersections(
      const arma::vec::fixed<3>& circleCentre,
      const double circleRadius,
      const arma::vec::fixed<3>& circleNormal,
      const arma::vec::fixed<3>& sphereCentre,
      const double sphereRadius) {
    if (circleRadius < 0) {
      throw std::domain_error("circleSphereIntersections: The circle radius must be positive (including 0).");
    } else if (sphereRadius < 0) {
      throw std::domain_error("circleSphereIntersections: The sphere radius must be positive (including 0).");
    }

    /* The circle sphere intersection points are calculated as follows:
     * 0. We assume that both centers are on the x-axis, ´circleNormal` is perpendicular to the x- and y-axis and `circleCentre` is at (0, 0, 0). 
     *    **Note:** This assumptions are lifted later on.
     * 1. Calculate the shortest distance between the sphere's centre and circle's plane.
     * 2. Calculate the centre of the sphere's inner circle segment, placed an the same plane as the given circle.
     * 3. Calculate the radius of the sphere's inner circle segment.
     *    Given the provided circle and the sphere's inner circle, the problem is now reduced to a circle circle intersection.
     * 4. Calculate the x-part of the intersection point.
     *    At this point, we either got 2 intersection points, when the distance between the x-part of the intersection point and `circleCentre` is less than `circleRadius`, 1 intersection point if both are equal and no intersection if its greater.
     *
     * A: We got 1 intersection
     *   5. Generate a unit vector pointing from `circleCentre` to `innerCentre`.
     *   6. Set the y-part to 0.
     *   7. Scale and translate the unit vector to be within the actual coordinate system (removing the assumption from 1.).
     *   
     * B: We got 2 intersections
     *   5. Generate a unit vector pointing from `circleCentre` to `innerCentre`.
     *   6. Calculate the y-part.
     *   7. Generate a second unit vector, perpendicular to `xUnitVector` and `circleNormal`.
     *   8. Scale and translate both unit vectors to be within the actual coordinate system (removing the assumption from 1.).
     */

    const double innerDistance = arma::dot(circleNormal, circleCentre - sphereCentre);
    if (std::abs(innerDistance) - sphereRadius > ::mant::machinePrecision) {
      // The circle's plane does not intersect with the sphere.
      return {};
    };

    const arma::vec::fixed<3>& innerCentre = sphereCentre + innerDistance * circleNormal;
    // Due to rounding errors, the inner result might be negative instead of being 0.
    const double innerRadius = std::sqrt(std::max(0.0, std::pow(sphereRadius, 2.0) - std::pow(innerDistance, 2.0)));

    const double distance = arma::norm(innerCentre - circleCentre);
    if (distance < ::mant::machinePrecision && std::abs(innerRadius - circleRadius) < ::mant::machinePrecision) {
      // Both circles are identical ...
      if (circleRadius > ::mant::machinePrecision) {
        // ... but have a non-zero radius.
        throw std::invalid_argument("circleSphereIntersections: Both centers and radii (> 0) are identical, resulting in infinite intersections.");
      }

      // ... and dots.
      return {circleCentre};
    }

    if (distance - circleRadius - innerRadius > ::mant::machinePrecision || std::abs(circleRadius - innerRadius) - distance >= ::mant::machinePrecision) {
      // Both circles are either to far away or to close.
      return {};
    };

    const double x = (std::pow(circleRadius, 2.0) - std::pow(innerRadius, 2.0) + std::pow(distance, 2.0)) / (2.0 * distance);
    const arma::vec::fixed<3>& xUnitVector = (innerCentre - circleCentre) / distance;

    if (std::abs(circleRadius - std::abs(x)) < ::mant::machinePrecision) {
      // One intersection
      return std::vector<arma::vec::fixed<3>>({circleCentre + x * xUnitVector});
    } else {
      // Two intersections
      const double y = std::sqrt(std::pow(circleRadius, 2.0) - std::pow(x, 2.0));
      const arma::vec::fixed<3>& yUnitVector = arma::normalise(arma::cross(xUnitVector, circleNormal));
      return std::vector<arma::vec::fixed<3>>({circleCentre + x * xUnitVector + y * yUnitVector, circleCentre + x * xUnitVector - y * yUnitVector});
    }
  }
}

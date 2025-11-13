#ifndef FIX_EIGEN5_ISFINITE_H
#define FIX_EIGEN5_ISFINITE_H

#include <cmath>
#include <Eigen/Core>

namespace Eigen {
namespace internal {

// Fix for Eigen 5.0.1 missing isfinite_impl for real numbers
template<typename Scalar>
EIGEN_DEVICE_FUNC inline bool isfinite_impl(const Scalar& x) {
  return std::isfinite(x);
}

// Explicit specializations for float and double
template<>
EIGEN_DEVICE_FUNC inline bool isfinite_impl<float>(const float& x) {
  return std::isfinite(x);
}

template<>
EIGEN_DEVICE_FUNC inline bool isfinite_impl<double>(const double& x) {
  return std::isfinite(x);
}

} // namespace internal
} // namespace Eigen

#endif // FIX_EIGEN5_ISFINITE_H
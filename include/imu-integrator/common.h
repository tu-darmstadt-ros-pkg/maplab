#ifndef IMU_INTEGRATOR_COMMON_H_
#define IMU_INTEGRATOR_COMMON_H_

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <glog/logging.h>

namespace imu_integrator {

static const int kStateSize = 16;
static const int kErrorStateSize = 15;

static const int kStateOrientationOffset = 0;
static const int kStateGyroBiasOffset = 4;
static const int kStateVelocityOffset = 7;
static const int kStateAccelBiasOffset = 10;
static const int kStatePositionOffset = 13;

static const int kErrorStateOrientationOffset = 0;
static const int kErrorStateGyroBiasOffset = 3;
static const int kErrorStateVelocityOffset = 6;
static const int kErrorStateAccelBiasOffset = 9;
static const int kErrorStatePositionOffset = 12;

static const int kStateOrientationBlockSize = 4;
static const int kErrorOrientationBlockSize = 3;
static const int kGyroBiasBlockSize = 3;
static const int kVelocityBlockSize = 3;
static const int kAccelBiasBlockSize = 3;
static const int kPositionBlockSize = 3;
static const int kStatePoseBlockSize = 7;

// New Eigen based error terms combine IMU biases
static const int kImuBiasBlockSize = 6;

static const int kImuReadingSize = 6;
static const int kAccelReadingOffset = 0;
static const int kGyroReadingOffset = 3;

static const double kNanoSecondsToSeconds = 1e-9;

}  // namespace imu_integrator



namespace common {
    template<typename Derived>
    inline Eigen::Matrix<typename Derived::Scalar, 3, 3> skew(
            const Eigen::MatrixBase<Derived> &vector) {
        EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, 3);

        Eigen::Matrix<typename Derived::Scalar, 3, 3> output;
        typename Derived::Scalar zero = static_cast<typename Derived::Scalar>(0);
        output << zero, -vector[2], vector[1], vector[2], zero, -vector[0],
                -vector[1], vector[0], zero;
        return output;
    }

    template<typename Derived, typename OtherDerived>
    void skew(
            const Eigen::MatrixBase<Derived> &vector,
            Eigen::MatrixBase<OtherDerived> const &matrix_const) {
        EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
        EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(OtherDerived, 3, 3);

        typedef typename OtherDerived::Scalar Scalar;

        Eigen::MatrixBase<OtherDerived> &matrix =
                const_cast<Eigen::MatrixBase<OtherDerived> &>(matrix_const);

        Scalar zero = static_cast<Scalar>(0.0);

        matrix.derived() << zero, -vector[2], vector[1], vector[2], zero, -vector[0],
                -vector[1], vector[0], zero;
    }


    template <typename Derived, typename ScalarType>
    void toRotationMatrixJPL(
            const Eigen::MatrixBase<Derived>& q,
            Eigen::Matrix<ScalarType, 3, 3>* rot_matrix) {
        CHECK_NOTNULL(rot_matrix);
        EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, 4);

        ScalarType one = static_cast<ScalarType>(1.0);
        ScalarType two = static_cast<ScalarType>(2.0);

        (*rot_matrix)(0, 0) = one - two * (q(1) * q(1) + q(2) * q(2));
        (*rot_matrix)(0, 1) = two * (q(0) * q(1) + q(2) * q(3));
        (*rot_matrix)(0, 2) = two * (q(0) * q(2) - q(1) * q(3));

        (*rot_matrix)(1, 0) = two * (q(0) * q(1) - q(2) * q(3));
        (*rot_matrix)(1, 1) = one - two * (q(0) * q(0) + q(2) * q(2));
        (*rot_matrix)(1, 2) = two * (q(1) * q(2) + q(0) * q(3));

        (*rot_matrix)(2, 0) = two * (q(0) * q(2) + q(1) * q(3));
        (*rot_matrix)(2, 1) = two * (q(1) * q(2) - q(0) * q(3));
        (*rot_matrix)(2, 2) = one - two * (q(0) * q(0) + q(1) * q(1));
    }

// Conversion from quaternion to euler angles.
    template <typename ScalarType>
    Eigen::Matrix<ScalarType, 3, 1> getRollPitchYawFromQuaternion(
            const Eigen::Quaternion<ScalarType>& q) {
        Eigen::Matrix<ScalarType, 3, 1> rpy_rad;
        rpy_rad(0) = std::atan2(
                2.0 * (q.w() * q.x() + q.y() * q.z()),
                1.0 - 2.0 * (q.x() * q.x() + q.y() * q.y()));
        rpy_rad(1) = std::asin(2.0 * (q.w() * q.y() - q.z() * q.x()));
        rpy_rad(2) = std::atan2(
                2.0 * (q.w() * q.z() + q.x() * q.y()),
                1.0 - 2.0 * (q.y() * q.y() + q.z() * q.z()));
        return rpy_rad;
    }

    namespace eigen_quaternion_helpers {
        const int kLocalSize = 3;
        const int kGlobalSize = 4;

        template<typename Scalar>
        inline Eigen::Matrix<Scalar, 3, 3> Gamma(
                const Eigen::Matrix<Scalar, 3, 1> &phi) {
            const Scalar phi_squared_norm = phi.squaredNorm();

            if (phi_squared_norm < 1e-6) {
                Eigen::Matrix<Scalar, 3, 3> gamma;
                gamma.setIdentity();
                gamma += 0.5 * common::skew(phi);
                return gamma;
            }
            const Scalar phi_norm = sqrt(phi_squared_norm);
            const Eigen::Matrix<Scalar, 3, 3> phi_skew(common::skew(phi));

            Eigen::Matrix<Scalar, 3, 3> gamma;
            gamma.setIdentity();
            gamma += ((1.0 - cos(phi_norm)) / phi_squared_norm) * phi_skew;
            const Scalar phi_cubed = (phi_norm * phi_squared_norm);
            gamma += ((phi_norm - sin(phi_norm)) / phi_cubed) * phi_skew * phi_skew;
            return gamma;
        }

        template<typename Scalar>
        inline Eigen::Quaternion<Scalar> ExpMap(
                const Eigen::Matrix<Scalar, 3, 1> &theta) {
            const Scalar theta_squared_norm = theta.squaredNorm();

            if (theta_squared_norm < 1e-6) {
                Eigen::Quaternion<Scalar> q(
                        1, theta(0) * 0.5, theta(1) * 0.5, theta(2) * 0.5);
                q.normalize();
                return q;
            }

            const Scalar theta_norm = sqrt(theta_squared_norm);
            const Eigen::Matrix<Scalar, 3, 1> q_imag =
                    sin(theta_norm * 0.5) * theta / theta_norm;
            Eigen::Quaternion<Scalar> q(
                    cos(theta_norm * 0.5), q_imag(0), q_imag(1), q_imag(2));
            return q;
        }

        inline Eigen::Quaterniond ExpMap(const Eigen::Vector3d &theta) {
            return ExpMap<double>(theta);
        }

    }
}


#endif  // IMU_INTEGRATOR_COMMON_H_

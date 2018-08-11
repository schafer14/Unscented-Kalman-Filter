#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
// if this is false, laser measurements will be ignored (except during init)
use_laser_ = true;

// if this is false, radar measurements will be ignored (except during init)
use_radar_ = true;

n_x_ = 5;
n_aug_ = 7;
n_sig_ = 2 * n_aug_ + 1;

// initial state vector
x_ = VectorXd(n_x_);

// initial covariance matrix
P_ = MatrixXd(n_x_, n_x_);
P_ << 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 0, 1;

// Process noise standard deviation longitudinal acceleration in m/s^2
std_a_ = 1.5;

// Process noise standard deviation yaw acceleration in rad/s^2
std_yawdd_ = 0.5;

//DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
// Laser measurement noise standard deviation position1 in m
std_laspx_ = 0.15;

// Laser measurement noise standard deviation position2 in m
std_laspy_ = 0.15;

// Radar measurement noise standard deviation radius in m
std_radr_ = 0.3;

// Radar measurement noise standard deviation angle in rad
std_radphi_ = 0.03;

// Radar measurement noise standard deviation radius change in m/s
std_radrd_ = 0.3;
//DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

/**
TODO:

Complete the initialization. See ukf.h for other member properties.

Hint: one or more values initialized above might be wildly off...
*/

is_initialized_ = false;
lambda_ = 3 - n_x_;
weights_ = VectorXd(n_sig_);
weights_.fill(0.5 / (n_aug_ + lambda_));
weights_(0) = lambda_ / (lambda_ + n_aug_);

R_radar_ = MatrixXd(3, 3);
R_radar_<< std_radr_ * std_radr_, 0, 0,
    0, std_radphi_ * std_radphi_, 0,
    0, 0, std_radrd_ * std_radrd_;

R_lidar_ = MatrixXd(2, 2);
R_lidar_ << std_laspx_ * std_laspx_, 0,
            0, std_laspy_ * std_laspy_;
}

UKF::~UKF() {}

/**
* @param {MeasurementPackage} meas_package The latest measurement data of
* either radar or laser.
*/
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
/**
TODO:

Complete this function! Make sure you switch between lidar and radar
measurements.
*/
  cout << "Starting Loop" << endl;
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double px = meas_package.raw_measurements_[0] * cos(meas_package.raw_measurements_[1]);
      double py = meas_package.raw_measurements_[0] * sin(meas_package.raw_measurements_[1]);
      x_ << px, py, 0, 0, 0;
    }
    else {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
  }

  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  } 
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  }
  cout << "Ending loop" << endl;
}

/**
* Predicts sigma points, the state, and the state covariance matrix.
* @param {double} delta_t the change in time (in seconds) between the last
* measurement and this one.
*/
void UKF::Prediction(double delta_t) {
/**
TODO:

Complete this function! Estimate the object's location. Modify the state
vector, x_. Predict sigma points, the state, and the state covariance matrix.
*/
  cout << "starting Predict" << endl;
  // GENERATE SIGMA POINTS
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  x_aug.fill(0.0);
  x_aug.head(5) = x_;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  MatrixXd L = P_aug.llt().matrixL();

  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  Xsig_pred_ = MatrixXd(n_x_, n_sig_);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double dyaw = Xsig_aug(4, i);
    double n_a = Xsig_aug(5, i);
    double n_yaw = Xsig_aug(6, i);


    double ppx, ppy;
    if (fabs(dyaw) > 0.001) {
      ppx = px + v / dyaw * (sin(yaw + dyaw * delta_t) - sin(yaw));
      ppy = py + v / dyaw * (cos(yaw) - cos(yaw + dyaw * delta_t));
    } else {
      ppx = px + v * cos(yaw) * delta_t;
      ppy = py + v * sin(yaw) * delta_t;
    }

    ppx = ppx + .5 * delta_t * delta_t * cos(yaw) * n_a;
    ppy = ppy + .5 * delta_t * delta_t * sin(yaw) * n_a;
    double pv = v + delta_t * n_a;
    double pyaw = yaw + dyaw * delta_t + .5 * delta_t * delta_t * n_yaw;
    double pdyaw = dyaw + delta_t * n_yaw;

    Xsig_pred_(0, i) = ppx;
    Xsig_pred_(1, i) = ppy;
    Xsig_pred_(2, i) = pv;
    Xsig_pred_(3, i) = pyaw;
    Xsig_pred_(4, i) = pdyaw;
  }

  x_ = Xsig_pred_ * weights_;

  P_.fill(0.0);
  for (int i = 0; i < n_sig_; i ++) {
    VectorXd xdiff = Xsig_pred_.col(i) - x_;

    if (xdiff(3) > M_PI) {
      xdiff(3) -= 2. * M_PI;
    }
    else if (xdiff(3) < -M_PI) {
      xdiff(3) += 2. * M_PI;
    }
        
    P_ = P_ + weights_(i) * xdiff * xdiff.transpose();
  }
}

/**
* Updates the state and the state covariance matrix using a laser measurement.
* @param {MeasurementPackage} meas_package
*/
void UKF::UpdateLidar(MeasurementPackage meas_package) {
/**
TODO:

Complete this function! Use lidar data to update the belief about the object's
position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  cout << "Calculating Lidar" << endl;
  int n_z_ = 2;
  VectorXd z = meas_package.raw_measurements_;
  MatrixXd Zsig = MatrixXd(n_z_, n_sig_);
  VectorXd z_pred = VectorXd(n_z_);
  MatrixXd S = MatrixXd(n_z_, n_z_);
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  for (int i = 0; i < n_sig_; i++) {
    Zsig.col(i) = Xsig_pred_.col(i).head(2);    

    z_pred += weights_(i) * Zsig.col(i);
  }

  S.fill(0.0);
  for (int i = 0; i < n_sig_; i ++) {
    VectorXd zdiff = Zsig.col(i) = z_pred;

    while (zdiff(1) > M_PI) zdiff(1) -= 2. * M_PI;
    while (zdiff(1) < -M_PI) zdiff(1) += 2. * M_PI;

    S += weights_(i) * zdiff * zdiff.transpose();
  }

  S += R_lidar_;

  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; i ++) {
    VectorXd xdiff = Xsig_pred_.col(i) - x_;
    VectorXd zdiff = Zsig.col(i) - z_pred;

    while (zdiff(1)> M_PI) zdiff(1)-=2.*M_PI;
    while (zdiff(1)<-M_PI) zdiff(1)+=2.*M_PI;

    while (xdiff(3)> M_PI) xdiff(3)-=2.*M_PI;
    while (xdiff(3)<-M_PI) xdiff(3)+=2.*M_PI;

    Tc += weights_(i) * xdiff * zdiff.transpose();
  }

  MatrixXd K = Tc * S.inverse();

  VectorXd zdiff = z - z_pred;

  while (zdiff(1)> M_PI) zdiff(1)-=2.*M_PI;
  while (zdiff(1)<-M_PI) zdiff(1)+=2.*M_PI;

  x_ = x_ + K * zdiff;
  P_ = P_ - K * S * K.transpose();

  cout << "Ending Laser" << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  cout << "Starting Radar" << endl;
  int n_z_ = 3;
  VectorXd z = meas_package.raw_measurements_;
  MatrixXd Zsig = MatrixXd(n_z_, n_sig_);
  VectorXd z_pred = VectorXd(n_z_);
  MatrixXd S = MatrixXd(n_z_, n_z_);
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  Zsig.fill(0.0);
  z_pred.fill(0.0);
  for (int i = 0; i < n_sig_; i ++) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    double p = sqrt(px  * px + py * py);
    double row = atan2(py, px);
    double pdot = (px * v1 + py * v2) / p;

    Zsig(0, i) = p;
    Zsig(1, i) = row;
    Zsig(2, i) = pdot;

    z_pred += weights_(i) * Zsig.col(i);
  }

  S.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    VectorXd zdiff = Zsig.col(i) - z_pred;
    
    while (zdiff(1) > M_PI) zdiff(1) -= 2. * M_PI;
    while (zdiff(1) < -M_PI) zdiff(1) += 2. * M_PI;

    S += weights_(i) * zdiff * zdiff.transpose();
  }


  S += R_radar_;

  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; i ++) {
    VectorXd xdiff = Xsig_pred_.col(i) - x_;
    VectorXd zdiff = Zsig.col(i) - z_pred;

    while (zdiff(1)> M_PI) zdiff(1)-=2.*M_PI;
    while (zdiff(1)<-M_PI) zdiff(1)+=2.*M_PI;

    while (xdiff(3)> M_PI) xdiff(3)-=2.*M_PI;
    while (xdiff(3)<-M_PI) xdiff(3)+=2.*M_PI;

    Tc += weights_(i) * xdiff * zdiff.transpose();
  }

  MatrixXd K = Tc * S.inverse();

  VectorXd zdiff = z - z_pred;

  while (zdiff(1)> M_PI) zdiff(1)-=2.*M_PI;
  while (zdiff(1)<-M_PI) zdiff(1)+=2.*M_PI;

  x_ = x_ + K * zdiff;
  P_ = P_ - K * S * K.transpose();

  cout << "Ending Radar" << endl;
}

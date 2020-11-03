#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;


FusionEKF::FusionEKF() {
  
  // Set flag for initialization with first sensor measurements
  is_initialized_ = false;

  // Create the EKF object
  ekf_ = KalmanFilter(); 

  // Start counting the time
  previous_timestamp_ = 0;

  // Create tools object
  tools = Tools(); 

  // Initialize measurement variables 
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  H_radar_jacobian_ = MatrixXd(3, 4); 
  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225; 
  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
  // measurement matrix H laser
  H_laser_ << 1, 0, 0, 0, 
              0, 1, 0, 0; 
  // measurement matrix H radar 
  H_radar_jacobian_ << 0, 0, 0, 0, 
                       0, 0, 0, 0,
                       0, 0, 0, 0; 
}


FusionEKF::~FusionEKF() {}


void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  
  // I N I T I A L I Z A T I O N 
  // Initialize Kalman Filter once first measurement arrives
  if (!is_initialized_) {
    
    // Set timestamp
    previous_timestamp_ = measurement_pack.timestamp_; 
    
    // Initialize x
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    // Initialize P (x and y positions are considered to be more certain)
    ekf_.P_ = MatrixXd(4, 4); 
    ekf_.P_ << 1.0, 0, 0, 0, 
               0, 1.0, 0, 0, 
               0, 0, 100.0, 0,
               0, 0, 0, 100.0; 

    // Set first measurement as initial state values
    // Distinguish lidar and radar measurements
    // x and y velocities are assumed to be 1
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      float rho = measurement_pack.raw_measurements_[0];  
      float phi = measurement_pack.raw_measurements_[1]; 
      float px = rho * cos(phi); 
      float py = rho * sin(phi);
      ekf_.x_ << px, py, 1, 1;   
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      float px = measurement_pack.raw_measurements_[0]; 
      float py = measurement_pack.raw_measurements_[1]; 
      ekf_.x_ << px, py, 1, 1; 
    }

    // done initializing, no need to predict or update
    is_initialized_ = true; 
    return;
  }

  // P R E D I C T I O N  
  // Get deltaT and update current time stamp for next iteration
  long long dT_in_micro_sec = measurement_pack.timestamp_ - previous_timestamp_; 
  double dT = dT_in_micro_sec * 1e-6; 
  previous_timestamp_ = measurement_pack.timestamp_;
  // Set (i.e. update) the Transition matrix F (which depends on delta t)
  ekf_.F_ = MatrixXd(4, 4); 
  ekf_.F_ << 1, 0, dT, 0, 
             0, 1, 0, dT, 
             0, 0, 1, 0, 
             0, 0, 0, 1; 
  // Set (i.e. update) the Process noise matrix Q (which depends on delta t)
  double acc_x = 9.0; 
  double acc_y = 9.0; 
  double dT_4 = dT*dT*dT*dT;
  double dT_3 = dT*dT*dT; 
  double dT_2 = dT*dT; 
  ekf_.Q_ = MatrixXd(4, 4); 
  ekf_.Q_ << (dT_4/4.0)*acc_x, 0, (dT_3/2.0)*acc_x, 0, 
             0, (dT_4/4.0)*acc_y, 0, (dT_3/2.0)*acc_y, 
             (dT_3/2.0)*acc_x, 0, dT_2*acc_x, 0, 
             0, (dT_3/2.0)*acc_y, 0, dT_2*acc_y; 
  // Call predict function
  ekf_.Predict();

  // U P D A T E 
  // Radar measurement  
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Get the jacobian for radar measurments 
    H_radar_jacobian_ = tools.CalculateJacobian(ekf_.x_); 
    // Set radar related EKF matrices
    ekf_.R_ = R_radar_; 
    ekf_.H_ = H_radar_jacobian_;
    // Call corresponding update function
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);  
  // Lidar measurement
  } else {
    // Set lidar related EKF matrices
    ekf_.R_ = R_laser_; 
    ekf_.H_ = H_laser_;
    // Call corresponding update function
    ekf_.Update(measurement_pack.raw_measurements_); 
  }
  
  // print the outcome of the EKF iteration
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

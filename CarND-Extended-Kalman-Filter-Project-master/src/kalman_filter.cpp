#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;


KalmanFilter::KalmanFilter() {}


KalmanFilter::~KalmanFilter() {}


void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}


void KalmanFilter::Predict() {
  // Implement Kalman Filter equations for Prediction
  x_ = F_*x_; 
  P_ = F_* P_*F_.transpose() + Q_;  
}


void KalmanFilter::Update(const VectorXd &z) {
  // Implement Kalman Filter equations for LIDAR updates
  VectorXd y = z - H_*x_; 
  MatrixXd S = H_*P_*H_.transpose() + R_; 
  MatrixXd K = P_*H_.transpose()*S.inverse(); 
  MatrixXd I = MatrixXd::Identity(4, 4); 
  x_ = x_ + K*y; 
  P_ = (I - K*H_)*P_; 
}


void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // Transform state values px, py, vx, vy into radar based values rho, phi, rho_dot
  // Check against division by near zero
  float px = x_[0]; 
  float py = x_[1]; 
  float vx = x_[2];
  float vy = x_[3]; 
  float rho = sqrt(px*px + py*py); 
  float phi = atan2(py, px);
  float rho_dot; 
  if (rho > 0.01){
    rho_dot = (px*vx + py*vy)/rho; 
  }
  else{
    rho_dot = 0; 
  }
  VectorXd h_of_x = VectorXd(3);
  h_of_x << rho, phi, rho_dot;  

  // Implement Extended Kalman Filter equations for RADAR updates
  VectorXd y = z - h_of_x; 
  // Update phi value of residual if it's out of bounce, i.e greater or less +pi and -pi
  float y_phi = y[1];
  while( y_phi < - M_PI ){
    y_phi += M_PI; 
  }
  while( y_phi >  M_PI ){
    y_phi -= M_PI; 
  }
  // Write corrected phi value back into residual vector
  y(1) = y_phi; 
  // Continue with remaining Kalman Filter equations
  MatrixXd S = H_*P_*H_.transpose() + R_; 
  MatrixXd K = P_*H_.transpose()*S.inverse(); 
  MatrixXd I = MatrixXd::Identity(4, 4); 
  x_ = x_ + K*y; 
  P_ = (I - K*H_)*P_;
}

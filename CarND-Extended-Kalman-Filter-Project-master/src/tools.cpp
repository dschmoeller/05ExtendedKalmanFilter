#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout; 
using std::endl; 


Tools::Tools() {}


Tools::~Tools() {}


VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
   // Create rmse vector
   VectorXd rmse = VectorXd(4);
   rmse << 0, 0, 0, 0; 
   // Loop over estimation and ground truth values
   // Accumulate sqaured residuals
   for (unsigned int i=0; i < estimations.size(); i++){
      VectorXd residual = estimations[i] - ground_truth[i]; 
      residual = residual.array()*residual.array();
      rmse += residual;  
   }
   // Calculate mean 
   rmse = rmse/estimations.size(); 
   // Calculate squared root
   rmse = rmse.array().sqrt(); 
   // Return averaged root means squared error
  return rmse;  
}


MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
   // Store state values and terms which are used multiple times 
   float px = x_state[0]; 
   float py = x_state[1];
   float vx = x_state[2]; 
   float vy = x_state[3];
   float xx_plus_yy = px*px + py*py; 
   float sqrt_xy = sqrt(xx_plus_yy);  
   float sqrt_xy_cub = sqrt_xy*sqrt_xy*sqrt_xy; 
   MatrixXd Hj = MatrixXd(3, 4);

   // Check division by (near) zero 
   // If division by zero would apply, return the null matrix as jacobian 
   // --> This results in not updating the state covariance matrix P 
   if ((px == 0 && py == 0) || (xx_plus_yy < 0.01)) {
      Hj << 0, 0, 0, 0, 
            0, 0, 0, 0, 
            0, 0, 0, 0; 
      return Hj; 
   } 
   // Calculate and return the Jacobian matrix
   Hj << px/(sqrt_xy), py/(sqrt_xy), 0, 0, 
         -py/(xx_plus_yy), px/(xx_plus_yy), 0, 0, 
         (py*(vx*py - vy*px))/(sqrt_xy_cub), (px*(vy*px - vx*py))/(sqrt_xy_cub), px/(sqrt_xy), py/(sqrt_xy); 
   return Hj; 
}


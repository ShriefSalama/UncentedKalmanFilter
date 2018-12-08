#include "ukf.h"
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

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.0175;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.1;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  //set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;

  //define spreading parameter
  lambda_ = 3 - n_aug;

  //create sigma point matrix
  Xsig_pred_ = Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

  weights_ = VectorXd(2*n_aug_+1);

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
	  //create augmented mean vector
	  VectorXd x_aug = VectorXd(n_aug_);

	  //create augmented state covariance
	  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

	  //create sigma point matrix
	  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	  //create augmented mean state
	    x_aug.head(5) = x_;
	    x_aug(5) = 0;
	    x_aug(6) = 0;

	    //create augmented covariance matrix
	    P_aug.fill(0.0);
	    P_aug.topLeftCorner(5,5) = P_;
	    P_aug(5,5) = std_a_*std_a_;
	    P_aug(6,6) = std_yawdd_*std_yawdd_;

	    //create square root matrix
	    MatrixXd L = P_aug.llt().matrixL();

	    //create augmented sigma points
	    Xsig_aug.col(0)  = x_aug;
	    for (int i = 0; i< n_aug_; i++)
	    {
	      Xsig_aug.col(i+1)       = x_aug + sqrt(lambda+n_aug_) * L.col(i);
	      Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda+n_aug_) * L.col(i);
	    }
/////////////////////////////////////////////////////////////////////////
	    ///////////// END of Generating Sigma Points ////////////////
//////////////////////////////////////////////////////////////////////////
	    for(unsigned i = 0; i<2*n_aug+1 ; i++)
	    {
	      double px         = Xsig_aug.col(i)(0);
	      double py         = Xsig_aug.col(i)(1);
	      double speed      = Xsig_aug.col(i)(2);
	      double angle      = Xsig_aug.col(i)(3);
	      double angle_rate = Xsig_aug.col(i)(4);
	      double noise_a    = Xsig_aug.col(i)(5);
	      double noise_yaw  = Xsig_aug.col(i)(6);

	      if(angle_rate <= 0.001)
	      {
	          Xsig_pred_.col(i)(0) =  px + (delta_t * speed * cos(angle)) + (0.5 * delta_t * delta_t * cos(angle) *noise_a);
	          Xsig_pred_.col(i)(1) =  py + (delta_t * speed * sin(angle)) + (0.5 * delta_t * delta_t * sin(angle) *noise_a);
	      }
	      else
	      {
	          Xsig_pred_.col(i)(0) = px + ((speed / angle_rate) * (sin(angle+angle_rate*delta_t) - sin(angle))) + (0.5 * delta_t * delta_t * cos(angle) *noise_a);
	          Xsig_pred_.col(i)(1) = py + ((speed / angle_rate) * (cos(angle) - cos(angle+angle_rate*delta_t))) + (0.5 * delta_t * delta_t * sin(angle) *noise_a);
	      }

	      Xsig_pred_.col(i)(2) = speed + delta_t * noise_a;
	      Xsig_pred_.col(i)(3) = angle + delta_t * angle_rate + (0.5 * delta_t * delta_t * noise_yaw);
	      Xsig_pred_.col(i)(4) = angle_rate +delta_t * noise_yaw;
	    }
/////////////////////////////////////////////////////////////////////////
	    ///////////// END of Predicting Sigma Points ////////////////
/////////////////////////////////////////////////////////////////////////

	       //set weights
	       weights_(0) = lambda_/(lambda_+n_aug_);
	       for (unsigned int i = 1; i<2*n_aug_ +1; i++)
	       {
	           weights_(i) = 1/(2*(lambda_+n_aug_));
	       }

	       //predict state mean
	       for (unsigned int i = 0; i<2*n_aug_+1; i++)
	       {
	           x_ += weights_(i)*Xsig_pred_.col(i);
	       }

	       //predict state covariance matrix
	       for (unsigned int i = 0; i<2*n_aug+1; i++)
	       {
	           VectorXd C = Xsig_pred_.col(i) - x;

	           while (C(3) > M_PI || C(3) < -M_PI )
	           {
	         	if (C(3) > M_PI) C(3) -= 2*M_PI;
	         	else if (C(3) < -M_PI) C(3) += 2*M_PI;
	           }

	           P_ += weights_(i)*C*C.transpose();
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

	//SSalama , Define R , and make sure of H
	  long x_size = x_.size();
	  MatrixXd I = MatrixXd::Identity(x_size, x_size);

	  MatrixXd H_input(2,5);
	  H_input << 1,0,0,0,0,
			  	 0,1,0,0,0;

	  //add measurement noise covariance matrix
	  MatrixXd R_laser = MatrixXd(2,2);
	  R_laser <<  std_laspx_*std_laspx_, 0,
	          	  0, std_laspy_*std_laspy_;

	  MatrixXd z_pred = H_input * x_;
	  VectorXd s = H_input * P_ * H_input.transpose() + R_laser ;
	  VectorXd k = P_ * H_input.transpose() * s.inverse();
	  VectorXd y = z - z_pred ;

	  //new state
	  x_ = x_ + (k*y);
	  P_ = (I - k* H_input) * P_;

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
	  //set measurement dimension, radar can measure r, phi, and r_dot
	  int n_z = 3;

	  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

	      // extract values for better readibility
	      double p_x = Xsig_pred_(0,i);
	      double p_y = Xsig_pred_(1,i);
	      double v   = Xsig_pred_(2,i);
	      double yaw = Xsig_pred_(3,i);

	      double v1 = cos(yaw)*v;
	      double v2 = sin(yaw)*v;

	      // measurement model
	      Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
	      Zsig(1,i) = atan2(p_y,p_x);                                 //phi
	      Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
	    }
/////////////////////////////////////////////////////////////////////////
	/// END of Predicting Sigma points to Measurement space ///
/////////////////////////////////////////////////////////////////////////

	  //mean predicted measurement
	  VectorXd z_pred = VectorXd(n_z);
	  z_pred.fill(0.0);
	  for (int i=0; i < 2*n_aug_+1; i++) {
	      z_pred = z_pred + weights_(i) * Zsig.col(i);
	  }
	  //innovation covariance matrix S
	  MatrixXd S = MatrixXd(n_z,n_z);
	  S.fill(0.0);
	  for (int i = 0; i < 2 * n_aug_ + 1; i++)
	  {  //2n+1 simga points
	    //residual
	    VectorXd z_diff = Zsig.col(i) - z_pred;

	    //angle normalization
        while (z_diff(1) > M_PI || z_diff(1) < -M_PI )
        {
      	if (z_diff(1) > M_PI) z_diff(1) -= 2*M_PI;
      	else if (z_diff(1) < -M_PI) z_diff(1) += 2*M_PI;
        }

	    S = S + weights_(i) * z_diff * z_diff.transpose();
	  }

	  //add measurement noise covariance matrix
	  MatrixXd R_rader = MatrixXd(n_z,n_z);
	  R_rader <<  std_radr_*std_radr_, 0, 0,
	          	  0, std_radphi_*std_radphi_, 0,
				  0, 0,std_radrd_*std_radrd_;
	  S = S + R_rader;

/////////////////////////////////////////////////////////////////////////
	  /// END of Calculating S, Z_pred ///
/////////////////////////////////////////////////////////////////////////
	  //create matrix for cross correlation Tc
	  MatrixXd Tc = MatrixXd(n_x_, n_z);

	  //calculate cross correlation matrix
	  Tc.fill(0.0);
	  for (unsigned int i = 0; i < 2*n_aug_+1; i++)
	  {
	     VectorXd C = Zsig.col(i) - z_pred;
	        while (C(1) > M_PI || C(1) < -M_PI )
	        {
	      	if (C(1) > M_PI) C(1) -= 2*M_PI;
	      	else if (C(1) < -M_PI) C(1) += 2*M_PI;
	        }

	     VectorXd L = Xsig_pred.col(i) - x;
	        while (L(3) > M_PI || L(3) < -M_PI )
	        {
	      	if (L(3) > M_PI) L(3) -= 2*M_PI;
	      	else if (L(3) < -M_PI) L(3) += 2*M_PI;
	        }

	     Tc += weights_(i) * L * C.transpose();
	  }
	  //calculate Kalman gain K;
	  MatrixXd K = Tc*S.inverse();

	  VectorXd z_diff = z - z_pred;
      while (z_diff(1) > M_PI || z_diff(1) < -M_PI )
      {
    	if (z_diff(1) > M_PI) z_diff(1) -= 2*M_PI;
    	else if (z_diff(1) < -M_PI) z_diff(1) += 2*M_PI;
      }

	  //update state mean and covariance matrix
	  x_ = x_ + K*z_diff;
	  P_ = P_ - K*S*K.transpose();

}

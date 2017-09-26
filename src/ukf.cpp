#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
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
  std_a_ = 4;            // Change later, started with 30

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;      // Change later, started with 30

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // The measurement noise values should not be changed; these are provided by the sensor manufacturer.
  /*

  The values for the process noise std_a_ and std_yawdd_ are both set to 30. These will need to be adjusted in order to 
  get your Kalman filter working. Think about what a standard deviation of 30 means. For a Gaussian distribution, we expect 
  the acceleration to be between [−60,+60]m/s^2, or [-60,+60]rad/s^2 95% of the time.
  Too high! To put those values in perspective, the fastest measured linear acceleration for a street legal sports 
  car is currently 0 to 60 mph in 2.2 seconds. 0 to 60 mph in 2.2 seconds is about 12 m/s^2. The bike simulation probably tends 
  to have even lower acceleration.
  */

  /*
  Pay special attention to how you initialize x and P. For more information go back to the unscented Kalman filter lectures   
  notes titled "What to Expect from the Project".
  */

  const double pi = 3.14159265358979323846;

  is_initialized_ = false;

  n_x_ = 5;             // Number of state dimensions

  n_aug_ = n_x_ + 2;    // Two extra dimensions to account for nu_a and nu_psi_dot_dot

  lambda_ = 3 - n_x_;   // Best practice for lambda is 3-n_x


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

  // Gets called from main.cpp

  if (!is_initialized_) {

    x_ << 0,0,0,0,0;          // Using CTRV model -- 5 dimensions: px, py, v, psi, psi-dot (yaw rate)
    P_ << 1, 0, 0, 0, 0,      // 5x5 dimensions
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

  

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {   // Initiate RADAR state

      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];

      float px = rho * cos(phi);
      float py = rho * sin(phi);

      x_(0) = px;     // Initialized vx, vy with 0
      x_(1) = py;


    }
    else {    // Initiate LASER state

      // Set initial state and zero velocity directly from sensor, since it's laser
      // Load just px, py, and remainder=0 since it's a laser 
      
      x_(0) = meas_package.raw_measurements_[0];      // Initialized vx, vy with 0
      x_(1) = meas_package.raw_measurements_[1];

    }

  
    time_us_ = meas_package.timestamp_;     // Initialize current time
    is_initialized_ = true;                     
    return; // done (!is_initialized)

  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;   // dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    if (use_radar_==true) { 
      cout << "Use_radar is true" << endl;
      UpdateRadar(meas_package);    // Update for radar
    } 
    else {
      cout << "Skip radar update" << endl;
    }
  }
  else {  // Laser
    if (use_laser_==true) {
      cout << "Use_laser is true" << endl;
      UpdateLidar(meas_package);  // Update for lidar
    }
    else {
      cout << "Skip laser update" << endl;
    }
  
  }
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

  
  /*****************************************************************************
   *  Create augmented sigma points
   ****************************************************************************/

  // (i.e., P_ augmented with 2x2 noise process vector Q)
  // where Q = [[std_a^2, 0], [0, std_yawdd^2]]
  // Lesson 18

  VectorXd x_aug = VectorXd(n_aug_);                      // Create augmented mean vector
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);              // Create augmented state covariance
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);   // Create augmented sigma point matrix


  x_aug.head(5) = x_;                     // create augmented mean state
  x_aug(5) = 0;
  x_aug(6) = 0;

  P_aug.fill(0.0);                        // create augmented state covariance 
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  MatrixXd L = P_aug.llt().matrixL();     // create square root matrix

  Xsig_aug.col(0)  = x_aug;               // create augmented sigma points
  for (int i = 0; i< n_aug_; i++) {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }


  /*****************************************************************************
   *  Predict sigma points
   ****************************************************************************/

  // Process Prediction:
  // Store every predicted sigma point is a column of a matrix (calligraphic X_a,k|k matrix 7x15), where
  // each column is the augmented state in 7 dimensions.
  // Then run each point through process model function   X_k+1 = f(xk_,mu), where mu = process noise
  // Yields predicted X_k+1|k matrix 5x15....
  // That is:
  // Input to process model is 7 dimension augmented state vector;
  // Output is 5 dimension predicted state vector. 

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);    // Create matrix with predicted sigma points as columns
  
  for (int i = 0; i< 2 * n_aug_ + 1; i++)       // Predict sigma points through process model equations
  {
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }


  /*****************************************************************************
   *  Predict state mean and covariance
   ****************************************************************************/

  // Using predicted sigma points Xsig_pred_
  // See Lesson 23 for equations for weights, predicted mean, and predicted covariance.
  // See Lesson 24 for code

  
  // Set weights

  weights_ = VectorXd(2 * n_aug_ + 1);

  double weight_0 = lambda_/(lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {            // 2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  // predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {    // iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
  
  // predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {    // iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
    //vangle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;

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


  // Code is from EKF Project (Project 1)

  /*****************************************************************************
   *  Predict
   ****************************************************************************/


  VectorXd z = VectorXd(2);             // z is measurement package -- 2 inputs for laser (px, py)
  z = meas_package.raw_measurements_;

  MatrixXd H = MatrixXd(2, 5);
  H << 1, 0, 0, 0, 0,  // From lesson 5 section 10 but with 5th dimension since state space has 5 dimensions
      0, 1, 0, 0, 0;
  
  MatrixXd R_laser = MatrixXd(2, 2);
  R_laser << 0.0225, 0,
              0, 0.0225;

  VectorXd z_pred = H * x_;    // Note: z is from measurement pack
  VectorXd y = z - z_pred;   

  // R_ is set to R_laser_ in FusionEKF::ProcessMeasurement() 
  MatrixXd S = H * P_ * H.transpose() + R_laser; 
  MatrixXd K =  P_ * H.transpose() * S.inverse();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;

  cout << "After Lidar Update................" << endl;
  cout << "P_ = " << P_ << endl;
  cout << "x_ = " << x_ << endl;


  /*****************************************************************************
   *  Calculate Laser NIS (Epsilon)
   ****************************************************************************/

  double NIS_laser_ = (z - z_pred).transpose() * S.inverse() * (z - z_pred);  // epsilon -- scalar value

  cout << "NIS_laser_ = " << NIS_laser_ << endl;


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

  You'll also need to calculate the raxdar NIS.
  */


  int n_z = 3;                            // num dimensions in z measurement
  VectorXd z = VectorXd(n_z);             // z is measurement package with n_z dimensions (3 for radar)
  z = meas_package.raw_measurements_;


  /*****************************************************************************
   *  Predict radar sigma points
   ****************************************************************************/

  // Measurement prediction:
  // Use sigma points we already generated
  // Measurement model is linear, b/c noise is additive, unlike the non-linear
  // noise process in the process model, so we don't need to go through state augmentation again.
  // Just transform individual sigma points into measurement space (rho, phi, rho_dot), and use
  // points to predict mean z_k+1|k and covariance matrix S_k+1|k of measurement.
  // Go from predicted sigma points (5x15) via measurement model to measured sigma points (3x15)
  // See Lesson 27

  /*
  // set vector for weights
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }
  */


  // create matrix for sigma points in measurement space
  
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {        // 2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;


  /*****************************************************************************
   *  Perform Radar Update
   ****************************************************************************/

  // Lesson 30

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {        // 2n+1 simga points

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;
  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  
  cout << "After Radar Update................" << endl;
  cout << "P_ = " << P_ << endl;
  cout << "x_ = " << x_ << endl;


  /*****************************************************************************
   *  Calculate Radar NIS (Epsilon)
   ****************************************************************************/

  // Lesson 31

  double NIS_radar_ = (z - z_pred).transpose() * S.inverse() * (z - z_pred);  // epsilon -- scalar value

  cout << "NIS_radar_ = " << NIS_radar_ << endl;


}

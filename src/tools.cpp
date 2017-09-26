#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

  // Lifted from CalculateRMSE() from Project 1

	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check validity of inputs
	if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
		cout << "Invalid estimation or ground truth data." << endl;
		return rmse;
	}

	// Accumulate squared residuals
	for (unsigned int i=0; i < estimations.size(); i++) {

		VectorXd residual = estimations[i] - ground_truth[i];

		// coefficient-wise multiplication
		residual = residual.array() * residual.array();
		rmse += residual;
	}

	// Calculate the mean
	rmse = rmse/estimations.size();

	// Sqrt
	rmse = rmse.array().sqrt();

	// Done
	if (rmse(0) > 0.11 || rmse(1) > 0.11) {
		cout << "x/y RMSE > 0.11: " << rmse(0) << ", " << rmse(1) << endl;
	}
	if (rmse(2) > 0.52 || rmse(3) > 0.52) {
		cout << "vx/vy RMSE > 0.52: " << rmse(2) << ", " << rmse(3) << endl;
	}

	return rmse;


}
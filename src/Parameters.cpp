/*
 * Parameters.cc
 *
 *  Created on: Apr 26, 2014
 *      Author: wilfeli
 */

#include <math.h>
#include "Parameters.h"



Parameters::Parameters(){
	//production function parameters vector of length 2
//	F_F_theta(1,2);
	F_F_theta = Eigen::MatrixXd::Zero(1,2);
	F_F_theta << 1.0, 1.0;
	//utility function parameters
//	GOAL_T_theta(1,2);
	GOAL_T_theta = Eigen::MatrixXd::Zero(1,2);
	GOAL_T_theta << 3.0, 0.5;


};

void
Parameters::init(){
	//create vector of betas
	BETA_T = Eigen::MatrixXd::Zero(T_MAX, 1);

	//create vector of betas
	for (int i=0; i<T_MAX; i++){
		BETA_T(i,0) = pow(BETA, i);
	};
};


ParametersW::ParametersW(){
	//production function parameters vector of length 2
//	F_F_F_theta(1,2);
	F_F_F_theta = Eigen::MatrixXd::Zero(1,2);
	F_F_F_theta << 1.0, 1.0;
	//utility function parameters
//	H_GOAL_T_theta(1,2);
	H_GOAL_T_theta = Eigen::MatrixXd::Zero(1,2);
	H_GOAL_T_theta << 3.0, 0.5;
};

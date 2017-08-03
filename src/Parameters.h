/*
 * Parameters.h
 *
 *  Created on: Apr 25, 2014
 *      Author: wilfeli
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <Eigen/Dense>
#include <string>

class Parameters{
public:
	Parameters();

	int T_MAX = 20; //maximum length of forecasting
    int ADP_BF_N = 1; //dimensions of basis functions
    int ADP_API_LM_N = 5;
    int ADP_API_LM_M = 5;
	std::string GRID = "small"; //type of a grid
	double DIV_SHARE = 0.5; //share of dividends to pay
    std::string ACCOUNTING_TYPE  = "classic"; //type of accounting
	double BETA = 0.95; //discounting rate
	double WM_LENGTH = std::numeric_limits<double>::infinity(); //length of wm
	double wm_INSIDE_UPDATE = 0.0; //update inside opt search
    int opt_CS_N = 10;//number of random draws for runs
    std::string opt_TYPE = "EO_CS"; //type of optimization algorithm (could be EO_ADP, QL, mRE)
    double QL_L = 0.05;//learning speed for QL
    double QL_experimental_rate = 0.1; //experimental rate
    double mRE_EPSILON = 0.95; //0.05; //adjustment for mRE algorithm
    double mRE_PHI = 0.95; //learning speed
    double mRE_T = 1.0; //cooling parameter


	Eigen::MatrixXd F_F_theta; //parameters of production function
	Eigen::MatrixXd GOAL_T_theta; //parameters of utility
	Eigen::MatrixXd BETA_T; //vector of betas for time discounting

	void init(); //initialize service parameters

};


class ParametersW{
public:
	ParametersW();

	//money on hand
	double F_asCBDt0_Q = 100.0;
	double H_asCBDt0_Q = 1.0;

	double F_wm_LENGTH = std::numeric_limits<double>::infinity();
	double H_wm_LENGTH = std::numeric_limits<double>::infinity();

	int F_T_MAX = 20;
	int H_T_MAX = 20;
    
    int F_opt_CS_N = 10;
    int H_opt_CS_N = 10;

	std::string F_GRID = "small";
	std::string H_GRID = "small";
    
    std::string F_opt_TYPE = "EO_CS";
    std::string H_opt_TYPE = "EO_CS";
    
    std::string SIMULATION_MODE = "test";
    std::string F_ACCOUNTING_TYPE  = "classic";

	double F_wm_INSIDE_UPDATE = 0.0; //not implemented in the code
	double H_wm_INSIDE_UPDATE = 0.0; //not implemented in the code

    double F_QL_L = 0.95;//learning speed for QL
    double H_QL_L = 0.95;//learning speed for QL

    double F_mRE_EPSILON = 0.95;//0.05; //adjustment for mRE algorithm
    double F_mRE_PHI = 0.95; //learning speed
    double F_mRE_T = 1.0; //cooling parameter
    double H_mRE_EPSILON = 0.95;//0.05; //adjustment for mRE algorithm
    double H_mRE_PHI = 0.95; //learning speed
    double H_mRE_T = 1.0; //cooling parameter

	double W_SEED = 2013;

	Eigen::MatrixXd F_F_F_theta;
	Eigen::MatrixXd H_GOAL_T_theta;


};

#endif /* PARAMETERS_H_ */

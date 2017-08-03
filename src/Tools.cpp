/*
 * Tools.cpp
 *
 *  Created on: Apr 22, 2014
 *      Author: wilfeli
 */

#include <vector>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <Eigen/Dense>
#include <iostream>
#include "Tools.h"

using Eigen::MatrixXd;
using std::vector;

namespace Tools{

void
mgrid(MatrixXd grid, MatrixXd* meshgrid){
	//creates meshgrid from vector of points
	//define number of steps for each input
	long N_rows = grid.rows();
	int n_steps_i;
	double n_steps_i_int;
	int j_val;

	double val;
	vector<int> n_steps;
	vector<int> ac_n_steps;
	vector<vector<double>> values;
	vector<double> values_temp;
	ac_n_steps.push_back(1);


    //create list of values
	for (int i=0; i < N_rows; i++){
		//n_steps
		n_steps_i_int = (grid(i,1) - grid(i,0))/grid(i,2);

		n_steps_i = static_cast<int>(floor(n_steps_i_int)) + 1;
		n_steps.push_back(n_steps_i);
		ac_n_steps.push_back(ac_n_steps.back() * n_steps_i);
        
		for (int j=0; j<n_steps_i; j++){
			val = grid(i,0) + grid(i,2)*j;
			values_temp.push_back(val);
		};

		values.push_back(values_temp);
		values_temp.clear();
	};


	//create matrix and fill it with values
	(*meshgrid) = MatrixXd::Zero(N_rows, ac_n_steps.back());

	int step_length = ac_n_steps.back();
	for (int i = 0; i < N_rows; i++){
		step_length = step_length / n_steps[i];
		for (int j = 0; j < ac_n_steps.back(); j++){
			j_val = ((int)(j/step_length))%n_steps[i];
			(*meshgrid)(i,j) = values[i][j_val];

		};
	};

//	std::cout << (*meshgrid).transpose() << "\n";

};
    
void
mgrid_test(MatrixXd grid, MatrixXd* meshgrid){
    //creates meshgrid from vector of points
    //reverse to correspond to python
    long N_rows = grid.rows();
    int n_steps_i;
    double n_steps_i_int;
    int j_val;
    
    double val;
    vector<int> n_steps;
    vector<int> ac_n_steps;
    vector<vector<double>> values;
    vector<double> values_temp;
    ac_n_steps.push_back(1);
    
    
    
    for (int i=0; i < N_rows; i++){
        //n_steps
        n_steps_i_int = (grid(i,1) - grid(i,0))/grid(i,2);
        
        n_steps_i = static_cast<int>(floor(n_steps_i_int)) + 1;
        n_steps.push_back(n_steps_i);
        ac_n_steps.push_back(ac_n_steps.back() * n_steps_i);
        
        for (int j=0; j<n_steps_i; j++){
            val = grid(i,0) + grid(i,2)*j;
            values_temp.push_back(val);
        };
        
        values.push_back(values_temp);
        values_temp.clear();
    };
    
    
    //create matrix and fill it with values
    (*meshgrid) = MatrixXd::Zero(N_rows, ac_n_steps.back());
    
    int step_length = 1;
    for (int i = 0; i < N_rows; i++){
        
        for (int j = 0; j < ac_n_steps.back(); j++){
            j_val = ((int)(j/step_length))%n_steps[i];
            (*meshgrid)(i,j) = values[i][j_val];
            
        };
        step_length = step_length * n_steps[i];
    };
    
//    std::cout << (*meshgrid).transpose() << "\n";
    
};
    




void
print_vector(std::vector<double> vec){
	//print all q
	for (auto iter:vec){
		std::cout << iter << " ";
	};
    
    std::cout << std::endl;
};
    
void
print_vector(std::vector<std::vector<double>> vec){
    //print all q
    for (auto iter:vec){
        print_vector(iter);
    };
    
};
    
    

MyRNG::MyRNG(double seed_):state(seed_){
    m_w = seed_;
};
    
// Produce a uniform random sample from the open interval (0, 1).
// The method will not return either end point.
double
MyRNG::GetUniform(){
    // 0 <= u < 2^32
	uint64_t u = GetUint();
    // The magic number below is 1/(2^32 + 2).
    // The result is strictly between 0 and 1.
    return (u + 1.0) * 2.328306435454494e-10;
};
    
// This is the heart of the generator.
// It uses George Marsaglia's MWC algorithm to produce an unsigned integer.
// See http://www.bobwheeler.com/statistics/Password/MarsagliaPost.txt

uint64_t
MyRNG::GetUint(){
    m_z = 36969 * (m_z & 65535) + (m_z >> 16);
    m_w = 18000 * (m_w & 65535) + (m_w >> 16);
    return (m_z << 16) + m_w;
};

    
    
double
get_normal(double mean, double sigma, MyRNG& rng){
    if(rng.hasSpare)
	{
		rng.hasSpare = false;
		return sigma * sqrt(rng.rn1) * sin(rng.rn2) + mean;
	}
    
	rng.hasSpare = true;
    
	rng.rn1 = rng.GetUniform();
	if(rng.rn1 < 1e-100) {rng.rn1 = 1e-100;};
	rng.rn1 = -2 * log(rng.rn1);
	rng.rn2 = rng.GetUniform() * M_PI * 2;
    
	return sigma*sqrt(rng.rn1) * cos(rng.rn2) + mean;
    

        
};
    
int
get_int(int min, int max, MyRNG& rng){
	uint64_t x_raw = rng.GetUint();
    //clean state of rng to simplify code
    rng.hasSpare = false;
    
    int x = min + (x_raw % (int)(max - min + 1));
    
    return x;
    
    
};

};

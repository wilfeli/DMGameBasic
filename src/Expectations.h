/*
 * Expectations.h
 *
 *  Created on: Apr 24, 2014
 *      Author: wilfeli
 */

#ifndef EXPECTATIONS_H_
#define EXPECTATIONS_H_

#include <vector>
#include <map>

class ExpectationBackward{
public:
	ExpectationBackward(double _mu, int _n, double _v): mu(_mu), n(_n), v(_v){
		//push p, n into s_t_1 as backward ? data
		s_t_1.insert(std::pair<std::string, double> ("mu",_mu));
		s_t_1.insert(std::pair<std::string, double> ("n",_n));
		s_t_1.insert(std::pair<std::string, double> ("v",_v));
	};

	double mu;
	int n;
	double v;

	std::vector<std::vector<double>> s_t;
	std::map<std::string, double> s_t_1;



};






#endif /* EXPECTATIONS_H_ */

/*
 * MarketTools.h
 *
 *  Created on: Apr 14, 2014
 *      Author: wilfeli
 */

#ifndef MARKETTOOLS_H_
#define MARKETTOOLS_H_

#include <Eigen/Dense>
#include <list>
#include <tuple>
#include <vector>
#include <deque>

#include "Agent.h"
class IBuyer;
class ISeller;


class Bid{
public:
	Bid(double price = 0.0, double quantity = 0.0);

public:
	double p;
	double q;
	double pq;
	double pq_share;

	double pq_t_0; //initial value
	double p_t_0; //initial value
	double q_t_0; //initial value
	double pq_share_t_0; //initial value

	IBuyer* issuer;
	Eigen::MatrixXd demand_curve; //p,q structure

};


class Ask{
public:
	Ask(double price = 0.0, double quantity = 0.0);

public:
	double p;
	double q;
	double q_share;

	double p_t_0;  //initial value when the decision is made
	double q_t_0; //initial value when the decision is made
	double q_share_t_0; //initial value

	ISeller* issuer;
	Eigen::MatrixXd supply_curve;  //p,q structure




};


class DemandCurveBid{
public:
	double p = 0.0;
	double q = 0.0;
	std::deque<std::pair<int,Bid*>> demand_curve_point_bids; //int holds information about active row in the demand curve for this price point
};


class SupplyCurveAsk{
public:
	double p = 0.0;
	double q = 0.0;
	std::deque<std::pair<int,Ask*>> supply_curve_point_asks; //int holds information about active row in the supply curve for this price point
};


class MarketPrice{
public:
	MarketPrice(double, bool);
	double p = 0.0;
	bool FLAG_EQ = 0.0;
//	std::vector<std::pair<bool, double>> p_history; //history of market prices
};


#endif /* MARKETTOOLS_H_ */

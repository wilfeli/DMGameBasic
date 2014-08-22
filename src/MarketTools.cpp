/*
 * MarketTools.cc
 *
 *  Created on: Apr 14, 2014
 *      Author: wilfeli
 */


#include "MarketTools.h"


Bid::Bid(double price, double quantity): p(price), q(quantity), demand_curve(1,2) {
		issuer = NULL;
		demand_curve << price, quantity;
		pq = price*quantity;
};



Ask::Ask(double price, double quantity): p(price), q(quantity), supply_curve(1,2) {
		issuer = NULL;
		supply_curve << price, quantity;
};



MarketPrice::MarketPrice(double p_, bool flag_): p(p_), FLAG_EQ(flag_){};

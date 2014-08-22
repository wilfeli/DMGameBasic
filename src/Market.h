/*
 * Market.h
 *
 *  Created on: Apr 13, 2014
 *      Author: wilfeli
 */

#ifndef MARKET_H_
#define MARKET_H_

#include <vector>
#include <map>
#include <memory>
#include <iostream>
#include <string>
#include <Eigen/Dense>
#include "MarketTools.h"


class ISellerL;
class IBuyerL;
class ISellerC;
class IBuyerC;
class ContractHK;
class W;




class MarketL{
public:
	MarketL(W &);
	MarketL(W&, std::string);
	std::vector<ISellerL*> sellers;
	std::vector<IBuyerL*> buyers;
	std::vector<Bid*> bids; //list of pointers to bids
	std::vector<Ask*> asks; //list of pointers to asks

	std::map<double,DemandCurveBid*> demand_curve;
	std::map<double,SupplyCurveAsk*> supply_curve;
	void AddSellerL(ISellerL*);
	void AddBuyerL(IBuyerL*);

	void step();
	void _collect();
	void _match();
	bool _clear(ISellerL*, IBuyerL*, double, double);
	void _clear_C(ContractHK*);

	int random_shuffle(int);


	typedef decltype(supply_curve)::iterator SupplyIter;
	typedef decltype(demand_curve)::iterator DemandIter;

	SupplyIter _i_p_supply(double);
	DemandIter _i_p_demand(double);



	void _demand();
	void _supply();
	void _eq();
	std::vector<MarketPrice*> market_price;

	W& w;

 };

class MarketC : public MarketL{
public:
	MarketC(W &);
	void AddSellerC(ISellerC*);
	void AddBuyerC(IBuyerC*);
	std::vector<ISellerC*> sellers;
	std::vector<IBuyerC*> buyers;
	void _match();
	bool _clear(ISellerC*, IBuyerC*, double, double);
	void _collect();
	void step();

};






#endif /* MARKET_H_ */

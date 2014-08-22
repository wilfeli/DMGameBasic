/*
 * Market.cc
 *
 *  Created on: Apr 13, 2014
 *      Author: wilfeli
 */


#include <list>
#include <algorithm>
#include <functional>
#include <iostream>


#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "W.h"
#include "Market.h"
#include "ISeller.h"
#include "IBuyer.h"
#include "Contract.h"
#include "Message.h"
#include "Agent.h"
#include "Tools.h"
#include "Parameters.h"

//general code to find equilibrium


MarketL::MarketL(W &w_, std::string mode): w(w_){
};

MarketL::MarketL(W &w_): w(w_){
	w.AddMarketL(this);
};



void
MarketL::AddSellerL(ISellerL* agent) { sellers.push_back(agent); };

void
MarketL::AddBuyerL(IBuyerL* agent) { buyers.push_back(agent); };


void
MarketL::_collect(){
	//clears previous asks
	asks.clear();

	//collect asks
	//for each seller, call
	for(auto seller=sellers.begin();seller != sellers.end(); ++seller){
		//?collect them only once
		asks.push_back((**seller).ac_L());
	};

	bids.clear();
	for(auto buyer=buyers.begin();buyer != buyers.end(); ++buyer){
		//?collect them only once
		bids.push_back((**buyer).ac_L());
	};
};


void
MarketL::_demand(){


	//clean previous demand_curve
	for (auto dc=demand_curve.begin(); dc!=demand_curve.end(); ++dc){
		dc->second->demand_curve_point_bids.clear();
        delete dc->second;
	};
	demand_curve.clear();



	//form list of cut-off prices
	std::list<double> p_t;
	long N_steps = 0;

	typedef std::map<double,DemandCurveBid*>::iterator ITER_DC;
	typedef decltype(demand_curve.rend()) REV_ITER_DC;

	ITER_DC p_stop;
	REV_ITER_DC rev_p_stop;
	double p;

	double key;

	for (auto bid = bids.begin(); bid != bids.end(); ++bid){
		//only for active agents
//		if (reinterpret_cast<Agent>(*(**bid).issuer).status != 0){
//		if (reinterpret_cast<Agent*>((**bid).issuer)->status != 0){
		if ((**bid).issuer->get_status() != 0){
			//get size of demand curve
			N_steps = (**bid).demand_curve.rows();

			//iterate over each row
			for (int i = 0; i < N_steps; ++i){
				//check that p isn't in the map
				//key equal to p
				key = (**bid).demand_curve(i,0);
				auto dcbidin = demand_curve.find(key);
				if (dcbidin == demand_curve.end()){
					//create new entry
					DemandCurveBid* dcbidnew = new DemandCurveBid();
					dcbidnew->p = key;
					dcbidnew->q += (**bid).demand_curve(i,1);
					dcbidnew->demand_curve_point_bids.push_back(std::make_pair(i, *bid));

					//add new entry to the map and get a pointer to it
					std::pair<ITER_DC,bool>
					dcbid_pos = demand_curve.insert(std::pair<double, DemandCurveBid*>(key,dcbidnew));

					//add q to accommodate inserting point
					//create pointer to the current price for the bid
					dcbidin = dcbid_pos.first;

					if (++dcbid_pos.first != demand_curve.end()){
						dcbidnew->q += dcbid_pos.first->second->q;
						//add all other agents' bids
						dcbidnew->demand_curve_point_bids.insert(dcbidnew->demand_curve_point_bids.end(),
																dcbid_pos.first->second->demand_curve_point_bids.begin(),
																dcbid_pos.first->second->demand_curve_point_bids.end());
					};

				}else{
					dcbidin->second->q += (**bid).demand_curve(i,1);
					dcbidin->second->demand_curve_point_bids.push_back(std::make_pair(i, *bid));

				};

				//update curve for demand with q
				p_stop = demand_curve.end();
				rev_p_stop = demand_curve.rend();

				if (i > 0){
					p = (**bid).demand_curve(i-1,0);
					p_stop = demand_curve.find(p);
				}

				if (p_stop == demand_curve.end()){
					rev_p_stop = demand_curve.rend();
				}else{
					rev_p_stop = std::reverse_iterator<ITER_DC>(p_stop);
				};

				if (dcbidin == demand_curve.begin()){

				}else{
					dcbidin--;
				};



				//add q to past p until previous cut-off
				for (std::reverse_iterator<ITER_DC> rev_dcbid (dcbidin);
						rev_dcbid != rev_p_stop; ++rev_dcbid){
					rev_dcbid->second->q += (**bid).demand_curve(i,1);
					rev_dcbid->second->demand_curve_point_bids.push_back(std::make_pair(i, *bid));
				};

			};
		};
	};


//	//print all q
//	for (auto iter:demand_curve){
//		std::cout << iter.first << "," << iter.second->q << "\n";
//	};


};

void
MarketL::_supply(){

	//clean previous supply_curve
	for (auto sc=supply_curve.begin(); sc!=supply_curve.end(); ++sc){
		sc->second->supply_curve_point_asks.clear();
        delete sc->second;
	};
	supply_curve.clear();



	//form list of cut-off prices
	std::list<double> p_t;
	long N_steps = 0;

	typedef std::map<double,SupplyCurveAsk*>::iterator ITER_sc;
	typedef decltype(supply_curve.rend()) REV_ITER_sc;

	ITER_sc p_stop;
	REV_ITER_sc rev_p_stop;
	double p = 0.0;

	double key;

	for (auto ask = asks.begin(); ask != asks.end(); ++ask){
		//only for active agents
//		if (reinterpret_cast<Agent>(*(**ask).issuer).status != 0){
		if ((**ask).issuer->get_status() != 0){
			//get size of supply curve
			N_steps = (**ask).supply_curve.rows();


			//iterate over each row
			for (long i = N_steps - 1; i >= 0; --i){
				//check that p isn't in the map
				//key equal to p
				key = (**ask).supply_curve(i,0);
				auto scaskin = supply_curve.find(key);
				if (scaskin == supply_curve.end()){
					//create new entry
					SupplyCurveAsk* scasknew = new SupplyCurveAsk();
					scasknew->p = key;
					scasknew->q += (**ask).supply_curve(i,1);
					scasknew->supply_curve_point_asks.push_back(std::make_pair(i, *ask));

					//add new entry to the map and get a pointer to it
					std::pair<ITER_sc,bool>
					scask_pos = supply_curve.insert(std::pair<double, SupplyCurveAsk*>(key,scasknew));

					//add q to accommodate inserting point
					//create pointer to the current price for the ask
					scaskin = scask_pos.first;


					if (scaskin != supply_curve.begin()){
						scaskin--;

						scasknew->q += scaskin->second->q;
						//add all other agents' asks
						scasknew->supply_curve_point_asks.insert(scasknew->supply_curve_point_asks.end(),
																scaskin->second->supply_curve_point_asks.begin(),
																scaskin->second->supply_curve_point_asks.end());
					};

				}else{
					scaskin->second->q += (**ask).supply_curve(i,1);
					scaskin->second->supply_curve_point_asks.push_back(std::make_pair(i, *ask));

				};

				//update curve for supply with q
				p_stop = supply_curve.end();
				rev_p_stop = supply_curve.rend();

				if (i < N_steps - 1){
					p = (**ask).supply_curve(i + 1,0);
					p_stop = supply_curve.find(p);
				}

				if (p_stop == supply_curve.end()){
					rev_p_stop = supply_curve.rend();
				}else{
					rev_p_stop = std::reverse_iterator<ITER_sc>(p_stop);
				};

				std::reverse_iterator<ITER_sc> rev_until (scaskin);

				//move to the next element
				rev_until--;

				//add q to future p until next cut-off
				for (std::reverse_iterator<ITER_sc> rev_scask (p_stop);
						rev_scask != rev_until; ++rev_scask){
					rev_scask->second->q += (**ask).supply_curve(i,1);
					rev_scask->second->supply_curve_point_asks.push_back(std::make_pair(i, *ask));
				};

			};
		};
	};

//	//print all q
//	for (auto iter:supply_curve){
//		std::cout << iter.first << "," << iter.second->q << "\n";
//	};


};

void
MarketL::_eq(){
    bool FLAG_ASK = false;
    bool FLAG_BID = false;
	bool FLAG_EQ = false; //equilibrium
	bool FLAG_S_H_D = false; //supply higher than demand
	bool FLAG_END = false; //indicator that needs to stop iteration
	double p_eq = 0.0;

	std::vector<double> p_all;

	//make list of all prices
	for (auto ask=supply_curve.begin(); ask != supply_curve.end(); ++ask){
		p_all.push_back(ask->first);
	};

	for (auto bid=demand_curve.begin(); bid != demand_curve.end(); ++bid){
		p_all.push_back(bid->first);
	};

	std::sort(p_all.begin(), p_all.end());

	//remove duplicates
	auto p_end = std::unique(p_all.begin(), p_all.end());

	auto p_i = p_end;
	--p_i;
    
    //check that have demand and supply
    if (supply_curve.begin() != supply_curve.end()){
        FLAG_ASK = true;
    };

    if (demand_curve.begin() != demand_curve.end()){
        FLAG_BID = true;
    };

    if (FLAG_BID && FLAG_ASK){
        while (( !FLAG_EQ ) && (!FLAG_END)){

            double q_s = 0.0;
            double q_d = 0.0;

            auto d_i = _i_p_demand((*p_i));
            auto s_i = _i_p_supply((*p_i));

            if (d_i != demand_curve.end()){
                q_d = _i_p_demand((*p_i))->second->q;
            };

            if (s_i != supply_curve.end()){
                q_s = _i_p_supply((*p_i))->second->q;
            };


            if (q_s >= q_d){
                FLAG_S_H_D = true;
            }else{
                if ((FLAG_S_H_D) && (p_i != std::prev(p_end))){
                    FLAG_EQ = true;
                    FLAG_S_H_D = false;
                    p_eq = *p_i;
                    
    //                if (w.param->SIMULATION_MODE != "test"){

                        //but if supply is now zero ? - take next price
                        if (q_s <= 0.0){
                            //take previous price
                            p_eq = *std::next(p_i);
                        };
    //                };
                };
                if (p_i == std::prev(p_end)) {
                    FLAG_END = true;
                    FLAG_EQ = true; //call it equilibrium
                    FLAG_S_H_D = false;
                    p_eq = *p_i;
                };
            };

            if ((p_i == p_all.begin()) && (!FLAG_END && !FLAG_EQ)){
                if (!supply_curve.empty() && !demand_curve.empty()){
                    p_eq = (supply_curve.begin()->first + demand_curve.rbegin()->first)/2;
                    FLAG_EQ = true;
                }else{
                    p_eq = *p_all.begin();
                };
                FLAG_END = true;
            };

            --p_i;
        };
    }else{
//        p_eq = std::numeric_limits<double>::signaling_NaN();
        p_eq = 0.0;
    };

	//save price to market

	market_price.push_back(new MarketPrice(p_eq, FLAG_EQ));

//	//print all q
//	for (auto iter:supply_curve){
//		std::cout << iter.first << "," << iter.second->q << std::endl;
//	};
//
//	//print all q
//	for (auto iter:demand_curve){
//		std::cout << iter.first << "," << iter.second->q << std::endl;
//	};
//
//
//
//	std::cout << p_eq << std::endl;


};


void
MarketL::step(){

	_collect();
	_demand();
	_supply();
	_eq();

	if ((market_price.back()->FLAG_EQ) && (!supply_curve.empty() && !demand_curve.empty())) {
		_match();
	};

}

void
MarketC::step(){

	_collect();
	_demand();
	_supply();
	_eq();

	if ((market_price.back()->FLAG_EQ) && (!supply_curve.empty() && !demand_curve.empty())) {
		_match();
	};


}

void
MarketL::_match(){

	//get sellers at that price point
	//randomize them
	//iterate over them
	double price = market_price.back()->p;

	//find position for the price


    auto scask = _i_p_supply(price);
    //shuffle the suppliers
    using std::placeholders::_1;
    
    if (scask != supply_curve.end()) {
        
        if (w.param->SIMULATION_MODE != "test"){
            std::random_shuffle(scask->second->supply_curve_point_asks.begin(),
                    scask->second->supply_curve_point_asks.end(),
                    std::bind(&MarketL::random_shuffle,
                            this, std::placeholders::_1));
        }else{
            Tools::random_shuffle(scask->second->supply_curve_point_asks.begin(),
                                scask->second->supply_curve_point_asks.end(),
                                w.myrng);
        };

        for (auto ask:scask->second->supply_curve_point_asks){
            //check how much is selling
            //construct a message for an agent
            MessageMarketLCheckAsk* mes = new MessageMarketLCheckAsk (this, price);
            double q_ask = dynamic_cast<ISellerL*>(ask.second->issuer)->_mas_q_sell(mes);

            //clean memory
            delete mes;

            //find buyers for this price
            auto dcbid = _i_p_demand(price);

            if (dcbid != demand_curve.end()){
                
                
                if (w.param->SIMULATION_MODE != "test"){
                    //randomize buyers
                    std::random_shuffle(dcbid->second->demand_curve_point_bids.begin(),
                            dcbid->second->demand_curve_point_bids.end(),
                            std::bind(&MarketL::random_shuffle,
                                                this, std::placeholders::_1));
                }else{
                    //randomize buyers
                    Tools::random_shuffle(dcbid->second->demand_curve_point_bids.begin(),
                                        dcbid->second->demand_curve_point_bids.end(),
                                        w.myrng);

                };

                for (auto bid:dcbid->second->demand_curve_point_bids){
                    MessageMarketLCheckBid* mes = new MessageMarketLCheckBid (this, price);
                    double q_bid = dynamic_cast<IBuyerL*> (bid.second->issuer)->_mas_q_buy(mes);

                    //clean memory
                    delete mes;

                    if (q_bid > 0.0){
                        double q_bought = 0.0;
                        if (q_ask > 0.0){
                            if (q_bought < q_bid){

                                double q_sell = std::min(q_bid - q_bought, q_ask);

                                bool FLAG_CLEARED = _clear(dynamic_cast<ISellerL*>(ask.second->issuer),
                                                        dynamic_cast<IBuyerL*>(bid.second->issuer),
                                                    price,
                                                    q_sell);

                                if (FLAG_CLEARED){
                                    q_bought += q_sell;
                                    q_ask -= q_sell;
                                }else{
                                    break;
                                };
                            };
                        }else{
                            break;
                        };
                    };
                };
            };
        };
    };
};

bool
MarketL::_clear(ISellerL* seller, IBuyerL* buyer, double price, double q){

	bool FLAG_CLEARED = true;

	MessageMarketLSellC* mes = new MessageMarketLSellC(this, buyer, price, q);

	seller->_sell_asCHK(mes);

	delete mes;

	return FLAG_CLEARED;

};

void
MarketL::_clear_C(ContractHK* contract){

	contract->holder->_buy_asCHK(contract);

};


MarketL::SupplyIter
MarketL::_i_p_supply(double p){
	//handle no finding and last element

	MarketL::SupplyIter supply_iter ;//= supply_curve.begin();
//	bool FLAG_P = false;

	//find first price that is lower
	std::reverse_iterator<MarketL::SupplyIter> supply_rev (supply_curve.end());

	while ((supply_rev != supply_curve.rend()) && (supply_rev->first > p)){
		supply_rev++;
	};

	if ((supply_rev.base() == supply_curve.end()) && (supply_curve.end()!=supply_curve.begin())){
		supply_iter = std::prev(supply_curve.end());
	}else{
		if (supply_rev == supply_curve.rend()){
			supply_iter = supply_curve.end();
		}else{
			supply_iter = supply_rev.base();
		};
	};

	return supply_iter;
};


MarketL::DemandIter
MarketL::_i_p_demand(double p){

	//find first price that is higher


	MarketL::DemandIter demand_iter = demand_curve.begin();

	while ((demand_iter->first < p) && (demand_iter != demand_curve.end())){
		demand_iter++;
	};

	return demand_iter;
};



int
MarketL::random_shuffle(int i){
	typedef decltype(w.rng) RandomGeneratorType;

	//use world random to get generator
	boost::uniform_int<> distribution(0, i - 1);
	boost::variate_generator< RandomGeneratorType&, boost::uniform_int<> >
	    MarketL_shuffle( w.rng, distribution );

	return MarketL_shuffle();

};


MarketC::MarketC(W &w_):MarketL(w_, "derived"){
	w.AddMarketC(this);

};


void
MarketC::_collect(){
	//collect asks
	//for each seller, call
	asks.clear();
	for(auto seller=sellers.begin();seller != sellers.end(); ++seller){
		//?collect them only once
		asks.push_back((**seller).ac_C());
	};

	bids.clear();
	for(auto buyer=buyers.begin();buyer != buyers.end(); ++buyer){
		//?collect them only once
		bids.push_back((**buyer).ac_C());
	};
};



/**
 *
 * factor out who goes first and make it property of the market
 *
 *
 */

void
MarketC::_match(){

	//get sellers at that price point
	//randomize them
	//iterate over them
	double price = market_price.back()->p;

	//find position for the price
    auto scask = _i_p_supply(price);

    //find buyers for this price
	auto dcbid = _i_p_demand(price);

    using std::placeholders::_1;


    if ((scask != supply_curve.end()) && (dcbid != demand_curve.end())) {
        
        
        if (w.param->SIMULATION_MODE != "test"){
            //randomize buyers
            std::random_shuffle(dcbid->second->demand_curve_point_bids.begin(),
                    dcbid->second->demand_curve_point_bids.end(),
                    std::bind(&MarketL::random_shuffle,
                                        this, std::placeholders::_1));

            //shuffle the suppliers
            //check that have some suppliers
            std::random_shuffle(scask->second->supply_curve_point_asks.begin(),
                    scask->second->supply_curve_point_asks.end(),
                    std::bind(&MarketL::random_shuffle,
                            this, std::placeholders::_1));
        }else{
            //randomize buyers
            Tools::random_shuffle(dcbid->second->demand_curve_point_bids.begin(),
                                dcbid->second->demand_curve_point_bids.end(),
                                w.myrng);
            
            //shuffle the suppliers
            //check that have some suppliers
            Tools::random_shuffle(scask->second->supply_curve_point_asks.begin(),
                                scask->second->supply_curve_point_asks.end(),
                                w.myrng);

        };


		for (auto bid:dcbid->second->demand_curve_point_bids){
			MessageMarketCCheckBid* mes = new MessageMarketCCheckBid (this, price);
			double q_bid = dynamic_cast<IBuyerC*>(bid.second->issuer)->_mas_q_buy(mes);

			//clean memory
			delete mes;


			if (q_bid > 0.0)
			{
				double q_bought = 0.0;

				for (auto ask:scask->second->supply_curve_point_asks){
					//check how much is selling
					//construct a message for an agent
					MessageMarketCCheckAsk* mes = new MessageMarketCCheckAsk (this, price);
					double q_ask = dynamic_cast<ISellerC*>(ask.second->issuer)->_mas_q_sell(mes);
					//clean memory
					delete mes;

					if (q_ask > 0.0){
						if (q_bought < q_bid){

							double q_sell = std::min(q_bid - q_bought, q_ask);

							bool FLAG_CLEARED = _clear(dynamic_cast<ISellerC*>(ask.second->issuer),
													dynamic_cast<IBuyerC*>(bid.second->issuer),
												price,
												q_sell);

							if (FLAG_CLEARED){
								q_bought += q_sell;
								q_ask -= q_sell;
							}else{
                                break;
                            };
						}else{
                            break;
                        };
					};
				};
			};
		};
    };
};

bool
MarketC::_clear(ISellerC* seller, IBuyerC* buyer, double price, double q){
	bool FLAG_CLEARED = true;

	MessageMarketCSellGoodC* mes = new MessageMarketCSellGoodC(this, buyer, price, q);

	FLAG_CLEARED = seller->_sell_asGC(mes);

	delete mes;

	MessagePSSendPayment* mes_ps = new MessagePSSendPayment(seller, price*q, "buy_asGC");

	if (FLAG_CLEARED){
		FLAG_CLEARED = buyer->_PS_send_payment(mes_ps);
	};


	delete mes_ps;

	return FLAG_CLEARED;

};



void
MarketC::AddSellerC(ISellerC* agent) { sellers.push_back(agent); };

void
MarketC::AddBuyerC(IBuyerC* agent) { buyers.push_back(agent); };

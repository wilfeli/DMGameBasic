/*
 * W.h
 *
 *  Created on: Apr 12, 2014
 *      Author: wilfeli
 */

#ifndef W_H_
#define W_H_

#include <vector>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include "Tools.h"

class CB;
class H;
class F;
class MarketL;
class MarketC;
class MessageBankruptcy;
class MessageBankruptcyH;
class ParametersW;



class W{
public:
	W(double);
//    ~W();

	void init();
	int NH = 0;
	int NFGC = 0;
	int NCB = 0;

	int t = -1;

	void AddF(F*);
	void AddH(H*);
	void AddCB(CB*);
	void AddMarketL(MarketL*);
	void AddMarketC(MarketC*);
	void ac_LS_bankruptcy(MessageBankruptcy*);
    void ac_LS_bankruptcy(MessageBankruptcyH*);
	void step();
	void print_step();

    ParametersW* param;
	std::vector<F*> Fs;
	std::vector<H*> Hs;
	MarketL* ml = NULL;
	MarketC* mc = NULL;
	CB* CBs = NULL;

	boost::mt19937 rng;
    Tools::MyRNG myrng;
};




#endif /* W_H_ */

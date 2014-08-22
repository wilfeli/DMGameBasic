/*
 * ISeller.h
 *
 *  Created on: Apr 16, 2014
 *      Author: wilfeli
 */

#ifndef ISELLER_H_
#define ISELLER_H_


#include "PaymentSystem.h"

class Message;
class MessageMarketLCheckAsk;
class MessageMarketLSellC;



class ISeller: virtual public IPSAgent{
public:
	virtual ~ISeller() {};
	virtual double _mas_q_sell(Message*){return 0.0;};
	virtual double get_status() = 0;
};


class ISellerL: public ISeller{
public:
	virtual ~ISellerL() {};

	virtual double _mas_q_sell(MessageMarketLCheckAsk*) = 0;
	virtual void _sell_asCHK(MessageMarketLSellC*) = 0;
	virtual Ask* ac_L() = 0;
};

class ISellerC: public ISeller{
public:
	virtual ~ISellerC() {};
	virtual double _mas_q_sell(MessageMarketCCheckAsk*) = 0;
	virtual bool _sell_asGC(MessageMarketCSellGoodC*) = 0;
	virtual Ask* ac_C() = 0;
};

#endif /* ISELLER_H_ */

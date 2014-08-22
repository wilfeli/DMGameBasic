/*
 * IBuyer.h
 *
 *  Created on: Apr 16, 2014
 *      Author: wilfeli
 */

#ifndef IBUYER_H_
#define IBUYER_H_

#include "PaymentSystem.h"

class ContractHK;
class Message;
class MessageMarketLCheckBid;
class MessageGoodC;


class IBuyer: virtual public IPSAgent{
public:

	virtual ~IBuyer() {};
	virtual double _mas_q_buy(Message*){return 0.0;}
	virtual double get_status() = 0;
};


class IBuyerL: public IBuyer{
public:
	virtual Bid* ac_L() = 0;
	virtual double _mas_q_buy(MessageMarketLCheckBid*) = 0;
	virtual void _buy_asCHK(ContractHK*) = 0;
	virtual ~IBuyerL() {};


};

class IBuyerC: public IBuyer{
public:
	virtual Bid* ac_C() = 0;
	virtual double _mas_q_buy(MessageMarketCCheckBid*) = 0;
	virtual ~IBuyerC() {};
	virtual bool _PS_send_payment(MessagePSSendPayment*) = 0;
	virtual void _buy_asGC(MessageGoodC*) = 0;

};

#endif /* IBUYER_H_ */

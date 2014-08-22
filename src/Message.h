/*
 * Message.h
 *
 *  Created on: Apr 13, 2014
 *      Author: wilfeli
 */

#ifndef MESSAGE_H_
#define MESSAGE_H_

#include <string>

template <typename Ti, typename Th> class Contract;
class MarketL;
class MarketC;
class IBuyerL;
class IBuyerC;
class ISellerL;
class IPSAgent;
class F;
class H;


class Message{
public:
//	Contract<Agent, Agent>* contract = NULL;
	Message() = default;
	Message(std::string);

	std::string type;


};



class MessageMarketLSellC : public Message{
public:
	MessageMarketLSellC(MarketL*, IBuyerL*, double, double);

	MarketL* sender;
	IBuyerL* buyer;
	double p;
	double q;


};

class MessageMarketCSellGoodC : public Message{
public:
	MessageMarketCSellGoodC(MarketC*, IBuyerC*, double, double);

	MarketC* sender;
	IBuyerC* buyer;
	double p;
	double q;


};


class MessagePSSendPayment: public Message{
public:
	MessagePSSendPayment(IPSAgent*, double, std::string);
	MessagePSSendPayment(double, std::string);

	IPSAgent* getter;
	IPSAgent* sender = NULL;
	double q;


};


class MessageMarketLCheckAsk : public Message{
public:

	MessageMarketLCheckAsk(MarketL*, double);

	MarketL* sender;
	double p_eq;

};

class MessageMarketLCheckBid : public Message{
public:

	MessageMarketLCheckBid(MarketL*, double);

	MarketL* sender;
	double p_eq;

};


class MessageMarketCCheckAsk : public Message{
public:

	MessageMarketCCheckAsk(MarketC*, double);

	MarketC* sender;
	double p_eq;

};

class MessageMarketCCheckBid : public Message{
public:

	MessageMarketCCheckBid(MarketC*, double);

	MarketC* sender;
	double p_eq;

};


class MessageGoodC: public Message{
public:
	MessageGoodC(double, double);

	double q;
	double p;
};


class MessageBankruptcy: public Message{
public:
	MessageBankruptcy(F*);

	F* agent;
};

class MessageBankruptcyH: public Message{
public:
	MessageBankruptcyH(H*);
    
	H* agent;
};


class MessageMakeDec: public Message{
public:
	MessageMakeDec();
};


class MessageStatus: public Message{
public:
	MessageStatus(double);
	double status = 1.0;
};


class MessageF_F: public Message{
public:
	double q;
	MessageF_F():q(0.0){type = "F_F";};
};


#endif /* MESSAGE_H_ */

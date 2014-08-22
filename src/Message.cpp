/*
 * Message.cc
 *
 *  Created on: Apr 13, 2014
 *      Author: wilfeli
 */

#include "Message.h"



Message::Message(std::string s_):type(s_){};


MessagePSSendPayment::MessagePSSendPayment(IPSAgent* a_, double q_, std::string type_):
		Message(type_), getter(a_), q(q_){};

MessagePSSendPayment::MessagePSSendPayment(double q_, std::string type_):
		Message(type_), q(q_){};


//MessageSellC::MessageSellC(Message* inf_, Contract<Agent, Agent>* c_): inf(inf_), contract(c_){
//	type = "SellC";
//};


MessageMarketLSellC::MessageMarketLSellC
(MarketL* _sender, IBuyerL* _buyer, double _p, double _q):
sender(_sender), buyer(_buyer), p(_p), q(_q) {
	type = "MarketLSellC";
};


MessageMarketCSellGoodC::MessageMarketCSellGoodC
(MarketC* _sender, IBuyerC* _buyer, double _p, double _q):
sender(_sender), buyer(_buyer), p(_p), q(_q) {
	type = "MarketCSellGoodC";
};



MessageMarketLCheckAsk::MessageMarketLCheckAsk
(MarketL* _sender, double _p_eq): sender(_sender), p_eq(_p_eq){
	type = "MarketLCheckAsk";
};

MessageMarketLCheckBid::MessageMarketLCheckBid
(MarketL* _sender, double _p_eq): sender(_sender), p_eq(_p_eq){
	type = "MarketLCheckBid";
};

MessageMarketCCheckAsk::MessageMarketCCheckAsk
(MarketC* _sender, double _p_eq): sender(_sender), p_eq(_p_eq){
	type = "MarketCCheckAsk";
};

MessageMarketCCheckBid::MessageMarketCCheckBid
(MarketC* _sender, double _p_eq): sender(_sender), p_eq(_p_eq){
	type = "MarketCCheckBid";
};


MessageGoodC::MessageGoodC(double q_, double p_): q(q_), p(p_){
	type = "asGC";
};

MessageBankruptcy::MessageBankruptcy(F* agent_): agent(agent_){
	type = "Bankruptcy";
};
MessageBankruptcyH::MessageBankruptcyH(H* agent_): agent(agent_){
	type = "Bankruptcy";
};

MessageMakeDec::MessageMakeDec(){
	type = "MakeDec";
};


MessageStatus::MessageStatus(double status_): status(status_){
	type = "Status";
};

/*
 * PaymentSystem.h
 *
 *  Created on: Apr 30, 2014
 *      Author: wilfeli
 */

#ifndef PAYMENTSYSTEM_H_
#define PAYMENTSYSTEM_H_

#include "IAgent.h"
#include "Message.h"

class IPSAgent: virtual public IProducerConsumerContractBDt0{
public:

	virtual ~IPSAgent() {};
	virtual ContractBDt0*  _PS_contract() = 0;
	virtual void _PS_receive_payment(MessagePSSendPayment*) = 0;
	virtual bool _PS_accept_PO(MessagePSSendPayment*) = 0;

};



//class IIssuerPS: public IPSAgent{
//public:
//
//	virtual ~IIssuerPS();
//	virtual bool _PS_accept_PO(MessagePSSendPayment*) = 0;
//
//
//};
//
//
//class IHolderPS: public IPSAgent, virtual public Agent{
//
//};








#endif /* PAYMENTSYSTEM_H_ */

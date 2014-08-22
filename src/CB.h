/*
 * CB.h
 *
 *  Created on: Apr 16, 2014
 *      Author: wilfeli
 */

#ifndef CB_H_
#define CB_H_

#include "Agent.h"
#include "PaymentSystem.h"

class ContractBDt0;
class W;
class F;
class H;

class CB : public Agent, public IPSAgent,
			virtual public IProducerConsumerContractBDt0 {
public:
	CB(W&);


	bool _PS_accept_PO(MessagePSSendPayment*);
	void _PS_receive_payment(MessagePSSendPayment*);
	ContractBDt0*  _PS_contract();
	void _create_asCBDt0(F*);
	void _create_asCBDt0(H*);

	ContractBDt0* own_account = NULL;

private:
	W& w;
};



#endif /* CB_H_ */

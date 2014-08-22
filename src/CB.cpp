/*
 * CB.cc
 *
 *  Created on: Apr 16, 2014
 *      Author: wilfeli
 */


#include "CB.h"
#include "W.h"
#include "Message.h"
#include "Contract.h"
#include "Parameters.h"
#include "F.h"
#include "H.h"


CB::CB(W& w_):w(w_){
	w.AddCB(this);
};


bool
CB::_PS_accept_PO(MessagePSSendPayment* mes){
	//subtract from sender , add to getter money

	ContractBDt0* c_to = mes->getter->_PS_contract();
	ContractBDt0* c_from = mes->sender->_PS_contract();

	bool FLAG_CLEARED = false;

	if (c_from->q >= mes->q){
		c_from->q -= mes->q;
		c_to->q += mes->q;

		FLAG_CLEARED = true;
	};

	mes->getter->_PS_receive_payment(mes);

	return FLAG_CLEARED;
};

void
CB::_PS_receive_payment(MessagePSSendPayment* mes){
};


ContractBDt0*
CB::_PS_contract(){
	return own_account;
};

void
CB::_create_asCBDt0(F* a){
	//new contract
	//money from w parameter
	ContractBDt0* c = new ContractBDt0(this, a, w.param->F_asCBDt0_Q);
	ass_asCBDt0.push_back(c);

	a->ass_asCBDt0.push_back(c);

};


void
CB::_create_asCBDt0(H* a){
	//new contract
	//money from w parameter
	ContractBDt0* c = new ContractBDt0(this, a, w.param->H_asCBDt0_Q);
	ass_asCBDt0.push_back(c);

	a->ass_asCBDt0.push_back(c);

};

/*
 * IAgent.cc
 *
 *  Created on: Apr 30, 2014
 *      Author: wilfeli
 */

#include "IAgent.h"
#include "Contract.h"
#include "Good.h"

void
IProducerConsumerContractHK::get_asCHK(ContractHK* c){
	ass_asCHK.push_back(c);
};



void
IProducerConsumerContractBDt0::get_asCBDt0(ContractBDt0* c){
	ass_asCBDt0.push_back(c);
};


void
IProducerConsumerGoodC::get_asGoodC(GoodC* good){
	ass_asGC.push_back(good);
};




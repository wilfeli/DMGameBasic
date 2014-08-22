/*
 * IAgent.h
 *
 *  Created on: Apr 30, 2014
 *      Author: wilfeli
 */

#ifndef IAGENT_H_
#define IAGENT_H_

#include <vector>

class ContractHK;
class ContractBDt0;
class GoodC;


class IProducerConsumerContractHK{
public:
	void get_asCHK(ContractHK*);

	std::vector<ContractHK*> ass_asCHK; //holds labor contracts


};

class
IProducerConsumerContractBDt0{
public:
	void get_asCBDt0(ContractBDt0*);

	std::vector<ContractBDt0*> ass_asCBDt0; //holds deposit contracts
};



class IProducerConsumerGoodC{
public:
	void get_asGoodC(GoodC*);

	std::vector<GoodC*> ass_asGC; //holds goods list

};


#endif /* IAGENT_H_ */

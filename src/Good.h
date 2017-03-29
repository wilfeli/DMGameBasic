/*
 * Good.h
 *
 *  Created on: Apr 30, 2014
 *      Author: wilfeli
 */

#ifndef GOOD_H_
#define GOOD_H_

#include "IAgent.h"
#include "Asset.h"




class GoodC: public Asset<IProducerConsumerGoodC>{
public:
	GoodC(IProducerConsumerGoodC*, double);


};


class HK{
public:
	HK(double);
	double q;
};


#endif /* GOOD_H_ */

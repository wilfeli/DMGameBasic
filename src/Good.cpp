/*
 * Good.cpp
 *
 *  Created on: Apr 30, 2014
 *      Author: wilfeli
 */


#include "Good.h"
#include "Asset.h"

GoodC::GoodC(IProducerConsumerGoodC* a_, double q_): Asset<IProducerConsumerGoodC>(std::string("asGC"), a_, q_){};

HK::HK(double q_):q(q_){};

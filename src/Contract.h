/*
 * contract.h
 *
 *  Created on: Apr 12, 2014
 *      Author: wilfeli
 */

#ifndef CONTRACT_H_
#define CONTRACT_H_

#include "Asset.h"
#include "PaymentSystem.h"

class IBuyerL;
class ISellerL;
class IPSAgent;




//template<typename Ti, typename Th> class Contract: public Asset<Th>{
////	auto holder;
////	auto issuer;
//public:
//
//	Contract(Ti*, Th*);
//	Ti* issuer;
//
//
//};
//
//
//class ContractHK: public Contract<ISellerL, IBuyerL>{
//public:
//	ContractHK(ISellerL*, IBuyerL*,
//			double, double,
//			int, int);
//
//	double p;
//	int t_begin;
//	int t_end;
////	ISellerL* issuer;
////	IBuyerL* holder;
//
////	double q;
////	std::string type = "asCHK";
//
//
//
//};
//
//
//class ContractBDt0: public Contract<IPSAgent, IPSAgent>{
//public:
//	ContractBDt0(IPSAgent*, IPSAgent*, double);
//};

class ContractHK{
public:
	ContractHK(ISellerL*, IBuyerL*,
			double, double,
			int, int);

	double p;
	int t_begin;
	int t_end;
	ISellerL* issuer;
	IBuyerL* holder;

	double q;
	std::string type;



};


class ContractBDt0{
public:
	ContractBDt0(IPSAgent*, IPSAgent*, double);
	IPSAgent* issuer;
	IPSAgent* holder;

	double q;
	std::string type;

};



#endif /* CONTRACT_H_ */

/*
 * contract.cc
 *
 *  Created on: Apr 12, 2014
 *      Author: wilfeli
 */

#include "Contract.h"

//has holder
//issuer
//payment schedule? or parameters
//timing
//template <typename Ti, typename Th>
//Contract<Ti, Th>::Contract(Ti* issuer_, Th* holder_):Asset<Th>(holder_), issuer(issuer_){};
//
//
//
//ContractHK::ContractHK(ISellerL* issuer_, IBuyerL* holder_,
//						double p_, double q_,
//						int t_begin_, int t_end_):Contract(issuer, holder),
//													p(p_),
//													t_begin(t_begin_),
//													t_end(t_end_){
//	type = "asCHK";
//	q = q_;
//};
//
//ContractBDt0::ContractBDt0(IPSAgent* issuer_, IPSAgent* holder_, double q_):Contract(issuer, holder){
//	type = "asCBDt0";
//	q = q_;
//};

ContractHK::ContractHK(ISellerL* issuer_, IBuyerL* holder_,
						double p_, double q_,
						int t_begin_, int t_end_):issuer(issuer_), holder(holder_),
													p(p_),
													t_begin(t_begin_),
													t_end(t_end_){
	type = "asCHK";
	q = q_;
};

ContractBDt0::ContractBDt0(IPSAgent* issuer_, IPSAgent* holder_, double q_):issuer(issuer_), holder(holder_){
	type = "asCBDt0";
	q = q_;
};

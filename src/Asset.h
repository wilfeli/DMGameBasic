/*
 * asset.h
 *
 *  Created on: Apr 12, 2014
 *      Author: wilfeli
 */

#ifndef ASSET_H_
#define ASSET_H_

#include <string>


template<typename Th>class Asset{
public:

//	Asset() = default;
	Asset(Th*);
	Asset(std::string s_, Th *a_, double q_):
		           type(s_), holder(a_), q(q_) { };

	//by default holder and issuer? types

	std::string type;
	Th* holder;
	double q;



};



#endif /* ASSET_H_ */

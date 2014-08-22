/*
 * asset.cc
 *
 *  Created on: Apr 12, 2014
 *      Author: wilfeli
 */


#include "Asset.h"

template<typename Th>
Asset<Th>::Asset(Th* holder_):holder(holder_){
	q = 0.0;
};

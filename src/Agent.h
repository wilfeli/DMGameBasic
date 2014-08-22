/*
 * Agent.h
 *
 *  Created on: Apr 12, 2014
 *      Author: wilfeli
 */

#ifndef AGENT_H_
#define AGENT_H_

#include <string>

class Agent{
public:
	int status = 1;
	int t_begin_life = 0;
	double t_end_life = 1.0/0.0;
	std::string type = "Agent";

};



#endif /* AGENT_H_ */

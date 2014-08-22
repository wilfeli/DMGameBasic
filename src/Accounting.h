/*
 * Accounting.h
 *
 *  Created on: Apr 24, 2014
 *      Author: wilfeli
 */

#ifndef ACCOUNTING_H_
#define ACCOUNTING_H_

#include <list>
#include <vector>
#include <Eigen/Dense>

class Message;
class ContractHK;
class MessageMarketLSellC;
class MessagePSSendPayment;
class MessageMarketCSellGoodC;
class MessageStatus;
class MessageGoodC;
class MessageF_F;



class AccountF{
public:
	double cost_p = 0.0;
	double cost_q = 0.0;
    std::string ACCOUNTING_TYPE  = "classic";

	std::vector<double> expences_wage;
	std::vector<double> asCHK_q;
	std::vector<double> sales_q;
	std::vector<double> sales_p;
	std::vector<double> sales_pq;
	std::vector<double> production_q;
	std::vector<double> profit;

	void ac_get_inf(Message*);
	void ac_get_inf(Message*, ContractHK*);
	void ac_get_inf(MessageF_F*);
	void ac_get_inf(MessageMarketCSellGoodC*);

	void ac_W_initialize_step();
	void ac_W_initialize_step(MessageStatus*);
	void ac_W_end_step(Eigen::MatrixXd*);
};


class AccountH{
public:
	std::vector<double> asCHK_p;
	std::vector<double> asCHK_q;
	std::vector<double> asGC_p;
	std::vector<double> asGC_pq;
	std::vector<double> asGC_q;
	std::vector<double> g_t; //realized goal
	std::vector<double> g_GC_q_t; //realized consumed good
	std::vector<double> g_CHK_q_t; //realized supplied labor
	std::vector<double> FI_div_t; //received dividends
	std::vector<double> asCHK_p_t; //realized wage payment
	std::vector<double> asCBDt0_pq_t; //realized money at account





//	void ac_get_inf(MessageSellC);
	void ac_get_inf(MessageMarketLSellC*, ContractHK*);
	void ac_get_inf(MessagePSSendPayment*);
	void ac_get_inf(MessageGoodC*);

	void ac_W_initialize_step();
    void ac_W_initialize_step(MessageStatus*);
	void ac_W_end_step(Eigen::MatrixXd*, double);




};


#endif /* ACCOUNTING_H_ */

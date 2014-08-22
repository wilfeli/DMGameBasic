/*
  * F.h
 *
 *  Created on: Apr 13, 2014
 *      Author: wilfeli
 */

#ifndef F_H_
#define F_H_

#include <memory>
#include <Eigen/Dense>
#include <map>
#include <list>


//#include "W.h" //or simply class W; as forward declaration
#include "Agent.h"
#include "IAgent.h"
#include "Message.h"
#include "MarketTools.h"
#include "ISeller.h"
#include "IBuyer.h"
#include "Asset.h"
#include "Expectations.h"
#include "Accounting.h"
#include "Parameters.h"
#include "Contract.h"
#include "Good.h"

class W;




class F : virtual public Agent, public IBuyerL, public ISellerC,
virtual public IProducerConsumerContractBDt0,
public IProducerConsumerContractHK,
public IProducerConsumerGoodC{
public:
	F(W&, int);
    void init();
//	~F() = default;




	Bid* ac_L();
	Ask* ac_C();
	Eigen::MatrixXd _s0();

	double get_status();
	void opt_CS(Eigen::MatrixXd*, Eigen::MatrixXd*, int N = 10);
    void create_decision_grid(Eigen::MatrixXd&, int);
	void _step_opt_CS(Eigen::MatrixXd*, Eigen::MatrixXd*, Eigen::MatrixXd*, Eigen::MatrixXd*, Eigen::MatrixXd*, int);
    void _step_opt_CS(Eigen::MatrixXd&, Eigen::MatrixXd*, Eigen::MatrixXd*, Eigen::MatrixXd*, Eigen::MatrixXd*, Eigen::MatrixXd*, Eigen::MatrixXd&, Eigen::MatrixXd*, int);
    
    void _opt_CS_s_V(Eigen::Block<Eigen::MatrixXd>&,  Eigen::MatrixXd&, Eigen::MatrixXd&, int);
    double _c_theta_x(Eigen::MatrixXd&, Eigen::Block<Eigen::MatrixXd>&, Eigen::MatrixXd&);
    
    
    Eigen::MatrixXd _API_LM(Eigen::MatrixXd*, int N=5, int M=5); //optimal choice given ADP algorithm
    void _bf(Eigen::MatrixXd*, Eigen::MatrixXd*); //basis functions
    void _bf(Eigen::Block<Eigen::MatrixXd>&, Eigen::MatrixXd*); //basis functions if passed block in
    double _V(Eigen::MatrixXd*, Eigen::MatrixXd*, Eigen::MatrixXd&); //value function estimation
    double _V(Eigen::Block<Eigen::MatrixXd>&, Eigen::MatrixXd*, Eigen::MatrixXd&);
    
    //    template<class T>
    //    double _V(T s, Eigen::MatrixXd*, Eigen::MatrixXd&);
    Eigen::MatrixXd _RLS(Eigen::MatrixXd&, Eigen::MatrixXd&, Eigen::Block<Eigen::MatrixXd>&, int i, double c);
    
    //QL part
    void opt_QL(); //opt decisions using QL algorithm
    void _set_QL();//calls QL setup
    std::vector<double> _place_state_QL(Eigen::MatrixXd&); //place state in bin
    bool _is_bin(int, int, std::vector<double>&); //checks that part of the state is in this bin
    void wm_update_Q(); //update Q in QL
    void create_bin_values(Eigen::MatrixXd&, std::vector<std::vector<double>>&); //creates bins with values for grid
    
    //mRE part
    void opt_mRE();
    void _set_mRE(); //setup mRE matrix
    void wm_update_mRE(); //update mRE matrix in mRE
    
	double _mas_q_buy(MessageMarketLCheckBid*);
	double _mas_q_sell(MessageMarketCCheckAsk*);
	void _buy_asCHK(ContractHK*);
	bool _sell_asGC(MessageMarketCSellGoodC*);

	void _PS_receive_payment(MessagePSSendPayment*);
	bool _PS_accept_PO(MessagePSSendPayment*);

	void ac_W(MessageMakeDec*); //make decisions
	ContractBDt0* _PS_contract();

	void ac_W_begin_step(); //begin step
	void ac_W_begin_step(MessageStatus*); //begin step if bankrupt
	void ac_W_end_step(); //end step
	void ac_W(MessageF_F*); //call to produce

	double _F_F(Eigen::MatrixXd*); //production

	void wm_update_k_s(); //update wm
	Eigen::MatrixXd wm_w0_tbeg(); //form w0 for wm
	void wm_update_expectation_backward(ExpectationBackward*, int);
	Eigen::MatrixXd wm_mean_variance(std::vector<double> &);

	W& w;
	Bid l_bid;
	Ask c_ask;
	int id;

	std::map<std::string, ExpectationBackward*> wm;
	AccountF* account;
	Parameters* param;
    
    Eigen::MatrixXd theta_x_all; //all decision grid parameters
    Eigen::MatrixXd theta_V; //ADP algorithm parameters for V
    Eigen::MatrixXd B_n_1_RLS; //stores ADP _RLS matrix

    Eigen::MatrixXd Q; //Q matrix
    Eigen::MatrixXd Q_grid;//grid for Q matrix - ranges for bins and step
    std::vector<std::vector<double>> Q_bin; //values for bins
    std::vector<double> s_t_bin; //stores bin at decision time
    std::vector<int> theta_x_t_bin;//stores decision bin

    
    //mRE part
    Eigen::MatrixXd mRE; //mRE matrix
    Eigen::MatrixXd mRE_grid;//grid for mRE matrix - ranges for bins and step
    std::vector<std::vector<double>> mRE_bin; //values for bins
    
    
	Eigen::MatrixXd s0; //state storage
	Eigen::MatrixXd theta_x; //decision storage
    

};





#endif /* F_H_ */

//
//  H.h
//  minimacro
//
//  Created by Ekaterina Sinitskaya on 6/3/14.
//  Copyright (c) 2014 Ekaterina Sinitskaya. All rights reserved.
//

#ifndef H_H_
#define H_H_

#include <map>
#include "Agent.h"
#include "IAgent.h"
#include "MarketTools.h"
#include "IBuyer.h"
#include "ISeller.h"
#include "Message.h"
#include "Asset.h"
#include "Expectations.h"
#include "Accounting.h"
#include "Parameters.h"
#include "Contract.h"
#include "Good.h"

class W;




class H : virtual public Agent, public IBuyerC, public ISellerL,
virtual public IProducerConsumerContractBDt0,
public IProducerConsumerContractHK,
public IProducerConsumerGoodC{
    //			virtual public IPSAgent{
public:
	H(W&, int);
    void init();
    
	AccountH* account;
	std::map<std::string, ExpectationBackward*> wm;
	Parameters* param;
    
	double get_status();
	double _mas_q_buy(MessageMarketCCheckBid*);
	double _mas_q_sell(MessageMarketLCheckAsk*);
	void _sell_asCHK(MessageMarketLSellC*);
	void _buy_asGC(MessageGoodC*);
	void ac_W(MessageMakeDec*);
	Ask* ac_L();
	Bid* ac_C();
    
	bool _PS_send_payment(MessagePSSendPayment*);
	void _PS_receive_payment(MessagePSSendPayment*);
	bool _PS_accept_PO(MessagePSSendPayment*);
	ContractBDt0*  _PS_contract();
    
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

    
    
    
    
	Eigen::MatrixXd wm_w0_tbeg();
	void wm_update_k_s();
	void wm_update_expectation_backward(ExpectationBackward*, int, bool);
	Eigen::MatrixXd wm_mean_variance(std::vector<double> &);
    
	void ac_W_end_step();
	void ac_W_begin_step();
    void ac_W_begin_step(MessageStatus*); //begin step if dead
    void ac_W_end_step(MessageStatus*); //end step if dead
    
	double _goal_t(Eigen::MatrixXd*);
    
	Eigen::MatrixXd _s0();
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

    
    
    
    
    
	void get_asHK(HK*);
    
	std::vector<HK*> ass_asHK; //holds goods list
    
private:
	W& w;
	Ask l_ask;
	Bid c_bid;
    
	int id;
	Eigen::MatrixXd theta_x; //decision storage
	Eigen::MatrixXd s0;
    
    
    
    
};





#endif

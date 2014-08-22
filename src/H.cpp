/*
 * H.cpp
 *
 *  Created on: Apr 11, 2014
 *      Author: wilfeli
 *
 *
 *
 *
 *
 *
 *
 *
 *
 */


#include <Eigen/Dense>
#include <math.h>

#include "H.h"
#include "Asset.h"
#include "W.h"
#include "Contract.h"
#include "Market.h"
#include "Tools.h"


using Eigen::MatrixXd;



H::H(W &w_, int id_):w(w_),
l_ask(0.0, 0.0),
c_bid(0.0, 0.0),
id(id_),
theta_x(1,3),
s0(1,18){
	// create reference to the world
	w.AddH(this);
    
	//create storage for goods
	ass_asGC.push_back(new GoodC(this, 0.0));
    
	//create HK
	ass_asHK.push_back(new HK(1.0));
    
    
	type = "H";
    
	c_bid.issuer = this;
	l_ask.issuer = this;
    
	//create wm
	ExpectationBackward* exp = new ExpectationBackward(1.0, 10, 0.5);
	wm["MasCHK"] = exp;
    
	//create wm
	exp = new ExpectationBackward(1.0, 10, 0.5);
	wm["MasGC"] = exp;
    
	//create wm
	exp = new ExpectationBackward(1.0, 0, 0.0);
	wm["api"] = exp;
    
    
	//create wm
	exp = new ExpectationBackward(0.0, 0, 0.01);
	wm["fi_ti"] = exp;
    
	//parameters
	param = new Parameters();
	param->init();
    
	//accounting
	account = new AccountH();
    
    
};

void
H::init(){
    //ADP setup
    theta_V = Eigen::MatrixXd::Constant(1, 1, param->ADP_BF_N);
    B_n_1_RLS = Eigen::MatrixXd::Zero(param->ADP_BF_N, param->ADP_BF_N);
    
    //QL setup
    _set_QL();
    
    
    //mRe setup
    _set_mRE();

    
};

void
H::ac_W(MessageMakeDec* mes){
    
	//calls opt choice
	_s0();
	//calls opt choice
    if (param->opt_TYPE == "EO_CS"){
        opt_CS(&s0,&theta_x, param->opt_CS_N);
    };
    
    if (param->opt_TYPE == "EO_ADP"){
        theta_V = _API_LM(&s0, param->ADP_API_LM_N, param->ADP_API_LM_M);

        Eigen::Block<Eigen::MatrixXd> s = s0.block(0,0,s0.rows(),s0.cols());
        _opt_CS_s_V(s,theta_x, theta_V, param->opt_CS_N);
    };
    
    if (param->opt_TYPE == "QL"){
        opt_QL();
    };
    
    if (param->opt_TYPE == "mRE"){
        opt_mRE();
        
    };
    
    
	//transforms it into ask and bid
	l_ask.q = theta_x(0,0) * s0(0,1);
	l_ask.p = theta_x(0,1) * s0(0,4);
    
	c_bid.pq_share = theta_x(0,2);
    
    
    
};


Ask*
H::ac_L(){
	//update numbers in bid
	l_ask.supply_curve = Eigen::MatrixXd::Zero(1,2);
	l_ask.supply_curve << l_ask.p , l_ask.q;
    
	l_ask.p_t_0 = l_ask.p;
	l_ask.q_t_0 = l_ask.q;
    
    
    
	//return reference to the ask
	return &l_ask;
};

Bid*
H::ac_C(){
	//update numbers in bid
	c_bid.pq = _s0()(0,0) * c_bid.pq_share;
    
	c_bid.pq_t_0 = c_bid.pq;
	c_bid.pq_share_t_0 = c_bid.pq_share;
    
    
	//update numbers in bid
    
    
//    if (param->GRID == "test"){
//        c_bid.demand_curve = Eigen::MatrixXd::Zero(1,2);
//        c_bid.demand_curve.row(0) << wm["MasGC"]->mu , c_bid.pq/wm["MasGC"]->mu ;
//        return &c_bid;
//    };
    
    
	//number of steps
	int N_STEPS_dc = 30;
	double P_MAX_dc = wm["MasGC"]->mu * 10;
	double P_MIN_dc = wm["MasGC"]->mu * 0.01;
	double STEP_SIZE = (P_MAX_dc - P_MIN_dc)/N_STEPS_dc;
    
	c_bid.demand_curve = Eigen::MatrixXd::Zero(N_STEPS_dc,2);
	for (int i = 0; i < N_STEPS_dc; i++){
        
		c_bid.demand_curve.row(i) << P_MIN_dc + i * STEP_SIZE, c_bid.pq/(P_MIN_dc + i * STEP_SIZE);
        
	};
    
    
    
    
	//return reference to the bid
	return &c_bid;
};




double
H::_mas_q_sell(MessageMarketLCheckAsk* inf){
	//checks how much to sell
	//return q from ask
    
	double q_sell;
    
	if (inf->p_eq >= l_ask.p){
		q_sell = l_ask.q;
	};
    
	return q_sell;
};

double
H::_mas_q_buy(MessageMarketCCheckBid* inf){
	//checks how much to sell
	//return q from ask
    
    
	double q_buy;
    
	int i=0;
    
	while ((c_bid.demand_curve(i,0) < inf->p_eq) && (i<c_bid.demand_curve.rows())){
		i++;
	};
    
	if (i<c_bid.demand_curve.rows()){
		q_buy = c_bid.demand_curve(i,1);
	}else{
		//buys zero
		q_buy = 0.0;
	};
    
    
    
	return q_buy;
};


void
H::_sell_asCHK(MessageMarketLSellC* inf){
	//new contract
	ContractHK* cHK = new ContractHK(this,
                                     inf->buyer,
                                     inf->p,
                                     inf->q,
                                     w.t,
                                     w.t);
    
	type = "H_EO";
    
	ass_asCHK.push_back(cHK);
    
	//?update ask
	l_ask.q -= cHK->q;
    
	//tell market about the contract
	inf->sender->_clear_C(cHK);
    
	//call accounting
	//create accounting message
    //	MessageSellC* mes = new MessageSellC(inf, cHK);
	account->ac_get_inf(inf, cHK);
    
    //	delete mes;
    
    
};

double
H::_goal_t(Eigen::MatrixXd* s){
    
	double l;
	double c;
    
	if ((*s)(0,3) <= 0.0){
		c = -0.5;
//        c = 0.0;
	}else{
        c = (*s)(0,3);
    };
    
	l = (*s)(0,2);
    
	return (param->GOAL_T_theta(0,0) * log(1+c) + param->GOAL_T_theta(0,1) * (1-l));
//    return (param->GOAL_T_theta(0,0) * log(0.5+c) + param->GOAL_T_theta(0,1) * (1-l));
};



void
H::_buy_asGC(MessageGoodC* mes){
    
	ass_asGC.back()->q += mes->q;
	account->ac_get_inf(mes);
    
};


void
H::_PS_receive_payment(MessagePSSendPayment* mes){
    
	account->ac_get_inf(mes);
    
};


bool
H::_PS_send_payment(MessagePSSendPayment* inf){
    
	inf->sender = this; //static_cast<IHolderPS*>(const_cast<H*>(this));
    
	bool FLAG_CLEARED = ass_asCBDt0.front()->issuer->_PS_accept_PO(inf);
    
    
    
	if(FLAG_CLEARED){
		account->ac_get_inf(inf);
        
	};
    
	return FLAG_CLEARED;
    
    
};

void
H::ac_W_begin_step(){
    
	//update accounting
	account->ac_W_initialize_step();
    
	//
	for (auto c:ass_asCHK){
		if (c->t_end <= w.t){
			delete c;
			c = NULL;
		};
	};
    
    
	//clear contracts
	ass_asCHK.erase(std::remove_if(ass_asCHK.begin(), ass_asCHK.end(),
                                   [&](ContractHK* x) -> bool { return (x); }),
                    ass_asCHK.end());
    
    
    
};

void
H::ac_W_begin_step(MessageStatus* mes){
	if (mes->status <= 0.0){
		account->ac_W_initialize_step(mes);
	};
};

void
H::ac_W_end_step(){
	//updates state
	_s0();
	double ut = _goal_t(&s0);
	account->ac_W_end_step(&s0, ut);
    
	//eats good
	ass_asGC.front()->q = 0.0;
    
};

void
H::ac_W_end_step(MessageStatus* mes){
    if (mes->status <= 0.0){
        ac_W_end_step();
    };
};



void
H::wm_update_k_s(){
	if (w.t > 0.0){
        
		MatrixXd w0 = wm_w0_tbeg();
		ExpectationBackward* exp_i;
		std::vector<double> s_t_i;
        
		//update k_s if conditions are met
//		if ((w0(0,1) > 0.0) || ((w0(0,1) == 0.0) && (w0(0,0) != 0.0)) || ((w0(0,4) > 0.0) && (w0(0,0) == 0.0))){
			s_t_i.push_back(w0(0,0));
			s_t_i.push_back((w0(0,1) > 0.0)? w0(0,1):w0(0,7));
            
			exp_i = wm["MasCHK"];
            
			exp_i->s_t.push_back(s_t_i);
            
			wm_update_expectation_backward(exp_i, 1, true);
//		};
        
		s_t_i.clear();
        
		s_t_i.push_back(w0(0,2));
		s_t_i.push_back((w0(0,2) > 0.0)? w0(0,3)/w0(0,2):w0(0,8));
        
        
		exp_i = wm["MasGC"];
        
		exp_i->s_t.push_back(s_t_i);
        
		wm_update_expectation_backward(exp_i, 1, true);
        
		s_t_i.clear();
        
        
		if (w.t > 0.0){
			exp_i = wm["api"];
            
			s_t_i.push_back(w0(0,5));
			exp_i->s_t.push_back(s_t_i);
            
			wm_update_expectation_backward(exp_i, 0, false);
            
            
			s_t_i.clear();
            
			exp_i = wm["fi_ti"];
            
			s_t_i.push_back(w0(0,6));
			exp_i->s_t.push_back(s_t_i);
            
//			wm_update_expectation_backward(exp_i, 0, false);
            wm_update_expectation_backward(exp_i, 0, true);
            
		};
        
        
        if (param->opt_TYPE == "QL"){
            wm_update_Q();
            
        };
        
        if (param->opt_TYPE == "mRE"){
            wm_update_mRE();
            
        };
    
        
	};
    
    
    
};


void
H::wm_update_expectation_backward(ExpectationBackward* exp_i, int p_position, bool FLAG_UPDATE_V){
    
	std::vector<double> p;
    
	int i_s_t_1_max = 0;
	long i_s_t_max = 0;
    
	//depending on the length of the wm and accumulated number of prices
	if (std::isinf(param->WM_LENGTH)){
		i_s_t_1_max = exp_i->s_t_1["n"];
		i_s_t_max = exp_i->s_t.size();
	}else{
		i_s_t_1_max = std::min(std::max(param->WM_LENGTH - exp_i->s_t.size(), 0.0), exp_i->s_t_1["n"]);
		i_s_t_max = std::min(param->WM_LENGTH, (double)exp_i->s_t.size());
	};
    
	for (int i = 0; i < i_s_t_1_max; i++){
		p.push_back(exp_i->s_t_1["mu"]);
	};
    
	//push other prices
	for (std::size_t i = (exp_i->s_t.size() - i_s_t_max);i < exp_i->s_t.size();i++){
		p.push_back(exp_i->s_t[i][p_position]);
	};
    
	MatrixXd mean_variance = wm_mean_variance(p);
    
	exp_i->mu = mean_variance(0,0);
	exp_i->n += 1.0;
    
    if (FLAG_UPDATE_V){
        exp_i->v = mean_variance(1,0);
        
        if (exp_i->v <= 0.0){
            exp_i->v = exp_i->mu * 0.01;
        };
        
        
    };
    
    
    
};



Eigen::MatrixXd
H::wm_mean_variance(std::vector<double> &x){
	MatrixXd mean_variance(2,1);
    
	//gets mean and variance
	double sum = std::accumulate(x.begin(), x.end(), 0.0);
	double mean = sum / x.size();
	std::vector<double> diff(x.size());
	std::transform(x.begin(), x.end(), diff.begin(),
	               std::bind2nd(std::minus<double>(), mean));
	double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    //	double stdev = std::sqrt(sq_sum / x.size());
	double variance = sq_sum/x.size();
    
	mean_variance(0,0) = mean;
	mean_variance(1,0) = variance;
    
	return mean_variance;
};


Eigen::MatrixXd
H::wm_w0_tbeg(){
	MatrixXd w0(1,9);
    
	long life_length = account->asCHK_q.size();
    
	w0(0,0) = account->asCHK_q.at(life_length - 2);
	w0(0,1) = account->asCHK_p.at(life_length - 2);
	w0(0,2) = account->asGC_q.at(life_length - 2);
	w0(0,3) = account->asGC_pq.at(life_length - 2);
    
	w0(0,4) = l_ask.q_t_0;
	w0(0,5) = account->FI_div_t.at(life_length - 2) + account->asCHK_p_t.at(life_length - 2);
	w0(0,6) = account->FI_div_t.at(life_length - 2);
    
	w0(0,7) = w.ml->market_price.at(w.t - 1)->p;
	w0(0,8) = w.mc->market_price.at(w.t - 1)->p;
    
    //	std::cout << w0 << std::endl;
    
    
    
	return w0;
};


Eigen::MatrixXd
H::_s0(){
	//money holdings
	s0(0,0) = ass_asCBDt0.front()->q;
    
	s0(0,1) = ass_asHK.front()->q;
    
	s0(0,2) = 0.0;
	//labor to be supplied
	for (auto c:ass_asCHK){
		s0(0,2) += c->q;
	};
    
    
	//amount of good c
	s0(0,3) = ass_asGC.front()->q;
    
	//expected price on labor market
	s0(0,4) = wm["MasCHK"]->mu;
    
	//number of observations
	s0(0,5) = wm["MasCHK"]->n;
    
	//expected price on goods market
	s0(0,6) = wm["MasGC"]->mu;
    
	//number of observations
	s0(0,7) = wm["MasGC"]->n;
    
    
    
    
	//average past income
	s0(0,8) = wm["api"]->mu;
    
	s0(0,9) = wm["api"]->n;
    
	//current income
	s0(0,10) = account->FI_div_t.back() + account->asCHK_p_t.back();
    
	//to be used in max routinue
	//current signed labor contracts - q
	s0(0,11) = 0.0;
    
	//current signed labor contracts - p
	s0(0,12) = 0.0;
    
	//financial income
	s0(0,13) = wm["fi_ti"]->mu;
    
	s0(0,14) = wm["fi_ti"]->n;
    
    
	//variance of labor market
	s0(0,15) = wm["MasCHK"]->v;
    
	//variance of goods market
	s0(0,16) = wm["MasGC"]->v;
    
	//variance of financial income
	s0(0,17) = wm["fi_ti"]->v;
    
    
	//
    
    
	return s0;
};




void
H::opt_CS(MatrixXd* s0,  MatrixXd* theta_x, int N){
    
//	N = 1;
    
	//gets matrix of thetas
	//prepare initial matrix
	MatrixXd m_s0 = MatrixXd::Zero(N, s0->cols());
    
	m_s0.rowwise() += s0->row(0);
    
    
	//allocate space to the matrix
	MatrixXd n_w;
	MatrixXd n_x;
	MatrixXd n_s = MatrixXd::Zero(N, s0->cols());
	MatrixXd n_c;
    
    
	//create matrix for the grid - small size here
	MatrixXd grid(3,3);
    
    create_decision_grid(grid,0);

    
    
    if (w.param->SIMULATION_MODE == "test"){
        Tools::mgrid_test(grid, &theta_x_all);
        
    } else{
        
        Tools::mgrid(grid, &theta_x_all);
    };
    
    //	std::cout << theta_x_all.transpose() << "\n";
    
    
	//create random number generators for the implementation
    
	boost::normal_distribution<> nd_w1((*s0)(0,4), pow((*s0)(0,15),0.5));
	boost::normal_distribution<> nd_w2((*s0)(0,6), pow((*s0)(0,16),0.5));
	boost::normal_distribution<> nd_w3((*s0)(0,13), pow((*s0)(0,17),0.5));
    
	boost::variate_generator<boost::mt19937&,
    boost::normal_distribution<>> rng_w1(w.rng, nd_w1);
    
	boost::variate_generator<boost::mt19937&,
    boost::normal_distribution<>> rng_w2(w.rng, nd_w2);
    
	boost::variate_generator<boost::mt19937&,
    boost::normal_distribution<>> rng_w3(w.rng, nd_w3);
    

    
    if (w.param->SIMULATION_MODE != "test"){

        //draw random variables
        //fast realization
        n_w = MatrixXd::Zero(N, 3 * param->T_MAX);
        for (int i = 0; i < N; i++){
            for(int j=0; j< n_w.cols(); j++){
                if (j%3 == 0){
                    n_w(i,j) = rng_w1();
                } else {
                    if (j%3 == 1){
                        n_w(i,j) = rng_w2();
                    }else{
                        n_w(i,j) = rng_w3();
                    };
                    
                };
                n_w(i,j) = std::max(0.0, (double)n_w(i,j));
            };
        };
    };
    

    
    
    
    
    
	//matrix of all
	MatrixXd c_all(theta_x_all.cols(),1);
    
	for (int i = 0 ; i < theta_x_all.cols(); i++){
        if (w.param->SIMULATION_MODE == "test"){
            n_w = MatrixXd::Zero(N, 3 * param->T_MAX);
            for (int i = 0; i < N; i++){
                for(int j=0; j< n_w.cols(); j++){
                    if (j%3 == 0){
                        n_w(i,j) = Tools::get_normal((*s0)(0,4), pow((*s0)(0,15),0.5), w.myrng);
                    } else {
                        if (j%3 == 1){
                            n_w(i,j) = Tools::get_normal((*s0)(0,6), pow((*s0)(0,16),0.5), w.myrng);
                        }else{
                            n_w(i,j) = Tools::get_normal((*s0)(0,13), pow((*s0)(0,17),0.5), w.myrng);
                        };
                        
                    };
                    n_w(i,j) = std::max(0.0, (double)n_w(i,j));
                };
            };
        };

		//call estimator of results of steps, given theta_x
        theta_x->row(0) = theta_x_all.col(i).transpose();
		n_x = MatrixXd::Zero(N, 3);
		n_s = m_s0;
		n_c = MatrixXd::Zero(N, param->T_MAX);
		_step_opt_CS(theta_x, &n_w, &n_x, &n_s, &n_c, param->T_MAX);
        
		//after n_c
		//average over n runs
//        std::cout << n_w << std::endl;
//        std::cout << n_c << std::endl;
		c_all(i,0) = (n_c * param->BETA_T).mean();
        
	};

    
    
    
    //	std::cout << c_all << "\n";
    
	//find all max elements
	double max = c_all(0,0);
	std::list<int> i_max;
	for (int i=0; i<c_all.rows(); i++){
		if (c_all(i,0) == max){
			i_max.push_back(i);
		}else{
			if (c_all(i,0) > max){
				max = c_all(i,0);
				i_max.clear();
				i_max.push_back(i);
			};
		}
	};
    
	//if more than 1 max - randomly pick
	long theta_max_i;
	if (i_max.size()>1){
		boost::random::uniform_int_distribution<> i_max_dist(0, i_max.size()-1);
		theta_max_i = i_max_dist(w.rng);
	}else{
		theta_max_i = i_max.front();
	};
    
    
	theta_x->row(0) = theta_x_all.col(theta_max_i).transpose();
    
    //	std::cout << "theta" << *theta_x << "\n";
};


double
H::get_status(){
	return status;
};


void
H::_step_opt_CS(MatrixXd *theta_x,
                MatrixXd *w,
                MatrixXd *x,
                MatrixXd *s,
                MatrixXd *c,
                int T){
    
	double l_t;
	double c_t;
    
    //	std::cout <<"w" << *w << "\n";
    
	for (int i=0; i< T; i++){
        
		x->col(0) = (*theta_x)(0,0) * s->col(1).array();
		x->col(1) = (*theta_x)(0,1) * s->col(4).array();
        
		//labor market results
		//if wage is less than realized wage - get hired, otherwise zero
		s->col(11) = (w->col(i*3).array() >= x->col(1).array()).cast<double>().array() *x->col(0).array();
		s->col(12) = w->col(i*3);
        
        //		auto FLAG_UPDATE = ((s->col(12) > 0.0).array() ||
        //							((s->col(12) == 0.0).array() &&
        //							(s->col(11) != 0.0).array()) ||
        //							((x->col(0) > 0.0).array() &&
        //							(s->col(11) == 0.0).array()));
        
		//TODO update inside , not coded for now
        //		if (param->wm_INSIDE_UPDATE){
        //		};
        
		//made consumption good decision
		x->col(2) = (*theta_x)(0,2) * s->col(0).array();
        
        
		//buy consumption goods , assume all goods are consumed in the previous
		//period, so no addition here
		for (int j=0; j<s->rows(); j++){
			if ((*w)(j,i*3 + 1) > 0.0){
				(*s)(j,3) = (*x)(j,2) / (*w)(j,i*3 + 1);
                
                
                
				//for speed calculate goal in this cycle
				//stores realization of a goal
				//as couldn't pass block need to calculate in place?
                
                //				(*c)(j,i) = _goal_t(&(s->block(1,s->cols(),j,0)));
			}else{
				(*s)(j,3) = 0.0;
			};
            
			if ((*s)(j,3) <= 0.0){
				c_t = -0.5;
//                c_t = 0.0;
			} else {
				c_t = (*s)(j,3);
			};
            
			l_t = (*s)(j,11);
            
            
			(*c)(j,i) = param->GOAL_T_theta(0,0) * log(1+c_t) + param->GOAL_T_theta(0,1) * (1-l_t);
//            (*c)(j,i) = param->GOAL_T_theta(0,0) * log(0.5+c_t) + param->GOAL_T_theta(0,1) * (1-l_t);

            
		};
        
		//account for money being spent
		s->col(0) = s->col(0).array() - w->col(i*3 + 1).array() * s->col(3).array();
        
        
		//labor payments
		s->col(0) = s->col(0).array() + s->col(11).array() * s->col(12).array();
        
		//dividend payment
		s->col(0) += w->col(i*3 + 2);
        
        //		std::cout << s->col(0) << "\n";
//        std::cout << *s << "\n";
        //		std::cout << *c << "\n";
        //		std::cout << *x << "\n";
        
		//skipped, because they are assigned, not added in a cycle
		//zero labor contracts
		//zero consumption good
        
	};
    
};


    
    
void
H::_bf(Eigen::MatrixXd* s, Eigen::MatrixXd* ret){
    (*ret)(0,0) = (*s)(0,0);
};

void
H::_bf(Eigen::Block<Eigen::MatrixXd>& s, Eigen::MatrixXd* ret){
    (*ret)(0,0) = s(0,0);
};


double
H::_V(Eigen::MatrixXd* s, Eigen::MatrixXd* phi_f, Eigen::MatrixXd& theta_V){
    _bf(s, phi_f);
    
    return (theta_V*(*phi_f).transpose()).sum();
    
};


double
H::_V(Eigen::Block<Eigen::MatrixXd>& s, Eigen::MatrixXd* phi_f, Eigen::MatrixXd& theta_V){
    _bf(s, phi_f);
    
    return (theta_V*(*phi_f).transpose()).sum();
    
};

Eigen::MatrixXd
H::_API_LM(Eigen::MatrixXd* s0, int N, int M){
    //Approximate policy iteration using linear models.
    //
    //p.407 ADP
    //
    
    //fix basis functions
    long n_theta_t = param->ADP_BF_N;
    
    //inner theta
    //theta for now and future
    Eigen::MatrixXd theta_V_n = Eigen::MatrixXd::Constant(1, n_theta_t, 1);
    
    //policy theta
    Eigen::MatrixXd theta_V_pi = theta_V;
    
    
    Eigen::MatrixXd s_n_m;
    Eigen::MatrixXd v_n_m;
    Eigen::MatrixXd x_n_m(1,3);
    Eigen::MatrixXd phi(1, param->ADP_BF_N);
    Eigen::MatrixXd phi1(1, param->ADP_BF_N);
    
    
    
    for (int n=0; n<N; n++){
        theta_V_n = theta_V_pi;
        
        s_n_m = MatrixXd::Zero(M+1, s0->cols());
        s_n_m.row(0) += s0->row(0);
        //v_n_m_T1 = 0
        v_n_m = Eigen::MatrixXd::Constant(M, n_theta_t, 0.0);
        
        
        for (int m=0; m<M; m++){
            
            Eigen::Block<Eigen::MatrixXd> s = s_n_m.block(m,0,1,(*s0).cols());
            
            
            //choose action, takes reference to s
            _opt_CS_s_V(s, x_n_m, theta_V_pi, param->opt_CS_N);
            
            s_n_m.block(m+1,0,1,(*s0).cols()) = s;
            
            Eigen::Block<Eigen::MatrixXd> s1 = s_n_m.block(m+1,0,1,(*s0).cols());
            s = s_n_m.block(m,0,1,(*s0).cols());
            
            
            //make step given action, to update state and return value function
            double _c = _c_theta_x(x_n_m, s1, theta_V_pi);
            
            double c_n_m = _c - _V(s1,&phi,theta_V_pi);
            
            _bf(s,&phi);
            _bf(s1,&phi1);
            Eigen::MatrixXd v_n_m = phi - param->BETA*phi1;
            
            //            std::cout << s_n_m.row(m) << std::endl;
            //            std::cout << s << std::endl;
            //
            //            std::cout << s_n_m.row(m+1) << std::endl;
            //            std::cout << theta_V_pi << std::endl;
            //            std::cout << theta_V_n << std::endl;
            
            
            //update theta_V_n
            theta_V_n = _RLS(theta_V_n,
                             v_n_m,
                             s,
                             m,
                             c_n_m);
            //assume that only 1 element in v and restrict it
            if (theta_V_n(0,0) > 1000){
                theta_V_n(0,0) = 1000;
            };
            if (theta_V_n(0,0) < 0.0){
                theta_V_n(0,0) = 0.01;
            };
            
        };
        
        theta_V_pi = theta_V_n;
    };
    
    //std::cout << theta_V_pi << std::endl;
    return theta_V_pi;
};

Eigen::MatrixXd
H::_RLS(Eigen::MatrixXd& theta, Eigen::MatrixXd& v, Eigen::Block<Eigen::MatrixXd>& s, int i, double c){
    
    //i - iteration if i = 0 - initialize B
    //for each theta_t:
    Eigen::MatrixXd B_n_1;
    
    if (i == 0){
        //identity matrix
        double e_B_0 = 0.0005;
        
        //B(n-1)
        B_n_1 = e_B_0 * Eigen::MatrixXd::Identity(param->ADP_BF_N, param->ADP_BF_N);
        
    }else{
        B_n_1 = B_n_1_RLS;
    };
    
    //container for basis functions
    Eigen::MatrixXd phi(1, param->ADP_BF_N);
    _bf(s,&phi);
    
    double e_n = c - (v.transpose()*theta)(0,0);
    
    double gamma_n = 1 + ((v.transpose()*B_n_1)*phi)(0,0);
    
    
    Eigen::MatrixXd B_n = B_n_1 - 1/gamma_n *((B_n_1*(phi*v.transpose()))*B_n_1);
    
    Eigen::MatrixXd theta_n = theta + 1/gamma_n * ((e_n*B_n_1)*phi);
    
    //    std::cout << 1/gamma_n * ((e_n*B_n_1)*phi) << std::endl;
    
    //    std::cout << theta_n << std::endl;
    
    B_n_1_RLS = B_n;
    
    
    
    return theta_n;
    
    
};




void
H::_opt_CS_s_V(Eigen::Block<Eigen::MatrixXd>& s0,  Eigen::MatrixXd& theta_x, Eigen::MatrixXd& theta_V, int N){
    
    //	N = 1;
    
	//gets matrix of thetas
	//prepare initial matrix
	MatrixXd m_s0 = MatrixXd::Zero(N, s0.cols());
    
	m_s0.rowwise() += s0.row(0);
    
    
	//allocate space to the matrix
	MatrixXd n_w;
	MatrixXd n_x;
	MatrixXd n_s = MatrixXd::Zero(N, s0.cols());
	MatrixXd n_c;
    MatrixXd v_bf_ret = MatrixXd::Zero(1, param->ADP_BF_N);
    MatrixXd n_v;
    
    
	//create matrix for the grid - small size here
	MatrixXd grid(3,3);
    
    create_decision_grid(grid,0);
    
    if (w.param->SIMULATION_MODE == "test"){
        Tools::mgrid_test(grid, &theta_x_all);
        
    } else{
        
        Tools::mgrid(grid, &theta_x_all);
    };
    
	//create random number generators for the implementation
    
	boost::normal_distribution<> nd_w1(s0(0,4), pow(s0(0,15),0.5));
	boost::normal_distribution<> nd_w2(s0(0,6), pow(s0(0,16),0.5));
	boost::normal_distribution<> nd_w3(s0(0,13), pow(s0(0,17),0.5));
    
	boost::variate_generator<boost::mt19937&,
    boost::normal_distribution<>> rng_w1(w.rng, nd_w1);
    
	boost::variate_generator<boost::mt19937&,
    boost::normal_distribution<>> rng_w2(w.rng, nd_w2);
    
	boost::variate_generator<boost::mt19937&,
    boost::normal_distribution<>> rng_w3(w.rng, nd_w3);
    
    
    
    if (w.param->SIMULATION_MODE != "test"){
        
        //draw random variables
        //fast realization
        n_w = MatrixXd::Zero(N, 3 * param->T_MAX);
        for (int i = 0; i < N; i++){
            for(int j=0; j< n_w.cols(); j++){
                if (j%3 == 0){
                    n_w(i,j) = rng_w1();
                } else {
                    if (j%3 == 1){
                        n_w(i,j) = rng_w2();
                    }else{
                        n_w(i,j) = rng_w3();
                    };
                    
                };
                n_w(i,j) = std::max(0.0, (double)n_w(i,j));
            };
        };
    };
    
	//matrix of all
	MatrixXd c_all(theta_x_all.cols(),1);
    
	for (int i = 0 ; i < theta_x_all.cols(); i++){
        if (w.param->SIMULATION_MODE == "test"){
            n_w = MatrixXd::Zero(N, 3 * param->T_MAX);
            for (int i = 0; i < N; i++){
                for(int j=0; j< n_w.cols(); j++){
                    if (j%3 == 0){
                        n_w(i,j) = Tools::get_normal(s0(0,4), pow(s0(0,15),0.5), w.myrng);
                    } else {
                        if (j%3 == 1){
                            n_w(i,j) = Tools::get_normal(s0(0,6), pow(s0(0,16),0.5), w.myrng);
                        }else{
                            n_w(i,j) = Tools::get_normal(s0(0,13), pow(s0(0,17),0.5), w.myrng);
                        };
                        
                    };
                    n_w(i,j) = std::max(0.0, (double)n_w(i,j));
                };
            };
        };
        
		//call estimator of results of steps, given theta_x
		theta_x.row(0) = theta_x_all.col(i).transpose();
		n_x = MatrixXd::Zero(N, 3);
		n_s = m_s0;
		n_c = MatrixXd::Zero(N, 1);
        n_v = MatrixXd::Zero(N, 1);
		_step_opt_CS(theta_x, &n_w, &n_x, &n_s, &n_c, &n_v, theta_V, &v_bf_ret, 1);
        
		//after n_c
		//average over n runs
        c_all(i,0) = (n_c + param->BETA*n_v).mean();
        
	};
    
    //    std::cout << (n_c + param->BETA*n_v).mean() << std::endl;
    
    
	//find all max elements
	double max = c_all(0,0);
	std::list<int> i_max;
	for (int i=0; i<c_all.rows(); i++){
		if (c_all(i,0) == max){
			i_max.push_back(i);
		}else{
			if (c_all(i,0) > max){
				max = c_all(i,0);
				i_max.clear();
				i_max.push_back(i);
			};
		}
	};
    
	//if more than 1 max - randomly pick
	long theta_max_i;
	if (i_max.size()>1){
		boost::random::uniform_int_distribution<> i_max_dist(0, i_max.size()-1);
		theta_max_i = i_max_dist(w.rng);
	}else{
		theta_max_i = i_max.front();
	};
    
    
	theta_x.row(0) = theta_x_all.col(theta_max_i).transpose();
    
    //    std::cout << theta_x << std::endl;
};

void
H::_step_opt_CS(MatrixXd& theta_x,
                MatrixXd *w,
                MatrixXd *x,
                MatrixXd *s,
                MatrixXd *c,
                MatrixXd* v,
                MatrixXd& theta_v,
                MatrixXd* ret,
                int T){
    
	double l_t;
	double c_t;

	//zero vector
	auto s_zero = MatrixXd::Zero(s->rows(),1);
    Eigen::Block<Eigen::MatrixXd> s_block = s->block(0,0,1,s->cols());;

    
    
    //	std::cout <<"w" << *w << "\n";
    
	for (int i=0; i< T; i++){
        
		x->col(0) = theta_x(0,0) * s->col(1).array();
		x->col(1) = theta_x(0,1) * s->col(4).array();
        
		//labor market results
		//if wage is less than realized wage - get hired, otherwise zero
		s->col(11) = (w->col(i*3).array() >= x->col(1).array()).cast<double>().array() *x->col(0).array();
		s->col(12) = w->col(i*3);
        
		//made consumption good decision
		x->col(2) = theta_x(0,2) * s->col(0).array();
        
        
		//buy consumption goods , assume all goods are consumed in the previous
		//period, so no addition here
		for (int j=0; j<s->rows(); j++){
			if ((*w)(j,i*3 + 1) > 0.0){
				(*s)(j,3) = (*x)(j,2) / (*w)(j,i*3 + 1);
                
                
                
				//for speed calculate goal in this cycle
				//stores realization of a goal
				//as couldn't pass block need to calculate in place?
                
                //				(*c)(j,i) = _goal_t(&(s->block(1,s->cols(),j,0)));
			}else{
				(*s)(j,3) = 0.0;
			};
            
			if ((*s)(j,3) <= 0.0){
				c_t = -0.5;
//                c_t = 0.0;
			} else {
				c_t = (*s)(j,3);
			};
            
			l_t = (*s)(j,11);
            
            
			(*c)(j,i) = param->GOAL_T_theta(0,0) * log(1+c_t) + param->GOAL_T_theta(0,1) * (1-l_t);
//            (*c)(j,i) = param->GOAL_T_theta(0,0) * log(0.5+c_t) + param->GOAL_T_theta(0,1) * (1-l_t);

            
		};
        
		//account for money being spent
		s->col(0) = s->col(0).array() - w->col(i*3 + 1).array() * s->col(3).array();
        
        
		//labor payments
		s->col(0) = s->col(0).array() + s->col(11).array() * s->col(12).array();

    
		//dividend payment
		s->col(0) += w->col(i*3 + 2);
        
        
        //zero out labor contracts
        s->col(11) = s_zero;
        s->col(12) = s_zero;
        //zero out consumption goods
        s->col(3) = s_zero;
        
        //add value function estimation
        for (int j=0; j < s->rows(); j++){
            s_block  = s->block(j,0,1,s->cols());
            
            (*v)(j,i) = _V(s_block,ret, theta_V);
        };
        


        
	};
    
};


double
H::_c_theta_x(Eigen::MatrixXd& theta_x, Eigen::Block<Eigen::MatrixXd>& s0, Eigen::MatrixXd& theta_V){
    
    int N = 1;
    double _c = 0.0;
    
    //    std::cout << s0 << std::endl;
    
	//gets matrix of thetas
	//prepare initial matrix
	MatrixXd m_s0 = MatrixXd::Zero(N, s0.cols());
    
	m_s0.rowwise() += s0.row(0);
    
    
	//allocate space to the matrix
	MatrixXd n_w;
	MatrixXd n_x = MatrixXd::Zero(N, 3);
	MatrixXd n_s = MatrixXd::Zero(N, s0.cols());
    MatrixXd v_bf_ret = MatrixXd::Zero(1, param->ADP_BF_N);
	MatrixXd n_c = MatrixXd::Zero(N, 1);;
    MatrixXd n_v = MatrixXd::Zero(N, 1);;
    
    
	//create random number generators for the implementation
    
	boost::normal_distribution<> nd_w1(s0(0,4), pow(s0(0,15),0.5));
	boost::normal_distribution<> nd_w2(s0(0,6), pow(s0(0,16),0.5));
	boost::normal_distribution<> nd_w3(s0(0,13), pow(s0(0,17),0.5));
    
	boost::variate_generator<boost::mt19937&,
    boost::normal_distribution<>> rng_w1(w.rng, nd_w1);
    
	boost::variate_generator<boost::mt19937&,
    boost::normal_distribution<>> rng_w2(w.rng, nd_w2);
    
	boost::variate_generator<boost::mt19937&,
    boost::normal_distribution<>> rng_w3(w.rng, nd_w3);
    
    
    
    if (w.param->SIMULATION_MODE != "test"){
        
        //draw random variables
        //fast realization
        n_w = MatrixXd::Zero(N, 3 * param->T_MAX);
        for (int i = 0; i < N; i++){
            for(int j=0; j< n_w.cols(); j++){
                if (j%3 == 0){
                    n_w(i,j) = rng_w1();
                } else {
                    if (j%3 == 1){
                        n_w(i,j) = rng_w2();
                    }else{
                        n_w(i,j) = rng_w3();
                    };
                    
                };
                n_w(i,j) = std::max(0.0, (double)n_w(i,j));
            };
        };
    };

    
    
    
    //redraw random in comparative with Python mode
    if (w.param->SIMULATION_MODE == "test"){
        n_w = MatrixXd::Zero(N, 3 * param->T_MAX);
        for (int i = 0; i < N; i++){
            for(int j=0; j< n_w.cols(); j++){
                if (j%3 == 0){
                    n_w(i,j) = Tools::get_normal(s0(0,4), pow(s0(0,15),0.5), w.myrng);
                } else {
                    if (j%3 == 1){
                        n_w(i,j) = Tools::get_normal(s0(0,6), pow(s0(0,16),0.5), w.myrng);
                    }else{
                        n_w(i,j) = Tools::get_normal(s0(0,13), pow(s0(0,17),0.5), w.myrng);
                    };
                    
                };
                n_w(i,j) = std::max(0.0, (double)n_w(i,j));
            };
        };
    };
    
    //call estimator of results of steps, given theta_x
    n_s = m_s0;
    _step_opt_CS(theta_x, &n_w, &n_x, &n_s, &n_c, &n_v, theta_V, &v_bf_ret, 1);
    _c = (n_c + param->BETA*n_v).mean();
    
    s0.row(0) = n_s.row(0);
    
    
    //    std::cout << s0 << std::endl;
    
    return _c;
    
};


void
H::opt_QL(){
    //find best action for the state
    
    //find current bin for the state
    _s0();
    std::vector<double> s_bin = _place_state_QL(s0);
    s_t_bin = s_bin;
    std::vector<int> Q_s;
    
    
    //create list of all indexes that correspond to current state
    for (int i = 0; i < s_bin.size(); i++){
        for (int j = 0; j< Q.cols(); j++){
            if (_is_bin(j, i, s_bin)){
                Q_s.push_back(j);
            };
            
        };
    };
    
    //find index with max Q
    std::vector<int> Q_max;
    double max = Q(4,Q_s.at(0));
    
    for (auto j:Q_s){
        if (Q(4,j) > max){
            max = Q(4,j);
            Q_max.clear();
            Q_max.push_back(j);
            
        }else{
            if (Q(4,j) == max){
                Q_max.push_back(j);
            };
        };
        
    };
    
    
    
    
    boost::random::uniform_01<> _exp_dist;
    double exp_rate = _exp_dist(w.rng);
    
    boost::random::uniform_int_distribution<> Q_dist(0, Q.cols()-1);
    long theta_max_i;
    
    
    if (exp_rate < param->QL_experimental_rate){
        theta_max_i = Q_dist(w.rng);
    }else{
        //pick randomly best action
        if (Q_max.size()>1){
            boost::random::uniform_int_distribution<> Q_max_dist(0, Q_max.size()-1);
            theta_max_i = Q_max.at(Q_max_dist(w.rng));
            
        }else{
            theta_max_i = Q_max.front();
        };
    };
    
    theta_x_t_bin.push_back(theta_max_i);
    
    //transform from index to decision
    theta_x.row(0) = Q.block(s_bin.size(),theta_max_i,3,1).transpose();
    
};

bool
H::_is_bin(int j, int i, std::vector<double>& s_bin){
    bool is_bin = false;
    
    if (i < s_bin.size()){
        is_bin = (Q(i,j) == Q_bin.at(i).at(s_bin.at(i))) && _is_bin(j, i+1, s_bin);
        //        std::cout << (Q(i,j) == s_bin.at(i)) << std::endl;
        //        std::cout << _is_bin(j, i+1, s_bin) << std::endl;
    }else{
        is_bin = true;
    };
    
    return is_bin;
};



void
H::wm_update_Q(){
    //updates Q matrix
    //reward
    double R_t1 = account->g_t.at(account->g_t.size()-2);
    
    //new state
    _s0();
    std::vector<double> s_bin = _place_state_QL(s0);
    std::vector<int> Q_s;
    
    
    //create list of all indexes that correspond to current state
    for (int i = 0; i < s_bin.size(); i++){
        for (int j = 0; j< Q.cols(); j++){
            if (_is_bin(j, i+1, s_bin)){
                Q_s.push_back(j);
            };
            
        };
    };
    
    //find index with max Q
    std::vector<int> Q_max;
    double max = Q(4,Q_s.at(0));
    
    for (auto j:Q_s){
        if (Q(4,j) > max){
            max = Q(4,j);
            Q_max.clear();
            Q_max.push_back(j);
            
        }else{
            if (Q(4,j) == max){
                Q_max.push_back(j);
            };
        };
        
    };
    
    //pick best action
    long theta_max_i;
    long max_i;
    if (Q_max.size()>1){
        boost::random::uniform_int_distribution<> Q_max_dist(0, Q_max.size()-1);
        max_i = Q_max_dist(w.rng);
        theta_max_i = Q_max.at(max_i);
    }else{
        theta_max_i = Q_max.front();
    };
    
    
    //update Q value
    Q(4,theta_x_t_bin.back()) = Q(4,theta_x_t_bin.back())
    + Q(5,theta_x_t_bin.back())
    * (R_t1 + param->BETA*Q(4,theta_max_i) - Q(4,theta_x_t_bin.back()));
    //    std::cout << Q << std::endl;
    
};



void
H::_set_QL(){
    //setup QL matrix
    //create matrix for the grid - small size here
	Q_grid = Eigen::MatrixXd::Zero(6,3);
    
    //money and stock of goods on hand
    Q_grid(0,0) = 0.0;
    Q_grid(0,1) = 10.0;
    Q_grid(0,2) = ((Q_grid(0,1) - Q_grid(0,0))/2);
    
    create_decision_grid(Q_grid, 1);

    //Q factors
    Q_grid(4,0) = 0.0;
    Q_grid(4,1) = 0.1;
    Q_grid(4,2) = 1.0;
    
    Q_grid(5,0) = param->QL_L; //speed of learning
    Q_grid(5,1) = Q_grid(5,0) + 0.1;
    Q_grid(5,2) = Q_grid(5,0) + 1.0;
    
    if (w.param->SIMULATION_MODE == "test"){
        Tools::mgrid_test(Q_grid, &Q);
        
    } else{
        
        Tools::mgrid(Q_grid, &Q);
    };
    
    create_bin_values(Q_grid, Q_bin);
    
    //update Q to push to hiring all of them
    for (int j = 0; j < Q.cols(); j++){
        if (Q(1,j) == Q_bin.at(1).back()){
            Q(4,j) = 0.05;
        };
    };
    
    
    //    std::cout << Q << std::endl;
    
    //    Tools::print_vector(Q_bin);
    
};


void
H::create_bin_values(Eigen::MatrixXd& _grid, std::vector<std::vector<double>>& _bin){
    long N_rows = _grid.rows();
	int n_steps_i;
	double n_steps_i_int;
    
	double val;
    std::vector<double> values_temp;
    
    //create matrix of bin values
    //create list of values
	for (int i=0; i < N_rows; i++){
		//n_steps
		n_steps_i_int = (_grid(i,1) - _grid(i,0))/_grid(i,2);
        
		n_steps_i = static_cast<int>(floor(n_steps_i_int)) + 1;
        
        
		for (int j=0; j<n_steps_i; j++){
			val = _grid(i,0) + _grid(i,2)*j;
			values_temp.push_back(val);
		};
        
		_bin.push_back(values_temp);
		values_temp.clear();
	};
    
};


std::vector<double>
H::_place_state_QL(Eigen::MatrixXd& s){
    std::vector<double> s_bin;
    
    //pick first and compare
    //push values for bins ina  vector
    
    bool FLAG_PLACED = false;
    
    int j = Q_bin.at(0).size()-1;
    while ((j>=0) && !FLAG_PLACED){
        
        if (s(0,0) >= Q_bin.at(0).at(j)){
            s_bin.push_back(j);
            FLAG_PLACED = true;
        };
        j -= 1;
    };
    
    return s_bin;
};

void
H::create_decision_grid(Eigen::MatrixXd& grid, int begin_index){
    int i = begin_index;
    
    if (param->GRID == "small"){
		grid(i,0) = 0.0;
		grid(i,1) = 1.1;
		grid(i,2) = 1.0;
		grid(i+1,0) = 0.8;
		grid(i+1,1) = 1.25;
		grid(i+1,2) = 0.2;
		grid(i+2,0) = 0.0;
		grid(i+2,1) = 1.0;
		grid(i+2,2) = 0.5;
	};
    
    if (param->GRID == "big"){
		grid(i,0) = 0.0;
		grid(i,1) = 1.1;
		grid(i,2) = 1.0;
		grid(i+1,0) = 0.1;
		grid(i+1,1) = 2.1;
		grid(i+1,2) = 0.45;
		grid(i+2,0) = 0.0;
		grid(i+2,1) = 1.0;
		grid(i+2,2) = 0.5;
        
	};
    if (param->GRID == "test"){
        
        grid(i,0) = 1.0;
        grid(i,1) = 1.1;
        grid(i,2) = 1.0;
        
        grid(i+1,0) = 0.8;
        grid(i+1,1) = 1.2;
        grid(i+1,2) = 1.0;
        
        grid(i+2,0) = 0.5;
        grid(i+2,1) = 1.1;
        grid(i+2,2) = 1.0;
        
    };
    
    
};

void
H::_set_mRE(){
    //
    mRE_grid = Eigen::MatrixXd::Zero(4,3);
    create_decision_grid(mRE_grid, 0);
    mRE_grid(3,0) = 1.0;
    mRE_grid(3,1) = 1.1;
    mRE_grid(3,2) = 1.0;
    
    if (w.param->SIMULATION_MODE == "test"){
        Tools::mgrid_test(mRE_grid, &mRE);
        
    } else{
        
        Tools::mgrid(mRE_grid, &mRE);
    };
    
    
    
    create_bin_values(mRE_grid, mRE_bin);
    
    //update mRE to push to hiring all of them
    for (int j = 0; j < mRE.cols(); j++){
        if (mRE(0,j) == mRE_bin.at(0).back()){
            mRE(3,j) = 1.05;
        };
    };
    
};



void
H::opt_mRE(){
    //form matrix of probabilities
    double sum_prob = 0.0;
    
    for (int j=0; j<mRE.cols();j++){
        sum_prob += exp(mRE(3,j)/param->mRE_T);
    };
    
    std::vector<double> probs;
    
    for (int j=0; j<mRE.cols();j++){
        probs.push_back(exp(mRE(3,j)/param->mRE_T)/sum_prob);
    };
    
    
    //roll weighted dice
    boost::random::discrete_distribution<> dist(probs.begin(), probs.end());
    
    
    int theta_i = dist(w.rng);
    
    theta_x_t_bin.push_back(theta_i);
    
    //transform from index to decision
    theta_x.row(0) = mRE.block(0,theta_i,3,1).transpose();
    
};



void
H::wm_update_mRE(){
    //adjustment matrix
    Eigen::MatrixXd _E_re  = Eigen::MatrixXd::Zero(1, mRE.cols());
    
    //reward
    double R_t1 = account->g_t.at(account->g_t.size()-2);
    
    //adjustment
    _E_re = mRE.row(3) * (param->mRE_EPSILON/(mRE.cols()-1));
    
    //include reward
    _E_re(0,theta_x_t_bin.back()) = R_t1 * (1 - param->mRE_EPSILON);
    
    //update whole matrix
    mRE.row(3) = mRE.row(3) * (1 - param->mRE_PHI) + _E_re;
    

};









bool
H::_PS_accept_PO(MessagePSSendPayment* mes){
	return true;
};



ContractBDt0*
H::_PS_contract(){
	return ass_asCBDt0.front();
};




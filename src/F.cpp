/*
 * F.cc
 *
 *  Created on: Apr 13, 2014
 *      Author: wilfeli
 */
#include <list>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <math.h>
#include <string>

#include "Tools.h"
#include "F.h"
#include "W.h"
#include "Contract.h"
#include "Market.h"
#include "H.h"


//using std::shared_ptr;
//using std::unique_ptr;
using Eigen::MatrixXd;
using namespace boost;
//using Tools::mgrid;

//will pass w as itself, it will be converted to reference
F::F(W &w_, int id_):w(w_),
					l_bid(0.0, 0.0),
					c_ask(0.0, 0.0),
					id(id_),
					s0(1,19),
					theta_x(1,4){
	w.AddF(this);

	//create storage for goods
	ass_asGC.push_back(new GoodC(this, 0.0));

	type = "FGC";

	l_bid.issuer = this;
	c_ask.issuer = this;

	//create wm
	ExpectationBackward* exp = new ExpectationBackward(1.0, 10, 0.5);
	wm["MasCHK"] = exp;

	//create wm
	exp = new ExpectationBackward(1.0, 10, 0.5);
	wm["MasGC"] = exp;


	//parameters
	param = new Parameters();
	param->init();


	//accounting
	account = new AccountF();
                        
    

};

void
F::init(){
    //ADP setup
    theta_V = Eigen::MatrixXd::Constant(1, 1, param->ADP_BF_N);
    B_n_1_RLS = Eigen::MatrixXd::Zero(param->ADP_BF_N, param->ADP_BF_N);
    
    //QL setup
    _set_QL();
    
    
    //mRe setup
    _set_mRE();

};

double
F::get_status(){
	return status;
};

void
F::ac_W_begin_step(){
	//begin step - pay off obligations
	//check that has enough money to pay all
	double q_all = 0.0;
	double pq = 0.0;
	MessagePSSendPayment* inf;
	bool PAY_DIV = false;

	for (auto c:ass_asCHK){
		if(c->t_end == w.t - 1){
			q_all += c->q * c->p;
		};
	};

	bool PAY_ALL = false;

	if (q_all <= ass_asCBDt0.front()->q){
		PAY_ALL = true;
	};

	for (auto c:ass_asCHK){

		if(c->t_end == w.t - 1){

			pq = PAY_ALL?c->q * c->p:(c->q * c->p) * ass_asCBDt0.front()->q /q_all;


			//make payment
			inf = new MessagePSSendPayment(c->issuer, pq ,"asCHK_b");
			inf->sender = this;

			ass_asCBDt0.front()->issuer->_PS_accept_PO(inf);

			delete inf;

		};
	};


	if (w.t > 0.0){

		//check last period profit
		q_all = (account->profit.back() > 0.0)? account->profit.back() * param->DIV_SHARE: 0.0;

		if (q_all > 0.0) {
			PAY_DIV = true;
			if (q_all <= ass_asCBDt0.front()->q){
				PAY_ALL = PAY_ALL && true;
			} else{
				PAY_ALL = false;
			};
		};

		//pay off dividends
		if (PAY_ALL  && PAY_DIV){
			pq = q_all/w.NH;
			for (auto a:w.Hs){
				//make payment
				inf = new MessagePSSendPayment(pq ,static_cast<std::string>("FI_div"));
				inf->getter = a; //dynamic_cast<IPSAgent*>(a);
				inf->sender = this;

				ass_asCBDt0.front()->issuer->_PS_accept_PO(inf);

				delete inf;
			};

		};
	};



	//checks for bankruptcy
	if ((ass_asCBDt0.front()->q < 0.0) || !PAY_ALL){
		//calls bankruptcy
		MessageBankruptcy* mes = new MessageBankruptcy(this);
		w.ac_LS_bankruptcy(mes);

		delete mes;

	};


	//clear contracts
	ass_asCHK.erase(std::remove_if(ass_asCHK.begin(), ass_asCHK.end(),
	                       [&](ContractHK* x) -> bool { return (x->t_end <= w.t); }),
			ass_asCHK.end());

//	for (auto c:ass_asCHK){
//		std::cout << true << std::endl;
//	};

	account->ac_W_initialize_step();

};

void
F::ac_W_begin_step(MessageStatus* mes){
	if (mes->status <= 0.0){
		account->ac_W_initialize_step(mes);
	};
};


void
F::ac_W_end_step(){
	// calls accounting to calculate profit
	_s0();
	account->ac_W_end_step(&s0);

};


double
F::_mas_q_sell(MessageMarketCCheckAsk* inf){
	//checks how much to sell
	//return q from ask

	double q_sell = 0.0;

	if (inf->p_eq >= c_ask.p){
		q_sell = c_ask.q;
	};

	return q_sell;
};

double
F::_mas_q_buy(MessageMarketLCheckBid* inf){
	//checks how much to sell
	//return q from ask

	double q_buy = 0.0;

	if (inf->p_eq <= l_bid.p){
		q_buy = l_bid.q;
	};

	return q_buy;
};


void
F::_buy_asCHK(ContractHK* c){
	//add to contract
	ass_asCHK.push_back(c);


	//update ask
	l_bid.q -= c->q;

	//message for accounting
	Message* mes = new Message("buy_asCHK");

	//call accounting
	account->ac_get_inf(mes, c);


	delete mes;

};

void
F::wm_update_k_s(){

	if (w.t > 0.0){

		MatrixXd w0 = wm_w0_tbeg();
		ExpectationBackward* exp_i;
		std::vector<double> s_t_i;

		//update k_s if conditions are met
//		if ((w0(0,1) > 0.0) || ((w0(0,1) == 0.0) && (w0(0,0) != 0.0)) || ((w0(0,4) > 0.0) && (w0(0,0) == 0.0))){
			s_t_i.push_back(w0(0,0));
			s_t_i.push_back(w0(0,1));

			exp_i = wm["MasCHK"];

			exp_i->s_t.push_back(s_t_i);

			wm_update_expectation_backward(exp_i, 1);
//		};

		s_t_i.clear();

		s_t_i.push_back(w0(0,2));
		s_t_i.push_back((w0(0,2) > 0.0)? w0(0,3)/w0(0,2):w0(0,5));


		exp_i = wm["MasGC"];

		exp_i->s_t.push_back(s_t_i);

		wm_update_expectation_backward(exp_i, 1);

		s_t_i.clear();
        
        
        if (param->opt_TYPE == "QL"){
            wm_update_Q();
            
        };
        
        if (param->opt_TYPE == "mRE"){
            wm_update_mRE();
            
        };
        
	};

};


Eigen::MatrixXd
F::wm_w0_tbeg(){
	MatrixXd w0(1,8);

	long life_length = account->asCHK_q.size();

	w0(0,0) = account->asCHK_q.at(life_length - 2);
	w0(0,1) = w.ml->market_price.at(w.t - 1)->p;
	w0(0,2) = account->sales_q.at(life_length - 2);
	w0(0,3) = account->sales_pq.at(life_length - 2);

	w0(0,4) = l_bid.q_t_0;
	w0(0,5) = w.mc->market_price.at(w.t - 1)->p;
	w0(0,6) = c_ask.q_t_0;
	w0(0,7) = account->profit.at(life_length - 2);

//	std::cout << w0 << std::endl;

	return w0;
};

void
F::wm_update_expectation_backward(ExpectationBackward* exp_i, int p_position){


	std::vector<double> p;

	int i_s_t_1_max = 0;
	long i_s_t_max = 0;

	//depending on the length of the wm and accumulated number of prices
	if (std::isinf(param->WM_LENGTH)){
		i_s_t_1_max = exp_i->s_t_1["n"];
		i_s_t_max = exp_i->s_t.size();
	}else{
		i_s_t_1_max = std::min(std::max(param->WM_LENGTH - exp_i->s_t.size(), 0.0), exp_i->s_t_1["n"]);
		i_s_t_max = std::min(param->WM_LENGTH, static_cast<double>(exp_i->s_t.size()));
	};

	for (int i = 0; i < i_s_t_1_max; i++){
		p.push_back(exp_i->s_t_1["mu"]);
	};

	//push other prices
	for (std::size_t i = (exp_i->s_t.size() - i_s_t_max);i < exp_i->s_t.size();i++){
		p.push_back(exp_i->s_t[i][p_position]);
	};

//	Tools::print_vector(p);

	Eigen::MatrixXd mean_variance = wm_mean_variance(p);

	exp_i->mu = mean_variance(0,0);
	exp_i->n += 1.0;
	exp_i->v = mean_variance(1,0);
    
    if (exp_i->v <= 0.0){
        exp_i->v = exp_i->mu * 0.01;
    };


};



Eigen::MatrixXd
F::wm_mean_variance(std::vector<double> &x){
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
F::_s0(){
	//money holdings
	s0(0,0) = ass_asCBDt0.front()->q;

	s0(0,1) = 0.0;
	//human capital for production
	for (auto c:ass_asCHK){
		s0(0,1) += c->q;
	};

	//amount of good c
	s0(0,2) = ass_asGC.front()->q;

	//expected price on labor market
	s0(0,3) = wm["MasCHK"]->mu;

	//number of observations
	s0(0,4) = wm["MasCHK"]->n;

	//expected price on goods market
	s0(0,5) = wm["MasGC"]->mu;

	//number of observations
	s0(0,6) = wm["MasGC"]->n;

	//p*q - total costs
	s0(0,7) = account->cost_p * account->cost_q;

	//cost - q
	s0(0,8) = account->cost_q;

	//current period labor contracts
	s0(0,9) = 0.0;

	//current period labor price
	s0(0,10) = 0.0;

	//current period sales q
	s0(0,11) = 0.0;

	//current period sales p
	s0(0,12) = 0.0;

	//current period production
	s0(0,13) = 0.0;

	//current profit
	s0(0,14) = 0.0;

	//state dead/alive
	s0(0,15) = status;

	//share of dividends
	s0(0,16) = param->DIV_SHARE;

	//variance of labor market
	s0(0,17) = wm["MasCHK"]->v;

	//variance of goods market
	s0(0,18) = wm["MasGC"]->v;


	return s0;
};


Bid*
F::ac_L(){
	//update numbers in bid
	l_bid.demand_curve = Eigen::MatrixXd::Zero(1,2);
	l_bid.demand_curve << l_bid.p , l_bid.q;

	l_bid.p_t_0 = l_bid.p;
	l_bid.q_t_0 = l_bid.q;



	//return reference to the bid
	return &l_bid;
};

Ask*
F::ac_C(){
	//update numbers in ask
	c_ask.q = _s0()(0,2) * c_ask.q_share;

	c_ask.p_t_0 = c_ask.p;
	c_ask.q_t_0 = c_ask.q;
	c_ask.q_share_t_0 = c_ask.q_share;

	//update numbers in ask
	c_ask.supply_curve = Eigen::MatrixXd::Zero(1,2);
	c_ask.supply_curve << c_ask.p , c_ask.q;


//	std::cout << ass_asCHK.size() << "\n";
//
//	std::cout << c_ask.supply_curve << "\n";

	//return reference to the bid
	return &c_ask;
};


bool
F::_sell_asGC(MessageMarketCSellGoodC* mes){
	//decrease offer
	c_ask.q -= mes->q;

	//sell actual good
	ass_asGC.back()->q -= mes->q;


	MessageGoodC* mgc = new MessageGoodC(mes->q, mes->p);


	//send goods to the buyer
	mes->buyer->_buy_asGC(mgc);
    
    //delete message
    delete mgc;


	//account sale
	account->ac_get_inf(mes);

	return true;


};



void
F::ac_W(MessageMakeDec* mes){

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
	l_bid.q = floor(theta_x(0,0));
	l_bid.p = theta_x(0,1) * s0(0,3);

	c_ask.q_share = theta_x(0,2);
	c_ask.p = theta_x(0,3) * s0(0,5);



};


void
F::ac_W(MessageF_F* mes){
	//calls production function
	_s0();
	mes->q = _F_F(&s0);
	ass_asGC.back()->q += mes->q;
	account->ac_get_inf(mes);

};

double
F::_F_F(Eigen::MatrixXd* s){
	double l;

	l = (*s)(0,1);

	return param->F_F_theta(0,0) * pow(l, param->F_F_theta(0,1));
};



void
F::opt_CS(MatrixXd* s0,  MatrixXd* theta_x, int N){

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
	MatrixXd grid(4,3);

    create_decision_grid(grid,0);

    if (w.param->SIMULATION_MODE == "test"){
        Tools::mgrid_test(grid, &theta_x_all);

    } else{

        Tools::mgrid(grid, &theta_x_all);
    };

	//create random number generators for the implementation

	boost::normal_distribution<> nd_w1((*s0)(0,3), pow((*s0)(0,17),0.5));
	boost::normal_distribution<> nd_w2((*s0)(0,5), pow((*s0)(0,18),0.5));

	boost::variate_generator<boost::mt19937&,
	                           boost::normal_distribution<>> rng_w1(w.rng, nd_w1);

	boost::variate_generator<boost::mt19937&,
	                           boost::normal_distribution<>> rng_w2(w.rng, nd_w2);


    if (w.param->SIMULATION_MODE != "test"){

        //draw random variables
        n_w = MatrixXd::Zero(N, 2 * param->T_MAX);
        for (int i = 0; i < N; i++){
            for(int j=0; j< n_w.cols(); j++){
                if (j%2 == 0){
                    n_w(i,j) = rng_w1();
                } else {
                    n_w(i,j) = rng_w2();
                };
                n_w(i,j) = std::max(0.0, (double)n_w(i,j));
            };
        };
    };
    

	//matrix of all
	MatrixXd c_all(theta_x_all.cols(),1);

	for (int i = 0 ; i < theta_x_all.cols(); i++){
        //redraw random in comparative with Python mode
        if (w.param->SIMULATION_MODE == "test"){
            n_w = MatrixXd::Zero(N, 2 * param->T_MAX);
            for (int k = 0; k < N; k++){
                for(int j=0; j< n_w.cols(); j++){
                    if (j%2 == 0){
                        n_w(k,j) = Tools::get_normal((*s0)(0,3), pow((*s0)(0,17),0.5), w.myrng);
                    } else {
                        n_w(k,j) = Tools::get_normal((*s0)(0,5), pow((*s0)(0,18),0.5), w.myrng);
                    };
                    n_w(k,j) = std::max(0.0, (double)n_w(k,j));
                };
            };
        };
        
        
		//call estimator of results of steps, given theta_x
		theta_x->row(0) = theta_x_all.col(i).transpose();
		n_x = MatrixXd::Zero(N, 4);
		n_s = m_s0;
		n_c = MatrixXd::Zero(N, param->T_MAX);
		_step_opt_CS(theta_x, &n_w, &n_x, &n_s, &n_c, param->T_MAX);

//        std::cout << n_c << std::endl;
//        std::cout << theta_x->row(0) << std::endl;
        
		//after n_c
		//average over n runs
		c_all(i,0) = (n_c * param->BETA_T).mean();

	};

//	std::cout << n_c << "\n";

//	std::cout << c_all.transpose() << "\n";

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
};


void
F::_step_opt_CS(MatrixXd *theta_x,
		MatrixXd *w,
		MatrixXd *x,
		MatrixXd *s,
		MatrixXd *c,
                int T){

	//sets decisions for labor market
	//zero vector
	auto s_zero = MatrixXd::Zero(s->rows(),1);



	for (int i=0; i< T; i++){


		x->col(0) = Eigen::MatrixXd::Constant(x->rows(), 1, (*theta_x)(0,0));
		x->col(1) = (*theta_x)(0,1) * s->col(3).array();
        
//        std::cout << *x << std::endl;

		//labor market results
		//if wage is less than realized wage - zero, otherwise hire all desired number

		(*s).col(9) = ((*w).col(i*2).array() <= x->col(1).array()).cast<double>().array() *x->col(0).array();
		(*s).col(10) = (*w).col(i*2);
        
        //cast to whole numbers
        for (int j=0; j < s->rows(); j++){
            (*s)(j,9) = floor((*s)(j,9));
        };
        
//        std::cout << *s << std::endl;



		//current production
		s->col(13) = s->col(9).array().pow(param->F_F_theta(0,1))*param->F_F_theta(0,0);
		s->col(2) += s->col(13);

//        std::cout << *s << std::endl;
        
		//update for cost structure
		s->col(7) = s->col(7).array() + ((*s).col(9).array() * (*s).col(10).array());
		s->col(8) += s->col(13);

//        std::cout << *s << std::endl;

		//sets decision for goods market
		x->col(2) = (*theta_x)(0,2) * s->col(2).array();
		x->col(3) = Eigen::MatrixXd::Constant(x->rows(), 1, (*theta_x)(0,3) * (*s)(0,5));

//        std::cout << *x << std::endl;
        
		//sales
		(*s).col(11) = ((*w).col(i*2 + 1).array() >= x->col(3).array()).cast<double>().array() *x->col(2).array();
		(*s).col(12) = (*w).col(i*2 + 1);

//        std::cout << *s << std::endl;

		//money
		s->col(0) = s->col(0).array() + (s->col(11).array() * s->col(12).array());
        
//        std::cout << *s << std::endl;

		//inventory
		s->col(2) -= s->col(11);
        
//        std::cout << *s << std::endl;

		double atc = 0.0;

		//updates current profit
		for (int j=0; j < s->rows(); j++){

			if ((*s)(j,8)>0.0){
				atc = (*s)(j,7) / (*s)(j,8);
                
                
                if (param->ACCOUNTING_TYPE=="modern"){
                    //update profit to cost of production
                    (*s)(j,14) -= atc * (*s)(j,11);
                };
                if (param->ACCOUNTING_TYPE=="classic"){
                    //update profit to cost of production
                    (*s)(j,14) -= (*s)(j,9) * (*s)(j,10);
                };
               
                
				//update cost of production
				(*s)(j,7) -= atc * (*s)(j,11);
				(*s)(j,8) -= (*s)(j,11);
			};

			//update profit to sales
			(*s)(j,14) += (*s)(j,11) * (*s)(j,12);

			atc = 0.0;
		};

//        std::cout << *s << std::endl;

		//labor payment
		s->col(0) = s->col(0).array() -  (s->col(9).array() * s->col(10).array());
        
//        std::cout << *s << std::endl;

		//dividends
		s->col(0) = s->col(0).array() - (s->col(14).array() > s_zero.array()).cast<double>().array() * s->col(14).array() * s->col(16).array();

//        std::cout << *s << std::endl;

		//zero out labor contract
		s->col(9) = s_zero;
		s->col(10) = s_zero;


        //dead or alive
		s->col(15) = (s->col(0).array() >= s_zero.array()).cast<double>().array() * s->col(15).array();

//        std::cout << *s << std::endl;
        
		//stores realized goal
		c->col(i) = (s->col(15).array() > s_zero.array()).cast<double>().array() * s->col(14).array();

//        std::cout << *c << std::endl;

		//zero out current profit
		s->col(14) = s_zero;


	};


};



void
F::_bf(Eigen::MatrixXd* s, Eigen::MatrixXd* ret){
    (*ret)(0,0) = (*s)(0,0);
};

void
F::_bf(Eigen::Block<Eigen::MatrixXd>& s, Eigen::MatrixXd* ret){
    (*ret)(0,0) = s(0,0);
};


double
F::_V(Eigen::MatrixXd* s, Eigen::MatrixXd* phi_f, Eigen::MatrixXd& theta_V){
    _bf(s, phi_f);

    return (theta_V*(*phi_f).transpose()).sum();

};


double
F::_V(Eigen::Block<Eigen::MatrixXd>& s, Eigen::MatrixXd* phi_f, Eigen::MatrixXd& theta_V){
    _bf(s, phi_f);
    
    return (theta_V*(*phi_f).transpose()).sum();
    
};


Eigen::MatrixXd
F::_API_LM(Eigen::MatrixXd* s0, int N, int M){
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
    Eigen::MatrixXd x_n_m(1,4);
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
F::_RLS(Eigen::MatrixXd& theta, Eigen::MatrixXd& v, Eigen::Block<Eigen::MatrixXd>& s, int i, double c){

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




double
F::_c_theta_x(Eigen::MatrixXd& theta_x, Eigen::Block<Eigen::MatrixXd>& s0, Eigen::MatrixXd& theta_V){
    
    int N = 1;
    double _c = 0.0;
    
//    std::cout << s0 << std::endl;
    
	//gets matrix of thetas
	//prepare initial matrix
	MatrixXd m_s0 = MatrixXd::Zero(N, s0.cols());
    
	m_s0.rowwise() += s0.row(0);
    
    
	//allocate space to the matrix
	MatrixXd n_w;
	MatrixXd n_x = MatrixXd::Zero(N, 4);
	MatrixXd n_s = MatrixXd::Zero(N, s0.cols());
    MatrixXd v_bf_ret = MatrixXd::Zero(1, param->ADP_BF_N);
	MatrixXd n_c = MatrixXd::Zero(N, 1);;
    MatrixXd n_v = MatrixXd::Zero(N, 1);;

    
    //create random number generators for the implementation
    
	boost::normal_distribution<> nd_w1(s0(0,3), pow(s0(0,17),0.5));
	boost::normal_distribution<> nd_w2(s0(0,5), pow(s0(0,18),0.5));
    
	boost::variate_generator<boost::mt19937&,
    boost::normal_distribution<>> rng_w1(w.rng, nd_w1);
    
	boost::variate_generator<boost::mt19937&,
    boost::normal_distribution<>> rng_w2(w.rng, nd_w2);
    
    
    if (w.param->SIMULATION_MODE != "test"){
        
        //draw random variables
        n_w = MatrixXd::Zero(N, 2 * 1);
        for (int i = 0; i < N; i++){
            for(int j=0; j< n_w.cols(); j++){
                if (j%2 == 0){
                    n_w(i,j) = rng_w1();
                } else {
                    n_w(i,j) = rng_w2();
                };
                n_w(i,j) = std::max(0.0, (double)n_w(i,j));
            };
        };
    };


    

    //redraw random in comparative with Python mode
    if (w.param->SIMULATION_MODE == "test"){
        n_w = MatrixXd::Zero(N, 2 * param->T_MAX);
        for (int k = 0; k < N; k++){
            for(int j=0; j< n_w.cols(); j++){
                if (j%2 == 0){
                    n_w(k,j) = Tools::get_normal(s0(0,3), pow(s0(0,17),0.5), w.myrng);
                } else {
                    n_w(k,j) = Tools::get_normal(s0(0,5), pow(s0(0,18),0.5), w.myrng);
                };
                n_w(k,j) = std::max(0.0, (double)n_w(k,j));
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
F::_opt_CS_s_V(Eigen::Block<Eigen::MatrixXd>& s0,  Eigen::MatrixXd& theta_x, Eigen::MatrixXd& theta_V, int N){
    
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
	MatrixXd grid(4,3);
    
    create_decision_grid(grid,0);
    
    if (w.param->SIMULATION_MODE == "test"){
        Tools::mgrid_test(grid, &theta_x_all);
        
    } else{
        
        Tools::mgrid(grid, &theta_x_all);
    };
    
	//create random number generators for the implementation
    
	boost::normal_distribution<> nd_w1(s0(0,3), pow(s0(0,17),0.5));
	boost::normal_distribution<> nd_w2(s0(0,5), pow(s0(0,18),0.5));
    
	boost::variate_generator<boost::mt19937&,
    boost::normal_distribution<>> rng_w1(w.rng, nd_w1);
    
	boost::variate_generator<boost::mt19937&,
    boost::normal_distribution<>> rng_w2(w.rng, nd_w2);
    
    
    if (w.param->SIMULATION_MODE != "test"){
        
        //draw random variables
        n_w = MatrixXd::Zero(N, 2 * 1);
        for (int i = 0; i < N; i++){
            for(int j=0; j< n_w.cols(); j++){
                if (j%2 == 0){
                    n_w(i,j) = rng_w1();
                } else {
                    n_w(i,j) = rng_w2();
                };
                n_w(i,j) = std::max(0.0, (double)n_w(i,j));
            };
        };
    };
    
    
	//matrix of all
	MatrixXd c_all(theta_x_all.cols(),1);
    
	for (int i = 0 ; i < theta_x_all.cols(); i++){
        //redraw random in comparative with Python mode
        if (w.param->SIMULATION_MODE == "test"){
            n_w = MatrixXd::Zero(N, 2 * 1);
            for (int k = 0; k < N; k++){
                for(int j=0; j< n_w.cols(); j++){
                    if (j%2 == 0){
                        n_w(k,j) = Tools::get_normal(s0(0,3), pow(s0(0,17),0.5), w.myrng);
                    } else {
                        n_w(k,j) = Tools::get_normal(s0(0,5), pow(s0(0,18),0.5), w.myrng);
                    };
                    n_w(k,j) = std::max(0.0, (double)n_w(k,j));
                };
            };
        };
        
        
		//call estimator of results of steps, given theta_x
		theta_x.row(0) = theta_x_all.col(i).transpose();
		n_x = MatrixXd::Zero(N, 4);
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
F::create_decision_grid(Eigen::MatrixXd& grid, int begin_index){
    int i = begin_index;
    
    if (param->GRID == "small"){
		grid(i,0) = 0.0;
		grid(i,1) = w.NH;
		grid(i,2) = ((grid(i,1) - grid(i,0))/4);
		grid(i+1,0) = 0.8;
		grid(i+1,1) = 1.25;
		grid(i+1,2) = 0.2;
		grid(i+2,0) = 0.0;
		grid(i+2,1) = 1.0;
		grid(i+2,2) = 0.5;
		grid(i+3,0) = 0.8;
		grid(i+3,1) = 1.25;
		grid(i+3,2) = 0.2;
	};
    
    if (param->GRID == "big"){
		grid(i,0) = 0.0;
		grid(i,1) = w.NH;
		grid(i,2) = (grid(i,1) - grid(i,0))/4;
		grid(i+1,0) = 0.1;
		grid(i+1,1) = 2.1;
		grid(i+1,2) = 0.45;
		grid(i+2,0) = 0.0;
		grid(i+2,1) = 1.0;
		grid(i+2,2) = 0.5;
		grid(i+3,0) = 0.1;
		grid(i+3,1) = 2.1;
		grid(i+3,2) = 0.45;
        
	};
    if (param->GRID == "test"){
        
        //        grid(0,0) = w.NH;
        //        grid(0,1) = w.NH + 0.5;
        //        grid(0,2) = 1.0;
        //        grid(1,0) = 1.0;
        //        grid(1,1) = 1.1;
        //        grid(1,2) = 1.0;
        //        grid(2,0) = 1.0;
        //        grid(2,1) = 1.1;
        //        grid(2,2) = 1.0;
        //        grid(3,0) = 1.0;
        //        grid(3,1) = 1.1;
        //        grid(3,2) = 1.0;
        
        grid(i,0) = 2.5;
        grid(i,1) = 2.5 + 1.5;
        grid(i,2) = 1.0;
        
        grid(i+1,0) = 1.0;
        grid(i+1,1) = 1.1;
        grid(i+1,2) = 1.0;
        
        grid(i+2,0) = 0.5;
        grid(i+2,1) = 1.1;
        grid(i+2,2) = 1.0;
        
        grid(i+3,0) = 1.2;
        grid(i+3,1) = 1.3;
        grid(i+3,2) = 1.0;
        
        
        
    };

    
};

void
F::_step_opt_CS(MatrixXd& theta_x,
                MatrixXd *w,
                MatrixXd *x,
                MatrixXd *s,
                MatrixXd *c,
                MatrixXd* v,
                MatrixXd& theta_v,
                MatrixXd* ret,
                int T){
    
	//sets decisions for labor market
	//zero vector
	auto s_zero = MatrixXd::Zero(s->rows(),1);
    Eigen::Block<Eigen::MatrixXd> s_block = s->block(0,0,1,s->cols());;
    
    
    
	for (int i=0; i< T; i++){
        
        
		x->col(0) = Eigen::MatrixXd::Constant(x->rows(), 1, theta_x(0,0));
		x->col(1) = theta_x(0,1) * s->col(3).array();
        
        //        std::cout << *x << std::endl;
        
		//labor market results
		//if wage is less than realized wage - zero, otherwise hire all desired number
        
		(*s).col(9) = ((*w).col(i*2).array() <= x->col(1).array()).cast<double>().array() *x->col(0).array();
		(*s).col(10) = (*w).col(i*2);
        
        //cast to whole numbers
        for (int j=0; j < s->rows(); j++){
            (*s)(j,9) = floor((*s)(j,9));
        };
        
        //        std::cout << *s << std::endl;
        
        
        
		//current production
		s->col(13) = s->col(9).array().pow(param->F_F_theta(0,1))*param->F_F_theta(0,0);
		s->col(2) += s->col(13);
        
        //        std::cout << *s << std::endl;
        
		//update for cost structure
		s->col(7) = s->col(7).array() + ((*s).col(9).array() * (*s).col(10).array());
		s->col(8) += s->col(13);
        
        //        std::cout << *s << std::endl;
        
		//sets decision for goods market
		x->col(2) = theta_x(0,2) * s->col(2).array();
		x->col(3) = Eigen::MatrixXd::Constant(x->rows(), 1, theta_x(0,3) * (*s)(0,5));
        
        //        std::cout << *x << std::endl;
        
		//sales
		(*s).col(11) = ((*w).col(i*2 + 1).array() >= x->col(3).array()).cast<double>().array() *x->col(2).array();
		(*s).col(12) = (*w).col(i*2 + 1);
        
        //        std::cout << *s << std::endl;
        
		//money
		s->col(0) = s->col(0).array() + (s->col(11).array() * s->col(12).array());
        
        //        std::cout << *s << std::endl;
        
		//inventory
		s->col(2) -= s->col(11);
        
        //        std::cout << *s << std::endl;
        
		double atc = 0.0;
        
		//updates current profit
		for (int j=0; j < s->rows(); j++){
            
			if ((*s)(j,8)>0.0){
				atc = (*s)(j,7) / (*s)(j,8);

                
                if (param->ACCOUNTING_TYPE=="modern"){
                    //update profit to cost of production
                    (*s)(j,14) -= atc * (*s)(j,11);
                };
                if (param->ACCOUNTING_TYPE=="classic"){
                    //update profit to cost of production
                    (*s)(j,14) -= (*s)(j,9) * (*s)(j,10);
                };
                
				//update cost of production
				(*s)(j,7) -= atc * (*s)(j,11);
				(*s)(j,8) -= (*s)(j,11);
			};
            
			//update profit to sales
			(*s)(j,14) += (*s)(j,11) * (*s)(j,12);
            
			atc = 0.0;
		};
        
        //        std::cout << *s << std::endl;
        
		//labor payment
		s->col(0) = s->col(0).array() -  (s->col(9).array() * s->col(10).array());
        
        //        std::cout << *s << std::endl;
        
		//dividends
		s->col(0) = s->col(0).array() - (s->col(14).array() > s_zero.array()).cast<double>().array() * s->col(14).array() * s->col(16).array();
        
//        std::cout << *s << std::endl;
        

		
        //zero out labor contract
		s->col(9) = s_zero;
		s->col(10) = s_zero;
        
        
        //dead or alive
		s->col(15) = (s->col(0).array() >= s_zero.array()).cast<double>().array() * s->col(15).array();
        
        //        std::cout << *s << std::endl;
        
		//stores realized goal
		c->col(i) = (s->col(15).array() > s_zero.array()).cast<double>().array() * s->col(14).array();

        
        //add value function estimation
        for (int j=0; j < s->rows(); j++){
            s_block  = s->block(j,0,1,s->cols());
            
            (*v)(j,i) = _V(s_block,ret, theta_V);
        };
        
//        std::cout << *c << std::endl;
        
//        std::cout << *v << std::endl;
        
		//zero out current profit
		s->col(14) = s_zero;
        
        
	};
    
    
};





void
F::opt_QL(){
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
    double max = Q(6,Q_s.at(0));
    
    for (auto j:Q_s){
        if (Q(6,j) > max){
            max = Q(6,j);
            Q_max.clear();
            Q_max.push_back(j);
            
        }else{
            if (Q(6,j) == max){
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
    theta_x.row(0) = Q.block(s_bin.size(),theta_max_i,4,1).transpose();
    
};

bool
F::_is_bin(int j, int i, std::vector<double>& s_bin){
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
F::wm_update_Q(){
    //updates Q matrix
    //reward
    double R_t1 = account->profit.at(account->profit.size()-2);
    
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
    double max = Q(6,Q_s.at(0));
    
    for (auto j:Q_s){
        if (Q(6,j) > max){
            max = Q(6,j);
            Q_max.clear();
            Q_max.push_back(j);
            
        }else{
            if (Q(6,j) == max){
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
    Q(6,theta_x_t_bin.back()) = Q(6,theta_x_t_bin.back())
                                + Q(7,theta_x_t_bin.back())
                                * (R_t1 + param->BETA*Q(6,theta_max_i) - Q(6,theta_x_t_bin.back()));
//    std::cout << Q << std::endl;
    
};



void
F::_set_QL(){
    //setup QL matrix
    //create matrix for the grid - small size here
	Q_grid = Eigen::MatrixXd::Zero(8,3);

    //money and stock of goods on hand
    Q_grid(0,0) = 0.0;
    Q_grid(0,1) = 100.0;
    Q_grid(0,2) = ((Q_grid(0,1) - Q_grid(0,0))/2);
    Q_grid(1,0) = 0.0;
    Q_grid(1,1) = w.NH;
    Q_grid(1,2) = ((Q_grid(1,1) - Q_grid(1,0))/2);
    
    create_decision_grid(Q_grid, 2);
    
    Q_grid(6,0) = 0.0;
    Q_grid(6,1) = 0.1;
    Q_grid(6,2) = 1.0;

    Q_grid(7,0) = param->QL_L; //speed of learning
    Q_grid(7,1) = Q_grid(7,0) + 0.1;
    Q_grid(7,2) = Q_grid(7,0) + 1.0;
    
    if (w.param->SIMULATION_MODE == "test"){
        Tools::mgrid_test(Q_grid, &Q);
        
    } else{
        
        Tools::mgrid(Q_grid, &Q);
    };
    
    create_bin_values(Q_grid, Q_bin);

    //update Q to push to hiring all of them
    for (int j = 0; j < Q.cols(); j++){
        if (Q(2,j) == Q_bin.at(2).back()){
            Q(6,j) = 0.05;
        };
    };
    
    
//    std::cout << Q << std::endl;
    
//    Tools::print_vector(Q_bin);
    
};


void
F::create_bin_values(Eigen::MatrixXd& _grid, std::vector<std::vector<double>>& _bin){
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
F::_place_state_QL(Eigen::MatrixXd& s){
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
    
    j = Q_bin.at(1).size()-1;
    FLAG_PLACED = false;
    while ((j>=0) && !FLAG_PLACED){
        
        if (s(0,2) >= Q_bin.at(1).at(j)){
            s_bin.push_back(j);
            FLAG_PLACED = true;
        };
        j -= 1;
    };
    
    return s_bin;
};



void
F::_set_mRE(){
    //
    mRE_grid = Eigen::MatrixXd::Zero(5,3);
    create_decision_grid(mRE_grid, 0);
    mRE_grid(4,0) = 1.0;
    mRE_grid(4,1) = 1.1;
    mRE_grid(4,2) = 1.0;
    
    if (w.param->SIMULATION_MODE == "test"){
        Tools::mgrid_test(mRE_grid, &mRE);
        
    } else{
        
        Tools::mgrid(mRE_grid, &mRE);
    };
    
    
    
    create_bin_values(mRE_grid, mRE_bin);
    
    //update mRE to push to hiring all of them
    for (int j = 0; j < mRE.cols(); j++){
        if (mRE(0,j) == mRE_bin.at(0).back()){
            mRE(4,j) = 1.05;
        };
    };
    
    
//    //storage for actual probabilities
//    mRE_grid(5,0) = 0.0;
//    mRE_grid(5,1) = 0.1;
//    mRE_grid(5,2) = 1.0;
    
};



void
F::opt_mRE(){
    //form matrix of probabilities
    double sum_prob = 0.0;
    
    for (int j=0; j<mRE.cols();j++){
        sum_prob += exp(mRE(4,j)/param->mRE_T);
    };
    
    std::vector<double> probs;
    
    for (int j=0; j<mRE.cols();j++){
        probs.push_back(exp(mRE(4,j)/param->mRE_T)/sum_prob);
    };

    
    //roll weighted dice
    boost::random::discrete_distribution<> dist(probs.begin(), probs.end());
    
    
    int theta_i = dist(w.rng);
    
    theta_x_t_bin.push_back(theta_i);
    
    //transform from index to decision
    theta_x.row(0) = mRE.block(0,theta_i,4,1).transpose();
    
};



void
F::wm_update_mRE(){
    //adjustment matrix
    Eigen::MatrixXd _E_re  = Eigen::MatrixXd::Zero(1, mRE.cols());
    
    //reward
    double R_t1 = account->profit.at(account->profit.size()-2);
    
    //adjustment
    _E_re = mRE.row(4) * (param->mRE_EPSILON/(mRE.cols()-1));
    
    //include reward
    _E_re(0,theta_x_t_bin.back()) = R_t1 * (1 - param->mRE_EPSILON);
    
    //update whole matrix
    mRE.row(4) = mRE.row(4) * (1 - param->mRE_PHI) + _E_re;
    
//    std::cout << mRE << std::endl;
};



void
F::_PS_receive_payment(MessagePSSendPayment* mes){};

bool
F::_PS_accept_PO(MessagePSSendPayment* mes){
	return true;
};


ContractBDt0*
F::_PS_contract(){
	return ass_asCBDt0.front();
};

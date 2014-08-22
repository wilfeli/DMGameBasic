/*
 * W.cc
 *
 *  Created on: Apr 13, 2014
 *      Author: wilfeli
 */


#include "CB.h"
#include "W.h"
#include "Market.h"
#include "F.h"
#include "H.h"
#include "Parameters.h"
#include "Message.h"

W::W(double seed):rng(seed), myrng(seed){
	//parameters
	param = new ParametersW();
	param->W_SEED = seed;

};

void
W::init(){

	//make cb
	new CB(*this);

	//make number of firms
	for (int i=0; i!=NFGC; ++i){
		//create bunch of firms
		new F(*this, i);

//		//if want to update wm - create new object and replace old
//		//create wm
//		ExpectationBackward* exp = new ExpectationBackward(1.0, 10, 0.5);
//
//		delete Fs.back()->wm["MasCHK"];
//		Fs.back()->wm["MasCHK"] = exp;

//		//if want to update production function parameters
//		Fs.back()->param->F_F_theta = param->F_F_F_theta;
        
        
        //if want to change accounting type
        Fs.back()->param->ACCOUNTING_TYPE = param->F_ACCOUNTING_TYPE;
        Fs.back()->account->ACCOUNTING_TYPE = Fs.back()->param->ACCOUNTING_TYPE;

		//if want to update wm length
		Fs.back()->param->WM_LENGTH = param->F_wm_LENGTH;

		//if want to update forecasting depth
		Fs.back()->param->T_MAX = param->F_T_MAX;

        //if want to update number of random draws
		Fs.back()->param->opt_CS_N = param->F_opt_CS_N;

        
		//if want to update grid type
		Fs.back()->param->GRID = param->F_GRID;
        
        //if want to change decision type
        Fs.back()->param->opt_TYPE = param->F_opt_TYPE;
        
        //update QL parameters
        Fs.back()->param->QL_L = param->F_QL_L;
        
        //update mRE parameters
        Fs.back()->param->mRE_PHI = param->F_mRE_PHI;
        Fs.back()->param->mRE_EPSILON = param->F_mRE_EPSILON;
        Fs.back()->param->mRE_T = param->F_mRE_T;

        //update internals to new parameters
        Fs.back()->param->init();
        Fs.back()->init();
        

        
		//create accounts for firms
		CBs->_create_asCBDt0(Fs.back());

	};


	//make number of Humans
	for (int i=0; i!=NH; ++i){
		//create bunch of Humans
		new H(*this, i);

//		//if want to update wm - create new object and replace old
//		//create wm
//		ExpectationBackward* exp = new ExpectationBackward(1.0, 10, 0.5);
//
//		delete Hs.back()->wm["MasCHK"];
//		Hs.back()->wm["MasCHK"] = exp;

//		//if want to update utility parameters
//		Hs.back()->param->GOAL_T_theta = param->H_GOAL_T_theta;

		//if want to update wm length
		Hs.back()->param->WM_LENGTH = param->H_wm_LENGTH;

		//if want to update forecasting depth
		Hs.back()->param->T_MAX = param->H_T_MAX;

        //if want to update number of random draws
		Hs.back()->param->opt_CS_N = param->H_opt_CS_N;

        
		//if want to update grid type
		Hs.back()->param->GRID = param->H_GRID;
        
        //if want to change decision type
        Hs.back()->param->opt_TYPE = param->H_opt_TYPE;

        //update QL parameters
        Hs.back()->param->QL_L = param->H_QL_L;
        
        //update mRE parameters
        Hs.back()->param->mRE_PHI = param->H_mRE_PHI;
        Hs.back()->param->mRE_EPSILON = param->H_mRE_EPSILON;
        Hs.back()->param->mRE_T = param->H_mRE_T;
        
        //update internals to new parameters
        Hs.back()->param->init();
        Hs.back()->init();

		//create accounts for Humans
		CBs->_create_asCBDt0(Hs.back());
	};


	//make market L
	new MarketL(*this);

	//make market C
	new MarketC(*this);

	//add firms to markets
	for (auto agent=Fs.begin();agent!=Fs.end();++agent){
		ml->AddBuyerL(*agent);
		mc->AddSellerC(*agent);
	}

	//add consumers to markets
	for (auto agent=Hs.begin();agent!=Hs.end();++agent){
		ml->AddSellerL(*agent);
		mc->AddBuyerC(*agent);
	}


};

void
W::step(){
	//update timing
	t++;
    
//    double snan = std::numeric_limits<double>::signaling_NaN();
    //check that no out of bounds errors
    
    if ((t>0.0) && (std::isnan(ml->market_price.at(t - 1)->p) || std::isnan(mc->market_price.at(t - 1)->p))){
        
        //mark all firms as bankrupt
		MessageBankruptcy* mes;
        for (auto a:Fs){
            if (a->status > 0.0){
                //calls bankruptcy
                mes = new MessageBankruptcy(a);
                ac_LS_bankruptcy(mes);
                delete mes;
            };
        };
        
        MessageBankruptcyH* mesH;
        //mark all humans as bankrupt
        for (auto a:Hs){
            if (a->status > 0.0){
                //calls bankruptcy
                mesH = new MessageBankruptcyH(a);
                ac_LS_bankruptcy(mesH);
                delete mesH;
            };
        };
        
    };
    

	//update if alive
	MessageStatus* mes_status = new MessageStatus(0.0);

	for (auto a:Fs){
		if (a->status > 0.0){
			a->ac_W_begin_step();
		}else{
			a->ac_W_begin_step(mes_status);
		};
	};

	

    
    
	for (auto a:Hs){
		if (a->status > 0.0){
			a->ac_W_begin_step();
		}else{
            a->ac_W_begin_step(mes_status);
        };
	};

	//update knowledge
    
    if (t > 0){
        for (auto a:Fs){
            if (a->status > 0.0){
                a->wm_update_k_s();
            };
        };


        for (auto a:Hs){
            if (a->status > 0.0){
                a->wm_update_k_s();
            };
        };
    };

	//make decisions
	MessageMakeDec* mes = new MessageMakeDec();


	for (auto a:Fs){
		if (a->status > 0.0){
			a->ac_W(mes);
		};
	};


	for (auto a:Hs){
		if (a->status > 0.0){
			a->ac_W(mes);
		};
	};

	delete mes;

	//market function
	ml->step();

	MessageF_F* mes_prod = new MessageF_F();
	//production happens
	for (auto a:Fs){
		if (a->status > 0.0){
			a->ac_W(mes_prod);
		};
	};

	delete mes_prod;

	//goods market function
	mc->step();


	for (auto a:Fs){
		if (a->status > 0.0){
			a->ac_W_end_step();
		};
	};


	for (auto a:Hs){
		if (a->status > 0.0){
			a->ac_W_end_step();
		}else{
            a->ac_W_end_step(mes_status);
        };
	};
    
    delete mes_status;



};

void
W::AddF(F* agent) { Fs.push_back(agent); };

void
W::AddCB(CB* agent) { CBs = agent; };


void
W::AddH(H* agent) { Hs.push_back(agent); };



void
W::AddMarketL(MarketL* agent) { ml = agent; };

void
W::AddMarketC(MarketC* agent) { mc = agent; };

void
W::ac_LS_bankruptcy(MessageBankruptcy* mes){
	//mark agent as dead

	mes->agent->status = 0;

	double q_dist = std::max(mes->agent->ass_asCBDt0.front()->q, 0.0);

	MessagePSSendPayment* inf;
	double pq = 0.0;

	//distribute money back to humans
	for (auto a:Hs){
		pq = q_dist/NH;
		inf = new MessagePSSendPayment(a, pq ,"Bankruptcy_b");
		inf->sender = mes->agent;

		mes->agent->ass_asCBDt0.front()->issuer->_PS_accept_PO(inf);

		delete inf;

	};

};

void
W::ac_LS_bankruptcy(MessageBankruptcyH* mes){
	//mark agent as dead
    
	mes->agent->status = 0;
};


void
W::print_step(){
	Eigen::MatrixXd data_all(1, 4);
	data_all = Eigen::MatrixXd::Zero(1, 4);

//	std::cout << data_all << std::endl;


	for (auto a:Fs){
        
        if (a->status > 0.0){
            //add all production
            data_all(0,0) += a->account->production_q.back();


            //add all sales from F side
            data_all(0,1) += a->account->sales_q.back();

            //add all profit
            data_all(0,2) += a->account->profit.back();
        };
	};

	for (auto a:Hs){
		//add all utility and average
		data_all(0,3) += a->account->g_t.back();

	};

	data_all(0,3) = data_all(0,3)/NH;

	std::cout << data_all << std::endl;

};

//W::~W(){
//    //call destructor of firms and consumers
//    
//};


/*
 * Accounting.cc
 *
 *  Created on: Apr 29, 2014
 *      Author: wilfeli
 */

#include <Eigen/Dense>

#include "Accounting.h"
#include "Contract.h"
#include "Tools.h"





void
AccountH::ac_W_end_step(Eigen::MatrixXd* s, double ut){
	//stores utility
	g_t.back() = ut;

	//stores consumption of a good - assume consume all
	g_GC_q_t.back() = (*s)(0,3);
	g_CHK_q_t.back() = (*s)(0,2);

};



void
AccountF::ac_W_end_step(Eigen::MatrixXd* s){
	//update cost according to production
    
    
    if (ACCOUNTING_TYPE == "modern"){
    
        if ((cost_q + production_q.back()) > 0.0){
            cost_p = (cost_q * cost_p + expences_wage.back())/(cost_q + production_q.back());
        } else{
            cost_p = 0.0;
        };

        cost_q += production_q.back();

        //update profit
        //add income
        profit.back() += sales_pq.back();
        //subtract expenses
        profit.back() -= cost_p * sales_q.back();


    //	Tools::print_vector(profit);


        //update costs according to sales
        cost_q -= sales_q.back();

        //handle case of zero current sales
        if (sales_q.back() == 0.0){
            profit.back() -= expences_wage.back();
            if (cost_q > 0.0){
                cost_p -= expences_wage.back() / cost_q;
            }else{
                cost_p = 0.0;
            };
        };
    };
    
    if (ACCOUNTING_TYPE == "classic"){
        
        if ((cost_q + production_q.back()) > 0.0){
            cost_p = (cost_q * cost_p + expences_wage.back())/(cost_q + production_q.back());
        } else{
            cost_p = 0.0;
        };
        
        cost_q += production_q.back();
        
        //update profit
        //add income
        profit.back() += sales_pq.back();
        //subtract expenses
        profit.back() -= expences_wage.back();
        
        
        //	Tools::print_vector(profit);
        
        
        //update costs according to sales
        cost_q -= sales_q.back();
        
    };
    


};



void
AccountH::ac_get_inf(MessageMarketLSellC* mes, ContractHK* c){
	//
	asCHK_p.back() += c->p;
	asCHK_q.back() += c->q;

};


void
AccountH::ac_get_inf(MessagePSSendPayment* mes){
	if (mes->type == "buy_asGC"){
		asGC_pq.back() += mes->q;
	};

	if (mes->type == "asCHK_b"){
		//is stored in previous period
		asCHK_p_t.back() += mes->q;
	};

	if (mes->type == "FI_div"){
		//is stored in previous period
			FI_div_t.back() += mes->q;
	};
};


void
AccountF::ac_get_inf(Message* mes, ContractHK* c){
	//contract length
	int c_length = c->t_end - c->t_begin + 1;

	//adds wage expenses to current expenses account proportional to contract length
	//other expenses are calculated from active contracts later
	expences_wage.back() += (c->q * c->p )/c_length;
	asCHK_q.back() += c->q;

	//


};


void
AccountF::ac_get_inf(MessageF_F* mes){
	production_q.back() += mes->q;
};


void
AccountF::ac_get_inf(MessageMarketCSellGoodC* mes){
	//

	sales_q.back() += mes->q;
	sales_p.back() = mes->p;
	sales_pq.back() += mes->q * mes->p;
};


void
AccountF::ac_W_initialize_step(){
	//pushes zeros for current time period
	expences_wage.push_back(0.0);
	sales_q.push_back(0.0);
	sales_p.push_back(0.0);
	sales_pq.push_back(0.0);
	asCHK_q.push_back(0.0);
	production_q.push_back(0.0);
	profit.push_back(0.0);

	//checks for active contracts


};

void
AccountF::ac_W_initialize_step(MessageStatus* mes){

	if (mes->status <= 0.0){
		//pushes nan for current time period
		double snan = std::numeric_limits<double>::signaling_NaN();
		expences_wage.push_back(snan);
		sales_q.push_back(snan);
		sales_p.push_back(snan);
		production_q.push_back(snan);
		profit.push_back(snan);
	};


};


void
AccountH::ac_W_initialize_step(MessageStatus* mes){
    
	if (mes->status <= 0.0){
        //pushes nan for current time period
		double snan = std::numeric_limits<double>::signaling_NaN();
        asCHK_p.push_back(snan);
        asCHK_q.push_back(snan);
        asGC_p.push_back(snan);
        asGC_pq.push_back(snan);
        asGC_q.push_back(snan);
        g_t.push_back(0.0);
        g_GC_q_t.push_back(0.0);
        g_CHK_q_t.push_back(0.0);
        FI_div_t.push_back(snan);
        asCHK_p_t.push_back(snan);
        asCBDt0_pq_t.push_back(snan);
        
    };
    
};

void
AccountH::ac_W_initialize_step(){

	asCHK_p.push_back(0.0);
	asCHK_q.push_back(0.0);
	asGC_p.push_back(0.0);
	asGC_pq.push_back(0.0);
	asGC_q.push_back(0.0);
	g_t.push_back(0.0);
	g_GC_q_t.push_back(0.0);
	g_CHK_q_t.push_back(0.0);
	FI_div_t.push_back(0.0);
	asCHK_p_t.push_back(0.0);
	asCBDt0_pq_t.push_back(0.0);
};




void
AccountH::ac_get_inf(MessageGoodC* mes){
	asGC_q.back() += mes->q;
};



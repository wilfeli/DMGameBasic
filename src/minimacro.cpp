//============================================================================
// Name        : minimacro.cpp
// Author      : wilfeli
// Version     :
// Copyright   : Your copyright notice
// Description : Minimacro model in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <Eigen/Dense>

#include "W.h"
#include "Parameters.h"
#include "H.h"
#include "F.h"
#include "Market.h"


using namespace std;



void
save_history(W& w, int N, std::string file_path){
	//saves history
	//parameters of the simulation

	std::string file_path_default;
	std::string file_name;

    
    
    if (file_path.empty()){
        file_path = file_path_default;
    };

	int i = 0;
	//read current number of iterations
	file_name = file_path + "HF_i.txt";
	std::ifstream in_file(file_name);

	std::string I_NUMBER;
    
    if (in_file){
        std::getline(in_file, I_NUMBER);
    }else{
        I_NUMBER = "0";
    };
	in_file.close();

	std::ofstream out_file(file_name);
	out_file << std::to_string(std::stoi(I_NUMBER) + 1) << std::endl;
    out_file.close();


	file_name = file_path + "HF_param_" + I_NUMBER + ".txt";
	out_file.open(file_name);

	if (out_file){
		//save number of H
		//save number of F
		out_file << "number of humans: "<< w.NH << std::endl;
		out_file << "number of firms: "<< w.NFGC << std::endl;
		out_file << "F wm length: " << w.param->F_wm_LENGTH << std::endl;
		out_file << "H wm length: " << w.param->H_wm_LENGTH << std::endl;
		out_file << "F type: " << w.Fs.back()->param->opt_TYPE << std::endl;
		out_file << "H type: " << w.Hs.back()->param->opt_TYPE << std::endl;
		out_file << "seed: " << w.param->W_SEED << std::endl;
		out_file << "F forecasting length: " << w.param->F_T_MAX << std::endl;
		out_file << "H forecasting length: " << w.param->H_T_MAX << std::endl;
		out_file << "F production function parameters: " << w.param->F_F_F_theta << std::endl;
		out_file << "H utility function parameters: " << w.param->H_GOAL_T_theta << std::endl;
		out_file << "F mRe Phi: " << w.param->F_mRE_PHI << std::endl;
		out_file << "H mRe Phi: " << w.param->H_mRE_PHI << std::endl;
		out_file << "F QL L: " << w.param->F_QL_L << std::endl;
		out_file << "H QL L: " << w.param->H_QL_L << std::endl;
		out_file << "F Grid: " << w.param->F_GRID << std::endl;
		out_file << "H Grid: " << w.param->H_GRID << std::endl;
        

		i = 0;

		for (auto a:w.Fs){
			if (a->status > 0.0){
				i++;
			};
		};

		out_file << "End number of F:" << i << std::endl;

	};

	out_file.close();
    
    Eigen::IOFormat CleanFmt(Eigen::StreamPrecision, 0, ", ", "\n");

	file_name = file_path + "HF_data_" + I_NUMBER + ".txt";
	out_file.open(file_name);
	//save prices of the market
	Eigen::MatrixXd data_all(N, 5);
	data_all = Eigen::MatrixXd::Zero(N, 5);
	//employment, sales on goods market, production, price on labor market,
	//price on goods market


	for (auto a:w.Hs){
		data_all.col(0) += Eigen::VectorXd::Map(a->account->asCHK_q.data(),(long)a->account->asCHK_q.size());
		data_all.col(1) += Eigen::VectorXd::Map(a->account->asGC_q.data(),(long)a->account->asGC_q.size());
	};

	for (auto a:w.Fs){
        if (a->status > 0.0){
            data_all.col(2) += Eigen::VectorXd::Map(a->account->production_q.data(),(long)a->account->production_q.size());
        };
	};


	for (int j = 0; j < N; ++j){
		data_all(j,3) = w.ml->market_price[j]->p;
		data_all(j,4) = w.mc->market_price[j]->p;

	};
	if (out_file){
		out_file << data_all.format(CleanFmt);
	};

//	cout << data_all;
	out_file.close();


	file_name = file_path + "HF_data_all_H_" + I_NUMBER + ".txt";
	Eigen::MatrixXd h_all(N, w.NH);
	i = 0;
	for (auto a:w.Hs){
		h_all.col(i) = Eigen::VectorXd::Map(a->account->g_t.data(),(long)a->account->g_t.size());
		i++;

	};

    out_file.open(file_name);
	if (out_file){
		out_file << h_all.format(CleanFmt);;
	};
	out_file.close();

	file_name = file_path + "HF_data_all_F_" + I_NUMBER + ".txt";
	Eigen::MatrixXd f_all(N, w.NFGC);
	i = 0;
	for (auto a:w.Fs){
		f_all.col(i) = Eigen::VectorXd::Map(a->account->profit.data(),(long)a->account->profit.size());
		i++;

	};
    out_file.open(file_name);
	if (out_file){
		out_file << f_all.format(CleanFmt);
	};
	out_file.close();


};


void
read_init(double& seed_, int& N_, std::map <std::string, std::string> &ini, std::string file_name){
	std::string file_path_default;
	std::string file_name_default;
    

	//open ini file
	file_name_default = "minimacro.ini";
    
    if (file_name.empty()){
        file_name = file_name_default;
    };
    
	std::ifstream in_file(file_name);
    
	std::string s, key, value;

    
    while (std::getline( in_file, s )){
        // Extract the key value
        std::string::size_type begin = s.find_first_not_of( " \f\t\v" );
        std::string::size_type end = s.find( '=', begin );
        key = s.substr( begin, end - begin );
        
        // (No leading or trailing whitespace allowed)
        key.erase( key.find_last_not_of( " \f\t\v" ) + 1 );
        
        // No blank keys allowed
        if (key.empty()) continue;
        
        // Extract the value (no leading or trailing whitespace allowed)
        begin = s.find_first_not_of( " \f\n\r\t\v", end + 1 );
        end   = s.find_last_not_of(  " \f\n\r\t\v" ) + 1;
        
        value = s.substr( begin, end - begin );
        ini[key] = value;

//        std::cout << s << " " << key << " " << value << std::endl;
	};
    in_file.close();

    //set ini values
    N_ = std::stoi(ini["N"]);
    seed_ = std::stoi(ini["seed"]);
    
};


int main(int argc, char *argv[]) {
    //path to ini file
    //argv[1] path to ini file
    //argv[2] path to save directory
    
    std::string file_ini_name;
    std::string file_path;
    //path to save directory
    if (argc > 2){
        file_path = argv[2];
    };
    
    if (argc > 1){
        file_ini_name = argv[1];
    };
    
	//seed
	double seed = 2013;

    //number of steps
    int N = 100;
    std::map <std::string, std::string> simulation_ini;

    read_init(seed, N, simulation_ini, file_ini_name);
    

	//create world
	W w(seed);


	w.NH = 10;
	w.NFGC = 3;
	w.NCB = 1;

//	w.param->F_F_F_theta << ?, ?;
//	w.param->H_GOAL_T_theta << ?, ?;
    
    w.param->SIMULATION_MODE = simulation_ini["SIMULATION_MODE"];
    w.param->F_T_MAX = std::stoi(simulation_ini["F_T_MAX"]);
    w.param->H_T_MAX = std::stoi(simulation_ini["H_T_MAX"]);
    w.param->F_opt_CS_N = std::stoi(simulation_ini["F_opt_CS_N"]);
    w.param->H_opt_CS_N = std::stoi(simulation_ini["H_opt_CS_N"]);
    
    std::string wm_length;
    wm_length = simulation_ini["F_wm_LENGTH"];
    if (wm_length == "inf"){
        w.param->F_wm_LENGTH = std::numeric_limits<double>::infinity();;
    }else{
        w.param->F_wm_LENGTH = std::stod(simulation_ini["F_wm_LENGTH"]);
    };
    wm_length = simulation_ini["H_wm_LENGTH"];
    if (wm_length == "inf"){
        w.param->H_wm_LENGTH = std::numeric_limits<double>::infinity();;
    }else{
        w.param->H_wm_LENGTH = std::stod(simulation_ini["H_wm_LENGTH"]);
    };
        

    w.param->F_GRID = simulation_ini["F_GRID"];
    w.param->H_GRID = simulation_ini["H_GRID"];
    
    //decision type
    w.param->F_opt_TYPE = simulation_ini["F_opt_TYPE"];
    w.param->H_opt_TYPE = simulation_ini["H_opt_TYPE"];

    //QL parameters
    w.param->F_QL_L = std::stod(simulation_ini["F_QL_L"]);
    w.param->H_QL_L = std::stod(simulation_ini["H_QL_L"]);
    
    //mRe parameters
//    w.param->F_mRE_EPSILON = std::stod(simulation_ini["F_mRE_EPSILON"]);
//    w.param->H_mRE_EPSILON = std::stod(simulation_ini["H_mRE_EPSILON"]);
//    w.param->F_mRE_T = std::stod(simulation_ini["F_mRE_T"]);
//    w.param->H_mRE_T = std::stod(simulation_ini["H_mRE_T"]);
    w.param->F_mRE_PHI = std::stod(simulation_ini["F_mRE_PHI"]);
    w.param->H_mRE_PHI = std::stod(simulation_ini["H_mRE_PHI"]);
    
    w.param->F_ACCOUNTING_TYPE = simulation_ini["F_ACCOUNTING_TYPE"];
    
    
    if (simulation_ini.find("NH") != simulation_ini.end()){
        w.NH = std::stoi(simulation_ini["NH"]);
    };

    if (simulation_ini.find("NFGC") != simulation_ini.end()){
        w.NFGC = std::stoi(simulation_ini["NFGC"]);
    };

    if (simulation_ini.find("NCB") != simulation_ini.end()){
        w.NCB = std::stoi(simulation_ini["NCB"]);
    };

    
    
	w.init();

	for (int i=0; i<N; ++i){
		w.step();
//        std::cout << i << " ";
//        w.print_step();
	};


//	w.print_step();

	save_history(w, N, file_path);


	return 0;
};







//============================================================================
// Name        : mots_3.cpp
// Author      : Christos Tsotskas
// Version     :
// Copyright   : Apache 2
// Description : MOTS in C++, Ansi-style
//============================================================================


//library that contains the tabu search function

//#include <cstring>
#include <sstream>
#include <fstream>

#include "mots.h"
#include "test_suite.h"
#include "objective_function_formulae.h"
#include "objective_function.h"
#include "verification_code.h"

#include <fstream>

#include <vector>

#include "easylogging++.h"

INITIALIZE_EASYLOGGINGPP



int main_example(int argc, char *argv[]){
	std::cout << "MOTS2 main" << std::endl;

	LOG(INFO) << "MOTS2 info log using default logger";


	std::string external_name;
	int flag=0;
	for(int f=0; f<argc; ++f){
		std::cout << f << ". " << argv[f] << std::endl;
		if(!strcmp(argv[f],"-r")  )
			flag=1;
		if(!strcmp(argv[f],"-c")  )
			external_name=argv[f+1];

	}


	benchmark_test_with_internal_configuration("LBM2_sensitivity_", "take_8",flag, 3 );

/*
	std::ostringstream integration_test_casename;
	test_report tr("t1", Point2(1,1), 1.0);
	std::ofstream overall_tests_report_file("all_tests.txt");


	if(overall_tests_report_file.is_open()){

		//			for(int i=1;i<6; ++i){
		int i=5;

		//				integration_test_casename << "case_" << i  ;
		integration_test_casename << external_name.c_str();

		//			tr=benchmark_test_with_external_configuration(integration_test_casename.str(), "ZDT2",flag);
		tr=benchmark_test_with_internal_configuration(integration_test_casename.str(), "ZDT2",flag, i*i*30 );

		std::cout 			<< "test_loop " << i << ": " << tr.show()  << std::endl;
		overall_tests_report_file 	<< "test_loop " << i << ": " << tr.show()  << std::endl;

		integration_test_casename.clear();
		integration_test_casename.str(std::string());
		//			}

		overall_tests_report_file.close();
	}else
		std::cout << "cannot open file!" << std::endl;
*/
	std::cout << "optimiser tested!" << std::endl;
	return 0;
}

/**
 * For ease of implementation, it is useful to link together all the configuration
settings that are related to the size of variables used. So, all of them should be
a function of the number of decision variables.
 */




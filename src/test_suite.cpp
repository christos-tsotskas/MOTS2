/*
 * test_suite.cpp
 *
 *  Created on: Jan 18, 2014
 *      Author: ctsotskas
 */


#include "mots.h"

#include "objective_function_formulae.h"
#include "objective_function.h"
#include "verification_code.h"

#include "test_suite.h"
#include "external_configuration_file.h"
#include "configuration_settings.h"

test_suite::test_suite(){
	//optimisation


	// zdt1 mots

}




void test_suite::test_zdt1(){

	//prepares a fresh run of mots2 for zdt1 configurations

	//feeds the information to mots2

	//runs mots2

	//collect results

}



void test_suite::test_all(){
	test_zdt1();

}



test_report benchmark_test_with_internal_configuration(std::string test_casename, const std::string function_name, int flag, int nVar){
	std::ofstream report_file("tests_report_file.txt", std::ios::app);
	time_t start_time; time (&start_time);

	std::cout << "test case: " << test_casename << " starts with internal configuration" << std::endl;


	general_configuration internal_configuration_old("no_filename",
			10, //1 - diversify
			5, //2 - intensify
			15, //3 - reduce
			0.00, //4 - SS
			0.5, //5 - SSRF
			1, //6 - save step
			3, //7 - sampling
			nVar, //8 - nVar
			2, //9 - nObj
			0, //10 loop limit
			3000, //11 evaluations limit
			0, //12 Improvements limit , number of consecutive improvements
			4, //13 - number of regions
			6, //14 - STM size
			"HV", // 15
			"full", //16
			-0.05, //17 - starting point
			200, //18
			300); //19

	general_configuration internal_configuration("configuration.txt");

	ObjFunction2 test_reference_point=ObjFunction2(2, 22.0);
	ObjFunction2 test_penalty_point=ObjFunction2(2,33333.0);

	Point2 test_lower_bound=Point2(nVar,0.0);
	test_lower_bound[0]=5.0;
	test_lower_bound[1]=1.0;
	test_lower_bound[2]=11.0;

	Point2 test_upper_bound=Point2(nVar,1.0);
	test_upper_bound[0]=11.0;
	test_upper_bound[1]=200.0;
	test_upper_bound[2]=29.0;


	Point2 test_starting_point=Point2(nVar,0.5);
	test_starting_point[0]=7.0;
	test_starting_point[1]=100.0;
	test_starting_point[2]=24.5;

	Point2 test_current_step=Point2(nVar,0.05);
	test_current_step[0]=0.1666666;
	test_current_step[1]=0.050251256;
	test_current_step[2]=0.055555556;

	Overall_Optimisation_Configuration_Settings myConf2(test_casename,
			internal_configuration,
			test_reference_point,
			test_penalty_point,
			test_lower_bound,
			test_upper_bound,
			test_starting_point,
			test_current_step);

	const unsigned int n_of_variables=myConf2.getExternalConfigurationFile().getVar();
	const unsigned int n_of_objectives=myConf2.getExternalConfigurationFile().getObj();

	objective_function_formulae obj_function(n_of_objectives);

	ObjectiveFunctionBasic<double> TS_ObjFunc(myConf2, obj_function);

	Container2 MTM(n_of_variables, n_of_objectives, "MTM","./memories");
	Container2 IM(n_of_variables, n_of_objectives, "IM","./memories");
	Container2 HISTORY(n_of_variables, n_of_objectives, "HISTORY","./memories");

	STM_Container2 STM(myConf2.getExternalConfigurationFile().getStmSize(), n_of_variables, "STM", "./memories");
	LTM_Container2Basic2<double> LTM( n_of_variables ,  myConf2.getExternalConfigurationFile().getRegions(), myConf2.get_lower_bound(), myConf2.get_upper_bound(),"LTM", "./memories");
	std::cout << "Memories done!" << std::endl;

	TabuSearch TS(myConf2, flag, TS_ObjFunc, MTM, IM, HISTORY, STM, LTM);

	double hyper_volume_indicator=TS.search2();

	time_t end; time (&end);
	double dif = difftime (end,start_time);
	std::cout << "end in " << dif<< "seconds" << std::endl;
	report_file << test_casename << "\t" << n_of_variables << "\t" << dif << " seconds " << __DATE__ << "\t" << __TIME__  << std::endl;
	report_file.close();

	return test_report(myConf2.getCaseName(), TS.getDatumPnt(), hyper_volume_indicator) ;
}

test_report benchmark_test_with_external_configuration(std::string test_casename, const std::string function_name, int flag){
	std::cout << "test case: " << test_casename << " starts with external configuration" << std::endl;

	general_configuration ext_conf_file("configuration_"+function_name+".txt");

	Overall_Optimisation_Configuration_Settings myConf2(test_casename,
			ext_conf_file,
			"reference_point_"+function_name+".txt",
			"penalty_point_"+function_name+".txt",
			"lower_bound_"+function_name+".txt",
			"upper_bound_"+function_name+".txt",
			"starting_point_"+function_name+".txt",
			"current_step_"+function_name+".txt");

	const unsigned int n_of_variables=myConf2.getExternalConfigurationFile().getVar();
	const unsigned int n_of_objectives=myConf2.getExternalConfigurationFile().getObj();

	objective_function_formulae obj_function(n_of_objectives);

	ObjectiveFunctionBasic<double> TS_ObjFunc(myConf2, obj_function);

	Container2 MTM(n_of_variables, n_of_objectives, "MTM","./memories");
	Container2 IM(n_of_variables, n_of_objectives, "IM","./memories");
	Container2 HISTORY(n_of_variables, n_of_objectives, "HISTORY","./memories");

	STM_Container2 STM(myConf2.getExternalConfigurationFile().getStmSize(), n_of_variables, "STM", "./memories");
	LTM_Container2Basic2<double> LTM( n_of_variables ,  myConf2.getExternalConfigurationFile().getRegions(), myConf2.get_lower_bound(), myConf2.get_upper_bound(),"LTM", "./memories");
	std::cout << "Memories done!" << std::endl;

	TabuSearch TS(myConf2, flag, TS_ObjFunc, MTM, IM, HISTORY, STM, LTM);

	double hyper_volume_indicator=TS.search2();

	return test_report(myConf2.getCaseName(), TS.getDatumPnt(), hyper_volume_indicator) ;
}


void dummy_test(objective_function_formulae& dummy_objective_function_formulae){
	std::cout << "dummy case" << std::endl;
	Point2 d1=simple_read_vector_from_file("dummy_ZDT6.txt");
	std::cout << "objectives: "<<	dummy_objective_function_formulae(d1) << std::endl;
}


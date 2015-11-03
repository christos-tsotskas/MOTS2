/*
 * test_suite.h
 *
 *  Created on: Jan 18, 2014
 *      Author: ctsotskas
 */

#ifndef TEST_SUITE_H_
#define TEST_SUITE_H_

#include <sstream>
#include <ctime>

#include "global_defines.h"


class test_report{
	std::string __case_name;
	Point2 __initial_point;
	double __HV;


public:
	test_report ( std::string case_name,
			Point2  initial_point,
			double HV ):
				__case_name( case_name),
				__initial_point(initial_point),
				__HV(HV){


	}

	std::string show(){
		std::ostringstream string_line_report;
		string_line_report << __case_name << "\t"<< __initial_point << "\tHV:" << __HV;
		return string_line_report.str();
	}

	~test_report(){

	}
};


class test_suite{
	void test_zdt1();

public:
	test_suite();
	void test_all();



};


test_report benchmark_test_with_internal_configuration(std::string mycasename, const std::string function_name, int flag, int nVar);

test_report benchmark_test_with_external_configuration(std::string test_casename, const std::string function_name, int flag);

void dummy_test(objective_function_formulae& dummy_objective_function_formulae);
#endif /* TEST_SUITE_H_ */

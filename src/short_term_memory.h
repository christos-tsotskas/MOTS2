/*
 * short_term_memory.h
 *
 *  Created on: Jan 19, 2014
 *      Author: ctsotskas
 */

#ifndef SHORT_TERM_MEMORY_H_
#define SHORT_TERM_MEMORY_H_

#include <iostream>
#include <cstdlib>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <deque>
#include <numeric>
#include <sys/stat.h> //file_exists
#include <cmath>
#include <cassert>
#include <ctime>
#include <fstream>
#include <istream>
//#include <string>
#include <cstring>

#include "global_defines.h"

class STM_Container2 : public std::list< Point2 >{
public:
	unsigned int STM_size;
	unsigned int local_nVar;
	std::string container_name;
	std::string __save_path;

	STM_Container2(const int &input_STM_size, int const n, std::string given_name, std::string save_path);
	void import(Point2 const &P);
	void showTabu();
	bool isTabu(const Point2 &P);
	void save_stm_container(std::string case_qualifier, int evaluations_counter);
	int load_stm_container(const char* load_path);
};

#endif /* SHORT_TERM_MEMORY_H_ */

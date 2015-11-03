/*
 * temp_container.h
 *
 *  Created on: Jan 19, 2014
 *      Author: ctsotskas
 */

#ifndef TEMP_CONTAINER_H_
#define TEMP_CONTAINER_H_

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

class temp_Container2 : public std::set< Point2 >  {

public:
	temp_Container2();
	Point2 selectRandom();
	void displayContent();
	int checkduplicate();
};


#endif /* TEMP_CONTAINER_H_ */

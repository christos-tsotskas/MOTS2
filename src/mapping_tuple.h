/*
 * mapping_tuple.h
 *
 *  Created on: Jan 19, 2014
 *      Author: ctsotskas
 */

#ifndef MAPPING_TUPLE_H_
#define MAPPING_TUPLE_H_

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

struct mapping_tuple{
public:
	const Point2 Pnt;
	const ObjFunction2 ObjF;
	mapping_tuple(Point2 P_arg, ObjFunction2 O_arg);
};



#endif /* MAPPING_TUPLE_H_ */

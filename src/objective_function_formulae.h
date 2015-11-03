/*
 * objective_function_formulae.h
 *
 *  Created on: Jan 20, 2014
 *      Author: ctsotskas
 */

#ifndef OBJECTIVE_FUNCTION_FORMULAE_H_
#define OBJECTIVE_FUNCTION_FORMULAE_H_

#include "global_defines.h"

class objective_function_formulae{
	unsigned int __nObj;
public:
	objective_function_formulae(unsigned int nOfObjectives);
	ObjFunction2 operator() (const Point2 &input);
};

#endif /* OBJECTIVE_FUNCTION_FORMULAE_H_ */

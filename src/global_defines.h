/*
 * global_defines.h
 *
 *  Created on: Jan 19, 2014
 *      Author: ctsotskas
 */

#ifndef GLOBAL_DEFINES_H_
#define GLOBAL_DEFINES_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>

//#define GPU_evaluation
#define CPU_evaluation


//Hooke and Jeeves serial and serial
//#define HJ_SERIAL
#define HJ_PARALLEL


/*Objectives define:
 *activates respective models ( variables, objectives, constraints and objective function)
 *Currently the following models are hard-coded
 *ZDT1, ZDT2, ZDT3, ZDT4,  ZDT5, ZDT6, DTLZ1, DTLZ2, DTLZ3, DTLZ4, DTLZ5, DTLZ6, DTLZ7
 *Besides the above, the following options are available:
 *1)The "xfoil" should ONLY be used with the ffd-xfoil module!
 *2)The "user_defined" is used so as to write a user defined model within the code and compile it again (in case source code is provided).
 *3)the "external" calls the external evaluator whenever  it is needed (at will) and reads the interface-files straight away,
 *4)the "external2" creates a file and expects the external evaluator to create a response-file at the same directory.
 */

/* ATTENTION!!!!!!!!
 * if user_defined, xfoil, external, external2 options are selected, then the ranges need to be specified
 */
//#define external
//#include "interface.h" //for my framework

#define CONTROL_TEST_FUNCTIONS

/*
 *  the following files will be created (the ./monitor_data directory should exist in in advance)
./monitor_data/diversify.out
./monitor_data/evals.out
./monitor_data/hypervolume.out
./monitor_data/i_local.out
./monitor_data/im_size.out
./monitor_data/intensify.out
./monitor_data/step_size.out
./monitor_data/reduce.out
./monitor_data/basePoint.out
./monitor_data/update_memories.out
 */


/*
 * creates under the root directory the following files
 * optima_plot.txt
 * optima_current.txt
 * tabu_plot.txt
 */
//////#define main_report


//code operation functionalities
//#define DEBUG
//#define REPORT
//#define REPORT_CONTAINERS
//#define DEBUG_EVALUATING
#define QUICK_LOG
#define FULL_LOG

//STAND ALONE OR AS-EXTERNAL-OPTIMISER
//#define GATAC  //activate in order to work with external files/modules


////////////////OBJECTIVE FUNCTION OPTIONS


#define EPSILON 1e-6
#define MAGNIFIER 100000000

#define PI 3.141592653589793238462643383279;

//#define xfoil

#ifndef TYPE_DEF_POINT
#define TYPE_DEF_POINT

//defines the decision variables point
typedef std::vector<double> Point2;

//defines the objectives point
typedef std::vector<double> ObjFunction2;

typedef std::pair<Point2,ObjFunction2 > entry;
#endif //TYPE_DEF_POINT

template <typename T> std::ostream& operator<<(std::ostream &os, const std::vector<T>& v)	{
	for (std::size_t i = 0; i < v.size(); i++)
	{
		os << v[i] << " ";
	}
	return os;
}

template<typename T >
double euclidean_dist(const std::vector<T> &current_point , const std::vector<T> &reference_point){
	double output=0.0;

	for (unsigned int i=0 ; i< current_point.size() ; ++i)
		output+=pow( current_point[i] -reference_point[i] ,2);

	return pow(output,0.5);
}

template<typename T >
void copy_contents(std::vector<T>  &P, std::vector<T> const &A){
	for (std::size_t i = 0; i < P.size(); i++)
		P[i] = A[i];
}

template<typename T>
int compare_error ( std::vector<T> const &v1,  std::vector<T> const &v2){
	/*
	 * resolve whether two vectors are equal,
	 * return 1 if they are numerically equal
	 * 0, otherwise
	 */


	double const comp_epsilon=1e-4;
	unsigned int equality_counter=0;

	if (v1.size()!=v2.size()){
		std::cout << "compare_error, impcompatible size! (compare_error)" << std::endl;

		exit(-33);
	}else{


		for (unsigned int i = 0; i < v1.size(); ++i)
			if( std::abs(v1[i] - v2[i])<comp_epsilon  )
				++equality_counter;

		if (equality_counter==v1.size())
			return 1;



		//for (unsigned int i = 0; i < v1.size(); ++i)
		//	nv[i] = floor(mangifier*v1[i] - mangifier*v2[i])/mangifier;

	}

	return 0;
}


template <typename T>
T RandomNumber2(T min, T max){
	T r = rand() /  ((T)RAND_MAX );
	return (T)(min + r * (max - min));
}




#endif /* GLOBAL_DEFINES_H_ */

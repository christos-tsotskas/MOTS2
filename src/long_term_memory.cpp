/*
 * long_term_memory.cpp
 *
 *  Created on: Jan 19, 2014
 *      Author: ctsotskas
 */

#include <cmath>

#include "long_term_memory.h"

#ifndef OPERATORS
#define	OPERATORS

//todo the following bits should be replaced and refactored as a template within the header file
void correct_accuracy(Point2 &P){
	/*for ( Point::iterator it=P.begin(); it!=P.end(); ++it)
if( fabs(*it)<epsilon )
	 *it=0;*/

	for ( unsigned int i=0; i < P.size(); ++i)
		if(  std::abs( P[i] ) < EPSILON )
			P[i]=0.0;
}

#define round(r) ( (r) > 0.0) ? floor( (r) + 0.5) : ceil( (r) - 0.5);

Point2 operator- (const Point2 &v1, const Point2 &v2){
	Point2 nv(v1.size(),0.0);
//	initialise_vector(nv, 0);
	//double result=0.0;

	if (v1.size()!=v2.size()){
		std::cout << "error, incompatible size! (operator -)" << std::endl;
		std::exit(-30);
	}else{

		//calculate the absolute error
		//		for (unsigned int i = 0; i < v1.size(); ++i){
		//			result= (double) v1[i] - v2[i];
		//			if(  abs(result)<  EPSILON )
		//				result=0.0;
		//
		//			nv[i] = result;
		//
		//		}

		//2os tropos
		for (unsigned int i = 0; i < v1.size(); ++i){
			nv[i] = floor(MAGNIFIER*v1[i] - MAGNIFIER*v2[i])/MAGNIFIER ;
			if( std::abs(nv[i])<EPSILON )
				nv[i]=0.0;
		}

		//3os tropos
		//		for (unsigned int i = 0; i < v1.size(); ++i){
		//			nv[i]= (double) ( v1[i] - v2[i] );
		//			(nv[i] > 0.0) ? floor(nv[i] + 0.5) : ceil(nv[i] - 0.5);
		//
		//			if( abs(nv[i])<EPSILON )
		//				nv[i]=0.0;
		//		}

	}

	return nv;
}

Point2 operator+ (const Point2 &v1, const Point2 &v2){
	Point2 nv(v1.size(),0.0);
//	initialise_vector(nv, 0);
	//double result=0.0;

	if (v1.size()!=v2.size()){
		std::cout << "error, incompatible size! (operator +)" << std::endl;
		exit(-30);
	}else{

		//calculate the absolute error
		//		for (unsigned int i = 0; i < v1.size(); ++i){
		//			result= (double) v1[i] + v2[i];
		//			if( abs(result)<  EPSILON )
		//				result=0.0;
		//
		//			nv[i] = result;
		//
		//		}

		//2os tropos
		for (unsigned int i = 0; i < v1.size(); ++i){
			nv[i] = floor(MAGNIFIER*v1[i] + MAGNIFIER*v2[i])/MAGNIFIER ;
			if( std::abs(nv[i])<EPSILON )
				nv[i]=0.0;
		}

		//3os tropos
		//		for (unsigned int i = 0; i < v1.size(); ++i){
		//			nv[i]= (double) ( v1[i] + v2[i] );
		//			(nv[i] > 0.0) ? floor(nv[i] + 0.5) : ceil(nv[i] - 0.5);
		//
		//			if( abs(nv[i])<EPSILON )
		//				nv[i]=0.0;
		//		}

	}

	return nv;
}

Point2 operator/ (const Point2 &v1, const Point2 &v2){ //only needed for local regions in LTM
	Point2 nv(v1.size(),0.0);
//	initialise_vector(nv, 0);
	double result=0.0;

	if (v1.size()!=v2.size()){
		std::cout << "error, incompatible size! (operator /)" << std::endl;
		exit(-30);
	}else{

		//calculate the absolute error
		for (unsigned int i = 0; i < v1.size(); ++i){
			result= (double) v1[i] / v2[i];
			if( std::abs(result)<  EPSILON )
				result=0.0;

			nv[i] = result;

		}


		//for (unsigned int i = 0; i < v1.size(); ++i)
		//	nv[i] = floor(mangifier*v1[i] + mangifier*v2[i])/mangifier;

	}

	return nv;
}

Point2 operator/ (const Point2 &v1, const int N){
	Point2 nv(v1.size(),0.0);
//	initialise_vector(nv, 0);
	double result=0.0;

	if (v1.size()!=0){
		std::cout << "error, zero size vector! (operator /)" << std::endl;
		exit(-35);
	}else{

		//calculate the absolute error
		for (unsigned int i = 0; i < v1.size(); ++i){
			result= (double) v1[i] / N;
			if( std::abs(result)<  EPSILON )
				result=0.0;

			nv[i] = result;

		}


		//for (unsigned int i = 0; i < v1.size(); ++i)
		//	nv[i] = floor(mangifier*v1[i] + mangifier*v2[i])/mangifier;

	}

	return nv;
}




bool operator==(const Point2 &P1, const Point2 &P2) {
	// check lengths of lists
	if (P1.size() != P2.size()){
		std::cout << " comparison ==, incompatible size!" << std::endl;
		exit(-31);
	}
	for (unsigned int k=0; k<P1.size(); ++k)
		if (  fabs(P1[k]-P2[k]) > EPSILON ) return false;
	//if (  P1[k]!=P2[k] ) return false;

	// check contents of lists
	/*for (unsigned int k=0; k<P1.size(); ++k) {
if (  fabs(P1[k]-P2[k])>epsilon ) return false;
}*/
	return true;
}


bool operator>(const Point2 &P1, const Point2 &P2) {
	// check lengths of lists
	if (P1.size() != P2.size()){
		std::cout << " comparison > , incompatible size! (operator>)" << std::endl;
		exit(-32);
	}
	// check contents of lists
	/*for (unsigned int k=0; k<P1.size(); ++k) {
if (  fabs(P1[k]-P2[k])>epsilon ) return false;
}*/

	for(unsigned int i=0; i<P1.size(); ++i)
		if( fabs(P1[i]-P2[i]) < EPSILON  )
			return false;

	return true;
}


bool operator<(const Point2 &P1, const Point2 &P2) {
	// check lengths of lists
	if (P1.size() != P2.size()){
		std::cout << " comparison < , incompatible size! (operator<)" << std::endl;
		exit(-32);
	}
	// check contents of lists
	/*for (unsigned int k=0; k<P1.size(); ++k) {
if (  fabs(P1[k]-P2[k])>epsilon ) return false;
}*/

	for(unsigned int i=0; i<P1.size(); ++i)
		if(P1[i]>P2[i])
			return false;

	return true;
}



template<typename T >
std::vector<T> vector_power(std::vector<T> &P, const double &n){
	std::vector<T> output(P.size(),0.0);
	for (unsigned int i=0 ; i< P.size() ; ++i)
		output[i]=pow( P[i],n);

	return output;
}
template<typename T>
std::vector<T> vector_power(const std::vector<T> &P, const double &n){
	std::vector<T> output(P.size(),0.0);
	for (unsigned int i=0 ; i< P.size() ; ++i)
		output[i]=pow( P[i],n);

	return output;
}



int FileExists(char const *strFilename) { //// ORIGINAL was BOOL!!!!!
	struct stat stFileInfo;
	////bool blnReturn;
	int intStat;

	// Attempt to get the file attributes
	intStat = stat(strFilename,&stFileInfo);
	if(intStat == 0) {
		// We were able to get the file attributes
		// so the file obviously exists.
		////	blnReturn = true; //// ORIGINAL
		return 1;
	} else {
		// We were not able to get the file attributes.
		// This may mean that we don't have permission to
		// access the folder which contains this file. If you
		// need to do that level of checking, lookup the
		// return values of stat which will give you
		// more details on why stat failed.
		////	blnReturn = false; //// ORIGINAL
		return 0;
	}

	//// return(blnReturn); //// ORIGINAL
}


#endif

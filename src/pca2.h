/*
 * pca2.h
 *
 *  Created on: 26 Jan 2013
 *      Author: Christos
 */

#ifndef PCA2_H_
#define PCA2_H_

#include <iostream>
#include <vector>
//#include <fstream>
//#include <cmath>
//#include <cstdlib>

#define SIGN(a, b) ( (b) < 0 ? -fabs(a) : fabs(a) )
#define EPS 0.005

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)





void erhand(const std::string err_msg);

float pythag(float a,float b);

template <typename T>
T rescale(T min, T max, T ratio);

template <typename T>
void caclucate_mean(const int n_points, const int dimensionality, const std::vector< std::vector<T> > &data, std::vector<T> &mean);

template <typename T>
void calc_stddev(const int n_points, const int dimensionality, const std::vector< std::vector<T> > &data,  std::vector<T> &mean,  std::vector<T> &stddev);

template <typename T>
void triangular_decomposition(std::vector< std::vector<T> > &a,int n, std::vector<T> &d, std::vector<T> &e);

template <typename T>
void tridiagonal_QLi(std::vector<T> &d, std::vector<T> &e, int n, std::vector< std::vector<T> > &z);

template <typename T>
void PCA(std::vector< std::vector<T> > &data, const int n, const int m, std::vector<T> &importance_vector);

#include "pca2.hpp"


//template <typename T>
//void PCA2(T &data);

#endif /* PCA2_H_ */

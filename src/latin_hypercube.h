/*
 * latin_hypercube.h
 *
 *  Created on: 27 Jan 2013
 *      Author: Christos
 */

#ifndef LATIN_HYPERCUBE_H_
#define LATIN_HYPERCUBE_H_

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <fstream>
# include <cstring>
# include <sstream>
# include <ctime>
# include <vector>
#include <cassert>
#include <map>
#include <deque>

using namespace std;

int get_seed ( );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_uniform ( int ilo, int ihi, int *seed );
void latin_random (  int dim_num, int point_num, int *seed, double x[], const vector<double> &lower_bound, const vector<double> &upper_bound  );
int *perm_uniform ( int n, int base, int *seed );
float r4_abs ( float x );
int r4_nint ( float x );
double r8_uniform_01 ( int *seed );
void r8mat_write (  string output_filename, int m, int n, double table[]);
void r8mat_write2 (  string output_filename, int m, int n, double table[], vector<double> &lower_bound, vector<double> &upper_bound);
void timestamp ( );


#endif /* LATIN_HYPERCUBE_H_ */

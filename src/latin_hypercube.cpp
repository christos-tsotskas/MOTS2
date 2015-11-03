

#include "latin_hypercube.h"
#include "pca2.h"


using namespace std;

//defines the variables vector, type and length
typedef std::vector<double> Point2;

//defines the objectives vector, type and length
typedef std::vector<double> ObjFunction2;

typedef std::pair<Point2,ObjFunction2> tuple;



//****************************************************************************80
template <typename T>
std::ostream& operator<<(std::ostream &os, const std::vector<T>& v)	{
	for (std::size_t i = 0; i < v.size(); i++)
	{
		os << v[i] << " ";
	}
	return os;
}

template<typename T>
vector<T> eval_ZDT4(vector<T> &input)
{
	const double epsilon=1e-4;
	const double pi=3.141592653589793238462643383279;

	if( input.size() != 10){
		cout << " ZDT4 input has" << input.size() << "variables, instead of 10" << endl;
		exit(-100);
	}
	int i = 0;
	int n = input.size();
	double f1 = 0;
	double g = 0;
	double h = 0;
	vector<T> output(2,0.0);

	assert(n == 10); //variable vector length
	assert(output.size() == 2);

	f1 = input[0];

	for (i = 1; i < n; i++)
	{
		double x = input[i];
		g += x * x - 10 * cos(4 * pi * x);
	}
	g = 1 + 10 * (n - 1) + g;
	h = 1 - sqrt(f1 / g);

	output[0]=f1;
	output[1] = g * h;

	return output;
}

int dominates(const ObjFunction2 &O1, const ObjFunction2 &O2){ //container
	//the input arguments are the entries "as-is" on the map structure
	// MAKE SURE THAT BOTH ENTRIES ARE FULLY FILLED IN!!!!
	//dexetai ws argument olo to entry apo to antistoixo map!

	/*
	 * A point x1 , with an objective function vector F1, is said to dominate point x2,
	 * with an objective function vector F2, if no component of F1 is greater than its corresponding
	 * component in F2, and at least one component is smaller.
	 * <->
	 * F1's components are smaller than F2's
	 *
	 */

	//dominance counters
	unsigned int nGreater=0;
	unsigned int nLess=0;

	/*the names below correspond to map's iterator clause
	 *namely, first is the key
	 *		second is the value of the corresponding key
	 */



	for (unsigned int i=0 ; i< O1.size() ; i++) {
		//"gurnaei" tis times tou second,
		//plh8os epanalhpsewn oso megalo einai to objective function vector

		if ( (O1[i]-O2[i]) > 0 ) // xreiazetai kai ena shmeio
			nGreater++ ;
		else
			if ( (O1[i]-O2[i]) < 0 )
				nLess++ ;
		// (int) O1[i]-O2[i] . this is to be compared with size of the objective vector. if equals, then return 1

	}
	if (nGreater > 0 && nLess==0)
		return -1 ;
	else if (nLess > 0 && nGreater==0 )
		return 1 ;
	else
		return 0 ;

}

int removeDominatedPoints(map<Point2,ObjFunction2> &memory){
	int number_of_points_prior_removal=memory.size();
	map<Point2,ObjFunction2>::iterator itA, itB;
	std::deque<Point2> toBeRemoved;
	for(itA=memory.begin() ; itA!=memory.end() ; ++itA)
		for(itB=memory.begin() ; itB!=memory.end() ; ++itB)
			if( !(itA->first == itB->first) && dominates(    itA->second   ,     itB->second )==1 )
				toBeRemoved.push_back(itB->first);
	// remove point from container

	//updates in memory cannot be applied on the fly; that's why I log the points to be removed
	for(std::deque<Point2>::iterator itC=toBeRemoved.begin()  ; itC!=toBeRemoved.end() ; ++itC )
		memory.erase(*itC);

	return number_of_points_prior_removal-toBeRemoved.size();
}

int addIfNotDominated(map<Point2,ObjFunction2> &memory, const Point2 &V1,const ObjFunction2 &O1){
#ifdef DEBUG
	if(V1.empty()){
		cout << "variable 1 is empty!" << endl;
		exit(-111);
	}
	if(O1.empty()){
		cout << "objective 1 is empty!" << endl;
		exit(-222);
	}
#endif
	for (map<Point2,ObjFunction2>::iterator it=memory.begin() ; it!= memory.end() ; ++it)
		if( dominates(  it->second  ,  O1 )==1)  //first is the existing point
			return 0;
	memory.insert( tuple(V1, O1) );
	return 1;
}

void save_map(const string name, const map<Point2,ObjFunction2> &memory ){
	ofstream output(name.c_str());
	for(map<Point2,ObjFunction2>::const_iterator it=memory.begin(); it!=memory.end(); ++it )
		output << it->first << endl;
	output.close();
}

void save_full_map(const string name, const map<Point2,ObjFunction2> &memory ){
	ofstream output(name.c_str());
	for(map<Point2,ObjFunction2>::const_iterator it=memory.begin(); it!=memory.end(); ++it )
		output << it->first << " " << it->second << endl;
	output.close();
}

Point2 find_extreme(const map<Point2,ObjFunction2> &memory, int objective){
	map<Point2,ObjFunction2>::const_iterator it=memory.begin();

	Point2 proposed_vector=it->first;
	double min_objective=it->second[objective];
	++it;

	for(  ; it!=memory.end(); ++it ){
		if(it->second[objective]<min_objective ){
			min_objective=it->second[objective];
			proposed_vector=it->first;

		}
	}
	return proposed_vector;
}


int main4 ( int argc, char *argv[] ){

	const int n=2000; //number of points
	const int m=10; //dimensionality
	int seed=get_seed ( );

	double *lhs_data=new double[m*n]; //data for LHS

	map<Point2,ObjFunction2> initial_memory;
	Point2 penalty(2, 100000);
	Point2 temp_design(m,-1);

	vector< vector<double> > data3(n, vector<double>(m,-1) );//data for PCA
	vector<double> importance_vector(m,-1);

	//set variables' range
	vector<double> upper_bound(m,0.0), lower_bound(m,0.0);
	upper_bound[0]=1.0;
	lower_bound[0]=0.0;

	for(int i=1; i<m; i++){
		upper_bound[i]=5.0;
		lower_bound[i]=-5.0;
	}



	timestamp ( );


	cout << "  Spatial dimension M = " << m << "\n";
	cout << "  Number of points N = " << n << "\n";
	cout << "  The seed is = " << seed << "\n";


	latin_random ( m, n, &seed, lhs_data, lower_bound, upper_bound  );


	//prepare list for evaluation
	for(int i=0; i<n; ++i){
		for(int j=0; j<m; ++j){
			temp_design[j]=lhs_data[i*m+j];
		}
		addIfNotDominated(initial_memory, temp_design, eval_ZDT4(temp_design) );
	}

	save_full_map("mixed_solutions.txt", initial_memory);

	removeDominatedPoints(initial_memory);

	save_full_map("final_solutions.txt", initial_memory);

	cout << "suggestions:" << endl;
	temp_design=find_extreme(initial_memory, 0);
	cout << "A:"<< temp_design << " with " << initial_memory[temp_design] << endl;
	temp_design=find_extreme(initial_memory, 1);
	cout << "B:" << temp_design << " with " << initial_memory[temp_design] << endl;

	//prepare map data for PCA
	int i=0;
	for(map<Point2,ObjFunction2>::iterator it=initial_memory.begin(); it!=initial_memory.end(); ++it ){
		for(int j=0; j<m ;++j)
			data3[i][j]=it->first[j];
		++i;
	}


	//prepare LHS data for PCA
	//	cout << "copy to data" << endl;
	//	for(int i=0; i<n; ++i){
	//		for(int j=0; j<m; ++j){
	//			data3[i][j]=lhs_data[i*m+j];
	//			//			cout << data3[i][j]<< " ";
	//		}
	//		cout << endl;
	//	}

	//read from file
	//	double buffer;
	//	//read input matrix
	//	ifstream inputfile("MTM_xdat.txt");
	//
	//	for(int i=0; i<n; ++i){
	//		for(int j=0; j<m; ++j){
	//			inputfile>>buffer;
	//			data3[i][j]=buffer;
	//			cout << data3[i][j] << " ";
	//		}
	//		cout << endl;
	//	}
	//	inputfile.close();


	PCA(data3,n, m, importance_vector);

	cout << "Report:" << endl;
	for(int j=0; j<m; ++j){
		cout.precision(7.7);
		cout <<  importance_vector[j] << " ";
	}


	r8mat_write ( "data_set.txt", m, n, lhs_data);

	cout << endl <<"  Normal end of execution." << endl;

	timestamp ( );

	delete [] lhs_data;
	return 0;
}
//****************************************************************************80

int get_seed ( )

//****************************************************************************80
//
//  Purpose:
//
//    GET_SEED returns a random seed for the random number generator.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int GET_SEED, a random seed value.
//
{
# define I_MAX 2147483647
	time_t clock;
	//	int i;
	int ihour;
	int imin;
	int isec;
	int seed;
	struct tm *lt;
	time_t tloc;
	//
	//  If the internal seed is 0, generate a value based on the time.
	//
	clock = time ( &tloc );
	lt = localtime ( &clock );
	//
	//  Hours is 1, 2, ..., 12.
	//
	ihour = lt->tm_hour;

	if ( 12 < ihour )
	{
		ihour = ihour - 12;
	}
	//
	//  Move Hours to 0, 1, ..., 11
	//
	ihour = ihour - 1;

	imin = lt->tm_min;

	isec = lt->tm_sec;

	seed = isec + 60 * ( imin + 60 * ihour );
	//
	//  We want values in [1,43200], not [0,43199].
	//
	seed = seed + 1;
	//
	//  Remap SEED from [1,43200] to [1,IMAX].
	//
	seed = ( int )
    																( ( ( double ) seed )
    																		* ( ( double ) I_MAX ) / ( 60.0 * 60.0 * 12.0 ) );
	//
	//  Never use a seed of 0.
	//
	if ( seed == 0 )
	{
		seed = 1;
	}

	return seed;
# undef I_MAX
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
	if ( i2 < i1 )
	{
		return i1;
	}
	else
	{
		return i2;
	}

}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
	if ( i1 < i2 )
	{
		return i1;
	}
	else
	{
		return i2;
	}

}
//****************************************************************************80

int i4_uniform ( int a, int b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM returns a scaled pseudorandom I4.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
	int k;
	float r;
	int value;

	if ( *seed == 0 )
	{
		cerr << "\n";
		cerr << "I4_UNIFORM - Fatal error!\n";
		cerr << "  Input value of SEED = 0.\n";
		exit ( 1 );
	}

	k = *seed / 127773;

	*seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

	if ( *seed < 0 )
	{
		*seed = *seed + 2147483647;
	}

	r = ( float ) ( *seed ) * 4.656612875E-10;
	//
	//  Scale R to lie between A-0.5 and B+0.5.
	//
	r = ( 1.0 - r ) * ( ( float ) ( i4_min ( a, b ) ) - 0.5 )
    																+         r   * ( ( float ) ( i4_max ( a, b ) ) + 0.5 );
	//
	//  Use rounding to convert R to an integer between A and B.
	//
	value = r4_nint ( r );

	value = i4_max ( value, i4_min ( a, b ) );
	value = i4_min ( value, i4_max ( a, b ) );

	return value;
}
//****************************************************************************80

void latin_random ( int dim_num, int point_num, int *seed, double x[],const vector<double> &lower_bound,const vector<double> &upper_bound )

//****************************************************************************80
//
//  Purpose:
//
//    LATIN_RANDOM returns points in a Latin Random square.
//
//  Discussion:
//
//    In each spatial dimension, there will be exactly one
//    point whose coordinate value lies between consecutive
//    values in the list:
//
//      ( 0, 1, 2, ..., point_num ) / point_num
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input/output, int *SEED, a seed for UNIFORM.
//
//    Output, double X[DIM_NUM,POINT_NUM], the points.
//
{
	int base = 0;
	int i;
	int j;
	int k;
	int *perm;
	double r;
	//
	//  For spatial dimension I,
	//    pick a random permutation of 1 to POINT_NUM,
	//    force the corresponding I-th components of X to lie in the
	//    interval ( PERM[J]-1, PERM[J] ) / POINT_NUM.
	//
	k = 0;
	for ( i = 0; i < dim_num; i++ )
	{
		perm = perm_uniform ( point_num, base, seed );

		for ( j = 0; j < point_num; j++ )
		{
			r = r8_uniform_01 ( seed );
			x[k] = ( ( ( double ) perm[j] ) + r ) / ( ( double ) point_num );
			k = k + 1;
		}
		delete [] perm;
	}


	for ( j = 0; j < point_num; j++ )
	{
		for ( i = 0; i < dim_num; i++ )
		{
			x[i+j*dim_num]=rescale( lower_bound[i], upper_bound[i],x[i+j*dim_num]) ;
		}

	}





	return;
}
//****************************************************************************80

int *perm_uniform ( int n, int base, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_UNIFORM selects a random permutation of N objects.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the number of objects to be permuted.
//
//    Input, int BASE, is 0 for a 0-based permutation and 1 for
//    a 1-based permutation.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int PERM_UNIFORM[N], a permutation of (BASE, BASE+1, ..., BASE+N-1).
//
{
	int i;
	int j;
	int k;
	int *p;

	p = new int[n];

	for ( i = 0; i < n; i++ )
	{
		p[i] = i + base;
	}

	for ( i = 0; i < n; i++ )
	{
		j = i4_uniform ( i, n - 1, seed );
		k    = p[i];
		p[i] = p[j];
		p[j] = k;
	}

	return p;
}
//****************************************************************************80

float r4_abs ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ABS returns the absolute value of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the quantity whose absolute value is desired.
//
//    Output, float R4_ABS, the absolute value of X.
//
{
	float value;

	if ( 0.0 <= x )
	{
		value = x;
	}
	else
	{
		value = -x;
	}
	return value;
}
//****************************************************************************80

int r4_nint ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NINT returns the nearest integer to an R4.
//
//  Example:
//
//        X         R4_NINT
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the value.
//
//    Output, int R4_NINT, the nearest integer to X.
//
{
	int value;

	if ( x < 0.0 )
	{
		value = - ( int ) ( r4_abs ( x ) + 0.5 );
	}
	else
	{
		value =   ( int ) ( r4_abs ( x ) + 0.5 );
	}

	return value;
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      r8_uniform_01 = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate,
//    strictly between 0 and 1.
//
{
	int k;
	double r;

	if ( *seed == 0 )
	{
		cerr << "\n";
		cerr << "R8_UNIFORM_01 - Fatal error!\n";
		cerr << "  Input value of SEED = 0.\n";
		exit ( 1 );
	}

	k = *seed / 127773;

	*seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

	if ( *seed < 0 )
	{
		*seed = *seed + 2147483647;
	}
	//
	//  Although SEED can be represented exactly as a 32 bit integer,
	//  it generally cannot be represented exactly as a 32 bit real number!
	//
	r = ( double ) ( *seed ) * 4.656612875E-10;

	return r;
}
//****************************************************************************80

void r8mat_write ( string output_filename, int m, int n, double table[])

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the table data.
//
{
	int i;
	int j;
	ofstream output;
	//
	//  Open the file.
	//
	output.open ( output_filename.c_str ( ) );

	if ( !output )
	{
		cerr << "\n";
		cerr << "R8MAT_WRITE - Fatal error!\n";
		cerr << "  Could not open the output file.\n";
		return;
	}
	//
	//  Write the data.
	//
	for ( j = 0; j < n; j++ )
	{
		for ( i = 0; i < m; i++ )
		{
			//			output << "\t" << setw(24) << setprecision(16) << table[i+j*m] ;
			output << table[i+j*m] << " " ;
		}
		output << "\n";
	}
	//
	//  Close the file.
	//
	output.close ( );

	return;
}

void r8mat_write2 ( string output_filename, int m, int n, double table[], vector<double> &lower_bound, vector<double> &upper_bound)

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the table data.
//
{
	int i;
	int j;
	ofstream output;
	//
	//  Open the file.
	//
	output.open ( output_filename.c_str ( ) );

	if ( !output )
	{
		cerr << "\n";
		cerr << "R8MAT_WRITE - Fatal error!\n";
		cerr << "  Could not open the output file.\n";
		return;
	}
	//
	//  Write the data.
	//
	for ( j = 0; j < n; j++ )
	{
		for ( i = 0; i < m; i++ )
		{
			output << "  " << setw(24) << setprecision(16) << rescale( lower_bound[i], upper_bound[i],table[i+j*m]) ;
		}
		output << "\n";
	}
	//
	//  Close the file.
	//
	output.close ( );

	return;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

	static char time_buffer[TIME_SIZE];
	const struct tm *tm;
	size_t len;
	time_t now;

	now = time ( NULL );
	tm = localtime ( &now );

	len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

	cout << time_buffer << "\n";

	return;
# undef TIME_SIZE
}

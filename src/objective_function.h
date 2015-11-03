/*
 * objective_function.h
 *
 *  Created on: Jan 19, 2014
 *      Author: ctsotskas
 */

#ifndef OBJECTIVE_FUNCTION_H_
#define OBJECTIVE_FUNCTION_H_

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
#include "container.h"
#include "configuration_settings.h"
#include "objective_function_formulae.h"
#include "test.cuh"

extern int FileExists(char const *strFilename);
extern void initialise_vector(Point2 &P, double value);

template <typename T>
class ObjectiveFunctionBasic{

public:
	unsigned int nVar;
	unsigned int nObj;
	int n_of_successful_calculations;
	int violations;

	Point2 min_bound; //minimum values for each design variable
	Point2 max_bound; //maximum
	Point2 range; //values' respective range

	ObjFunction2 failedObjectiveFunctionVector;

#ifdef xfoil
	//used for external ffd
	int n_vara  ;
	T x0_ffd[16];
	//Point x0_ffd;
	int ncp[2];
	T ffd_dpx[500];
	T ffd_dpy[500];
	int np;
	int ffd_flag, ffdflag;
	Point vars;




	int xfoil_flag, xfoilflag;

	T cl  , cd ;



	//used for external xfoil
#endif

#ifdef external
	char *variables_name, *objectives_name, *evaluator;
#endif

#ifdef external2
	char variables_name[500], objectives_name[500];
#endif


#ifdef user_defined
	int set_user_defined_bounds(Point2 &low, Point2 &high);
	ObjFunction2 eval_user_defined_function(Point2 &input);

#endif



	void read_bounds_from_file(char const *name, unsigned int const n_of_entries, Point2 &low_bound, Point2 &high_bound);

	//constructors
	//ObjectiveFunctionBasic();

	ObjectiveFunctionBasic(Overall_Optimisation_Configuration_Settings &conf2);
	ObjectiveFunctionBasic(Overall_Optimisation_Configuration_Settings &conf2, objective_function_formulae &formulae);
	double sum(const Point2 &Pnt);

	ObjFunction2 get_penalty_objectives();
	ObjFunction2 penalise(Point2 &input);

	ObjFunction2 eval_external_evaluator(Point2 &input);

	ObjFunction2 eval_external_evaluator2(Point2 &input);
	void eval_external_evaluator2_parallel(std::map<Point2,ObjFunction2> &inCont);
	void CPU_evaluation_inteface(Container2 &inCont);
	void CUDA_evaluation_interface(Container2 &inCont);
	void eval_internal_parallel(Container2 &inCont);

	ObjFunction2 calculateObjFun(Point2 &input);
	void calculateObjFun_parallel(Container2 &C1);
	int isValid ( Point2 &P);

private:
	objective_function_formulae __formulae;
};

#ifndef OBJECTIVES_ACCURACY
#define OBJECTIVES_ACCURACY
typedef ObjectiveFunctionBasic<double> ObjectiveFunction;
#endif

template<typename T>
ObjectiveFunctionBasic<T>::ObjectiveFunctionBasic(Overall_Optimisation_Configuration_Settings &conf2) :
nVar(conf2.getExternalConfigurationFile().getVar()),
nObj(conf2.getExternalConfigurationFile().getObj()),
n_of_successful_calculations(0),
violations(0),
min_bound(conf2.get_lower_bound()),
max_bound(conf2.get_upper_bound()),
range(nVar,0.0),
failedObjectiveFunctionVector(conf2.get_penalty_point()),
__formulae(nVar)
{



#if defined(ZDT1) || defined(ZDT2) || defined(ZDT3) || defined(ZDT4) || defined(ZDT6) || defined(DTLZ1) || defined(DTLZ2) || defined(DTLZ3) || defined(DTLZ4) || defined(DTLZ5) || defined(DTLZ6) || defined(DTLZ7)
	std::cout << "one of the built-in options was selected, so the design_vector_ranges.txt is not read!!!" << std::endl;
	std::cout << "\tone of the following is activated: ZDT1, ZDT2, ZDT3, ZDT4,  ZDT5, ZDT6, DTLZ1, DTLZ2, DTLZ3, DTLZ4, DTLZ5, DTLZ6, DTLZ7" << std::endl;

#else

#endif

	//new ones for xfoil
#ifdef xfoil
	//	initialise_vector(max_bound,0.3);
	//	initialise_vector(min_bound,-0.4);

	n_vara = 8;

	//x0_ffd = Point( n_vara + 8, 0.0);
	xfoil_flag=0;
	xfoilflag=0;
	cl=0.;
	cd=0.;
	vars = Point(8, 0.0);
#endif

#ifdef external

	//Dortmund External
	//	variables_name="./variables.txt";
	//	objectives_name="./objectives.txt";
	//	evaluator="./main_pc64-opteron-linux-gcc4.x";
	//	initialise_vector(max_bound,50);
	//	initialise_vector(min_bound,-50);

	//nurul External
	//	variables_name="FuelComponent.inp";
	//	objectives_name="Results.out";
	//	evaluator="BioEvaluat.exe";
	//	min_bound[0]=0.0;
	//	min_bound[1]=0.0;
	//	min_bound[2]=0.51;
	//
	//	max_bound[0]=100.0;
	//	max_bound[1]=0.0;
	//	max_bound[2]=0.7;





	//todo: fix it to read the design_vector_ranges.txt


	//			initialise_vector(max_bound,50);
	//	initialise_vector(min_bound,-50);


#endif

#ifdef external2



	//GATAC External
	//variables_name="toBeEvaluated.txt";
	//objectives_name="toBeEvaluated.txt.out";
	strcpy(variables_name,"toBeEvaluated.txt");
	strcpy(objectives_name,"toBeEvaluated.txt.out");


#endif

	for(unsigned int i=0; i<range.size(); ++i)
		range[i]=max_bound[i]-min_bound[i];
}

template<typename T>
ObjectiveFunctionBasic<T>::ObjectiveFunctionBasic(Overall_Optimisation_Configuration_Settings &conf2, objective_function_formulae &formulae) :
nVar(conf2.getExternalConfigurationFile().getVar()),
nObj(conf2.getExternalConfigurationFile().getObj()),
n_of_successful_calculations(0),
violations(0),
min_bound(conf2.get_lower_bound()),
max_bound(conf2.get_upper_bound()),
range(nVar,0.0),
failedObjectiveFunctionVector(conf2.get_penalty_point()),
__formulae(formulae)
{



	//read design vector ranges - START

#if defined(ZDT1) || defined(ZDT2) || defined(ZDT3) || defined(ZDT4) || defined(ZDT6) || defined(DTLZ1) || defined(DTLZ2) || defined(DTLZ3) || defined(DTLZ4) || defined(DTLZ5) || defined(DTLZ6) || defined(DTLZ7)
	std::cout << "one of the built-in options was selected, so the design_vector_ranges.txt is not read!!!" << std::endl;
	std::cout << "\tone of the following is activated: ZDT1, ZDT2, ZDT3, ZDT4,  ZDT5, ZDT6, DTLZ1, DTLZ2, DTLZ3, DTLZ4, DTLZ5, DTLZ6, DTLZ7" << std::endl;

#else

#endif

	//read design vector ranges - END



	//new ones for xfoil
#ifdef xfoil
	//	initialise_vector(max_bound,0.3);
	//	initialise_vector(min_bound,-0.4);

	n_vara = 8;

	//x0_ffd = Point( n_vara + 8, 0.0);
	xfoil_flag=0;
	xfoilflag=0;
	cl=0.;
	cd=0.;
	vars = Point(8, 0.0);
#endif

#ifdef external

	//Dortmund External
	//	variables_name="./variables.txt";
	//	objectives_name="./objectives.txt";
	//	evaluator="./main_pc64-opteron-linux-gcc4.x";
	//	initialise_vector(max_bound,50);
	//	initialise_vector(min_bound,-50);

	//nurul External
	//	variables_name="FuelComponent.inp";
	//	objectives_name="Results.out";
	//	evaluator="BioEvaluat.exe";
	//	min_bound[0]=0.0;
	//	min_bound[1]=0.0;
	//	min_bound[2]=0.51;
	//
	//	max_bound[0]=100.0;
	//	max_bound[1]=0.0;
	//	max_bound[2]=0.7;

	//todo: fix it to read the design_vector_ranges.txt
#endif

#ifdef external2
	//GATAC External
	//variables_name="toBeEvaluated.txt";
	//objectives_name="toBeEvaluated.txt.out";
	strcpy(variables_name,"toBeEvaluated.txt");
	strcpy(objectives_name,"toBeEvaluated.txt.out");
#endif

	for(unsigned int i=0; i<range.size(); ++i)
		range[i]=max_bound[i]-min_bound[i];
}

template<typename T>
double ObjectiveFunctionBasic<T>::sum(const Point2 &Pnt){
	double result=0;

	result = std::accumulate(Pnt.begin(),Pnt.end()-1,result);
	return result;
}


#ifdef user_defined
template<typename T>
int ObjectiveFunctionBasic<T>::set_user_defined_bounds(Point2 &low, Point2 &high){
	//alternative definition
	/*for(int i=0; i< low.size(); ++i)
low[i]=0.0;

for(int i=0; i< high.size(); ++i)
high[i]=1.0;*/

	initialise_vector(high,5);
	initialise_vector(low,-5);

	high[0]=1;
	low[0]=0;

	return 1;

}

//implement the constraints and the objectives of the model
template<typename T>
ObjFunction2  ObjectiveFunctionBasic<T>::eval_user_defined_function(Point2 &input){
	ObjFunction2 output(failedObjectiveFunctionVector);


	//example of ZDT4 implementation
	//	if( input.size() != 10){
	//		std::std::cout << " ZDT4 input has" << input.size() << "variables, instead of 10" << std::std::endl;
	//		exit(-100);
	//	}
	//	int i = 0;
	//	int n = input.size();
	//	double f1 = 0;
	//	double g = 0;
	//	double h = 0;
	//
	//
	//	assert(n == 10); //variable vector length
	//	assert(output.size() == 2);
	//
	//	f1 = input[0];
	//
	//	for (i = 1; i < n; i++)
	//	{
	//		double x = input[i];
	//		g += x * x - 10 * std::cos(4 * pi * x);
	//	}
	//	g = 1 + 10 * (n - 1) + g;
	//	h = 1 - std::sqrt(f1 / g);
	//
	//	output[0]=f1;
	//	output[1] = g * h;


	//trajectory - start
	if( input.size() != 14){
		std::std::cout << " trajectory input has" << input.size() << "variables, instead of 14" << std::std::endl;
		exit(-100);
	}

	if(simulate_trajectory(input, output)!=0)
		output=failedObjectiveFunctionVector;

	//trajectory - end


	//LBM2 - start






	//LBM2 - end



	return output;
}

#endif //user_defined


template<typename T>
ObjFunction2 ObjectiveFunctionBasic<T>::calculateObjFun(Point2 &input){
	ObjFunction2 output(get_penalty_objectives());

#ifdef DEBUG_EVALUATING
	std::cout << " evaluating " << input << std::endl;
#endif

#ifdef xfoil
	double f1 , f2;

	if(isValid(input) ) {


		for(  int i=0; i < 16; ++i)
			x0_ffd[i]=0;


		ncp[0]=4;
		ncp[1]=2;

		for(  int i=4; i <input.size()+4; ++i)
			x0_ffd[i] = input[i-4];



		std::cout << "<!-- This is a comment: &#160;  &nbsp <XFOIL>calling FFD and xfoil" << std::endl;


		pid_t mypid=getpid();
		char extn[8];
		gcvt(mypid,7,extn);
		std::cout << "extn... = " << extn << std::endl;
		//     sprintf(ff,ff1,extn);
		string ff1 = string("ffd.") + string(extn);
		const char* ff = ff1.c_str();


		ffd_flag = ffd_(x0_ffd[0], ncp, ffd_dpx[0], ffd_dpy[0],np, ffdflag);

		if(ffd_flag!=0)
			std::cout << "FFD - OK" << std::endl;
		else
			std::cout << "FFD - FAILED!!!" << std::endl;

		xfoil_flag = xfoil_(cl, cd, ff, xfoilflag, ffd_dpx[0], ffd_dpy[0], np);

		++n_of_successful_calculations;

		if (xfoilflag==0 && xfoil_flag==0 && cl==cl && cd==cd) {
			std::cout << "XFOIL successful..." << std::endl;
			f1  = -(cl/1.4644);
			f2  = cd/0.03051;
		} else {

			std::cout << "XFOIL NOT successful..." << std::endl;
			std::cout << "about to call InvalidDesignVector" << std::endl;
			//throw cppot::InvalidGeometryException();
			++violations;
			return failedObjectiveFunctionVector;
		}
		std::cout << "cl = " << f1  << std::endl;
		std::cout << "cd = " << f2  << std::endl;
		std::cout << "interface xfoil flag = " << xfoilflag << std::endl;
		std::cout << "interface xfoil flag2 = " << xfoil_flag << std::endl;
		std::cout << "interface ffd flag = " << ffdflag << std::endl;


	}else{
		++violations;
		return failedObjectiveFunctionVector;
	}

	output[0]=f1;
	output[1]=f2;

#endif

	if(isValid(input)==0 ) {
		++violations;
		return failedObjectiveFunctionVector;
	}

#ifdef external
	//Dortmund and Nurul
	//		output=eval_external_evaluator(input );

	//xfoil
	//		output=eval_user_defined_function(input );

	std::cout << "evaluating design vector:" << input << std::endl;

	std::ofstream conf_file("LBM2_configuration.txt");

	conf_file << "4000" << std::endl;
	conf_file << "2500" << std::endl;
	conf_file << "0.0175" << std::endl;
	conf_file << input[0] << std::endl;
	conf_file << input[1] << std::endl;
	conf_file << input[2] << std::endl;
	conf_file << "59" << std::endl;
	conf_file << "512" << std::endl;
	conf_file << "interface_file" << std::endl;



	conf_file.close();


	(void)system("./LBM2");

	double buffer1, buffer2;
	std::ifstream objectives_file("LBM2_objectives.txt");

	objectives_file >> buffer1;
	objectives_file >> buffer2;

	std::cout << "read objectives: " << buffer1 << "\t" << buffer2 << std::endl;

	output[0]= (-1) * buffer1 / 6018.85; //vorticity
	output[1]= buffer2 / 38.9915; //pressure drop

	std::cout << "normalised objectives: " << output << std::endl;



	objectives_file.close();

#endif

#ifdef external2
	output=eval_external_evaluator2(input );
#endif

#ifdef user_defined
	output=eval_user_defined_function(input );
#endif

	++n_of_successful_calculations;
//	output=__formulae(input);
	return output;
}


template<typename T>
void ObjectiveFunctionBasic<T>::calculateObjFun_parallel(Container2 &C1){

#ifdef external2
	eval_external_evaluator2_parallel(C1);
#endif //external2
#ifndef external//ZDT
	eval_internal_parallel(C1);
#endif //ZDT
	//todo if there is a constraint based on objective function values, it goes HERE!

}

template<typename T>
int  ObjectiveFunctionBasic<T>::isValid ( Point2 &P){


	for(unsigned int i=0; i<P.size(); ++i){
		//std::cout << " checks margins!" << min_bound[i] << " " << P[i] << " " << max_bound[i] << std::endl;
		if( (P[i]>max_bound[i]) || (P[i]<min_bound[i])){
			//			std::cout << "validity:FAILED" << std::endl;


#ifdef DEBUG_EVALUATING
			std::cout << "design " << P << " presents a violation at the index" << i << std::endl;
#endif//DEBUG_EVALUATING



			std::cout<<"invalid, design variable" << i << " is out of range" << std::endl;


			return 0;
		}
	}

	//LBM constraints
	double r_small=P[0];
	double D_big=70.0;
	double s_min=2 * r_small +1.0 ;
	double s_max=  D_big- r_small - 1.0 ;
	double S_temp=P[2];
	if( S_temp<s_min || S_temp>s_max ){
		std::cout << "violated constraints" << std::endl;
		std::cout << S_temp <<  " < " << s_min << std::endl;
		std::cout << S_temp <<  " > " << s_max << std::endl;

		return 0;
	}

	//std::cout << "validity:OK !" << std::endl;
	return 1;
	//write up any constraints
}

///////////////////////////////////////////////////////////////////////
//////  Implementation of  ObjectiveFunctionBasic  END
///////////////////////////////////////////////////////////////////////


#ifndef OBJECTIVES_ACCURACY
#define OBJECTIVES_ACCURACY
typedef ObjectiveFunctionBasic<double> ObjectiveFunction;
#endif





template<typename T>
ObjFunction2 ObjectiveFunctionBasic<T>::get_penalty_objectives(){

#ifdef DEBUG_EVALUATING
	//	cout << "design vector " << input << " was penalised " << endl;
#endif//DEBUG_EVALUATING

	return failedObjectiveFunctionVector;
}


template<typename T>
ObjFunction2 ObjectiveFunctionBasic<T>::penalise(Point2 &input){

#ifdef DEBUG_EVALUATING
	std::cout << "design vector " << input << " was penalised " << std::endl;
#endif//DEBUG_EVALUATING

	return failedObjectiveFunctionVector;
}



#ifdef external
template<typename T>
ObjFunction2  ObjectiveFunctionBasic<T>::eval_external_evaluator(Point2 &input){
	//originally, before March 2012, ekana ena obj fun vector me 0.0 kai meta to evaza megalh timh

	//change on 22 March 2012
	ObjFunction2 output(get_penalty_objectives());


	//produce variable file
	std::ofstream file_out(variables_name);

	for(unsigned int i=0; i<input.size(); ++i)
		file_out << input[i] << " ";

	file_out << std::endl;
	file_out.close();

	//evaluate through external evaluator
	std::cout << "waiting external objective"<< std::endl;
	(void)system(evaluator);


	//wait
	do{
#if defined(_WIN32) || defined(_WIN64) || defined(__WIN32) || defined(WINDOWS) && !defined(__CYGWIN__)
		Sleep(READ_TIME); // micro-seconds
#endif

#if defined(__gnu_linux__)   || defined(__linux__)
		//TODO:activate the following line:
		usleep( 3000000 ); // micro-seconds
#endif
	}while( !FileExists(objectives_name) );


	//read objectives file
	std::ifstream file_in(objectives_name);

	for(unsigned int i=0; i<output.size(); ++i)
		file_in >> output[i];

	file_in.close();

	//delete objectives file
	remove(objectives_name);

	return output;
}
#endif //external

#ifdef xfoil
template<typename T>
ObjFunction2  ObjectiveFunctionBasic<T>::eval_user_defined_function(Point2 &input){
	ObjFunction2 output(failedObjectiveFunctionVector);


	/*	//example of ZDT4 implementation
if( input.size() != 10){
cout << " ZDT4 input has" << input.size() << "variables, instead of 10" << endl;
exit(-100);
}
	 */
	int i = 0;
	int n = input.size();
	/*	double f1 = 0;
double g = 0;
double h = 0;


assert(n == 10); //variable vector length
assert(output.size() == 2);

f1 = input[0];
cout << "var[0] = " << input[0] << endl;
for (i = 1; i < n; i++)
{
double x = input[i];
g += x * x - 10 * cos(4 * pi * x);
cout << "var[" << i << "] = " << input[i] << endl;
}
g = 1 + 10 * (n - 1) + g;
h = 1 - sqrt(f1 / g);

output[0]=f1;
output[1] = g * h;
	 */
	const int n_vara=8;
	//    const int n_vara=nreal;
	int ncp[2];
	//    double init_vars[n_vara];
	double x0_ffd[n_vara+8];
	//    double ffd_dpx[100];
	double ffd_dpx[500];
	//    double ffd_dpy[100];
	double ffd_dpy[500];
	//    float ffd_dp[100][2];
	int ffd_flag=0, ffdflag=0;
	int xfoil_flag=0, xfoilflag=0;
	double cl=0., cd=0.;
	int np;


	ncp[0]=4;
	ncp[1]=2;


	for (unsigned int i=0; i<n_vara+8; i++) {

		x0_ffd[i] = 0.0;
	}

	for (int i=0; i<n_vara; i++) {

		//                    x0_ffd[i+4] = init_vars[i];
		x0_ffd[i+4] = input[i];
		//                    cout << "x0_ffd[" << i+4 << "] = " << x0_ffd[i+4] << endl;
	}

	//     char ff1[30]="ffd%s.def";
	//     char ff[]="";
	//     char remov1[50]="rm ffd%s.def";
	//     char remov[]="";
	pid_t mypid=getpid();
	char extn[8];
	gcvt(mypid,7,extn);
	//     sprintf(ff,ff1,extn);
	//     sprintf(remov,remov1,extn);
	string ff1 = string("ffd.") + string(extn);
	const char* ff = ff1.c_str();

	ffd_flag = ffd_(x0_ffd[0], ncp, ffd_dpx[0], ffd_dpy[0],np, ffdflag);
	//	    ffd_flag = ffd_(x0_ffd, ncp, ffd_dp[0][0]);
	//            cout << "np from ffd... = " << np << endl;

	//    ofstream ffd_vector("ffd_deformed", ios::out);

	xfoil_flag = xfoil_(cl, cd, ff, xfoilflag, ffd_dpx[0], ffd_dpy[0], np);
	//     xfoil_flag = xfoil_(cl, cd, ff, xfoilflag);

	if (xfoilflag==0 && xfoil_flag==0 && isfinite(cl/cd) != 0) {
		cout << "XFOIL successful..." << endl;
		output[0] = -(cl/1.4644);
		output[1] = cd/0.03051;
	} else {
		cout << "XFOIL NOT successful..." << endl;
		output[0] = 100;
		output[1] = 100;

	}

	cout << "obj1 = " << output[0] << endl;
	cout << "obj2 = " << output[1] << endl;


	return output;
}


#endif //external

#ifdef external2




template<typename T>
ObjFunction2  ObjectiveFunctionBasic<T>::eval_external_evaluator2(Point2 &input){

	ObjFunction2 output(failedObjectiveFunctionVector);


	//produce variable file
	ofstream file_out(variables_name);

	for(unsigned int i=0; i<input.size(); ++i)
		file_out << input[i] << " ";

	file_out << endl;
	file_out.close();

	//evaluate through external evaluator
	cout << "EXT EVLTR2 waiting external objective"<<endl;



	//wait
	do{
#if defined(_WIN32) || defined(_WIN64) || defined(__WIN32) || defined(WINDOWS) && !defined(__CYGWIN__)
		Sleep(READ_TIME); // micro-seconds
#endif

#if defined(__gnu_linux__)   || defined(__linux__)
		usleep(READ_TIME); // micro-seconds
#endif
	}while( !FileExists(objectives_name) );


	//read objectives file
	ifstream file_in(objectives_name);

	for(unsigned int i=0; i<output.size(); ++i)
		file_in >> output[i];

	file_in.close();

	//delete objectives file
	remove(objectives_name);

	return output;
}


template<typename T>
void  ObjectiveFunctionBasic<T>::eval_external_evaluator2_parallel(std::map<Point2,ObjFunction2> &inCont){
#if defined(_WIN32) || defined(_WIN64) || defined(__WIN32) || defined(WINDOWS) && !defined(__CYGWIN__)
	Sleep(READ_TIME/2); // milliseconds
#endif
	//produce variable file
	ofstream file_out(variables_name);
	for(Container2::iterator it=inCont.begin(); it!=inCont.end(); ++it)
		file_out << it->first << endl;

	file_out << endl;
	file_out.close();

	//evaluate through external evaluator
	cout << "EXT EVLTR2 || waiting external objective"<<endl;



	//wait
	do{
#if defined(_WIN32) || defined(_WIN64) || defined(__WIN32) || defined(WINDOWS) && !defined(__CYGWIN__)
		Sleep(READ_TIME/2); // milliseconds
#endif

#if defined(__gnu_linux__)   || defined(__linux__)
		usleep(READ_TIME); // micro-seconds
#endif
	}while( !FileExists(objectives_name) );


	//read objectives file, once created
	ifstream file_in(objectives_name);

	for(Container2::iterator it=inCont.begin(); it!=inCont.end(); ++it)
		for(unsigned int i=0; i<nObj; ++i)
			file_in >> it->second[i];

	file_in.close();

	//delete objectives file
	remove(objectives_name);


}
#endif //external2

template<typename T>
void  ObjectiveFunctionBasic<T>::CPU_evaluation_inteface(Container2 &inCont){
	for(std::map<Point2,ObjFunction2>::iterator it=inCont.begin(); it!=inCont.end(); ++it){
#ifdef DEBUG_EVALUATING
		std::cout << " evaluating " << input << std::endl;
#endif //DEBUG_EVALUATING
		++n_of_successful_calculations;
		it->second=__formulae(it->first);

	}

}

template<typename T>
void  ObjectiveFunctionBasic<T>::CUDA_evaluation_interface(Container2 &inCont){
#ifdef GPU_evaluation
	CUDA_evaluation CE(16,inCont);
	CE.run();
	CE.update_container_objectives(inCont);
	n_of_successful_calculations+=inCont.size();
#endif //GPU_evaluation
}

template<typename T>
void  ObjectiveFunctionBasic<T>::eval_internal_parallel(Container2 &inCont){

#ifdef CPU_evaluation
	CPU_evaluation_inteface(inCont);
#endif //CPU_evaluation

#ifdef GPU_evaluation
	CUDA_evaluation_interface(inCont);
#endif //GPU_evaluation
}



#endif /* OBJECTIVE_FUNCTION_H_ */

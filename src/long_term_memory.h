/*
 * long_term_memory.h
 *
 *  Created on: Jan 19, 2014
 *      Author: ctsotskas
 */

#ifndef LONG_TERM_MEMORY_H_
#define LONG_TERM_MEMORY_H_

#include <iostream>
#include <cstdlib>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <deque>
#include <numeric>
#include <sys/stat.h> //file_exists
//#include <cmath>
#include <cassert>
#include <ctime>
#include <fstream>
#include <istream>
//#include <string>
#include <cstring>
#include <sstream>

#include "global_defines.h"

//todo the following bits should be replaced and refactored as a template within the header file
extern void correct_accuracy(Point2 &P) ;

#define round(r) ( (r) > 0.0) ? floor( (r) + 0.5) : ceil( (r) - 0.5);

extern Point2 operator- (const Point2 &v1, const Point2 &v2);
extern Point2 operator+ (const Point2 &v1, const Point2 &v2);
extern Point2 operator/ (const Point2 &v1, const Point2 &v2);
extern Point2 operator/ (const Point2 &v1, const int N);
extern bool operator==(const Point2 &P1, const Point2 &P2);
extern bool operator>(const Point2 &P1, const Point2 &P2);
extern bool operator<(const Point2 &P1, const Point2 &P2);
extern int FileExists(char const *strFilename);

template <typename T>
class LTM_Container2Basic2 :  public std::set< Point2 >{
	unsigned int local_nVar;
	unsigned int local_nRegions; //todo should be const
	Point2 local_lower_bound, local_upper_bound;
	Point2 tempVector;
	T offset;
	Point2 range_points;
	std::vector<Point2 > variable_range_points;
	std::string container_name;
	std::string __save_path;
	LTM_Container2Basic2 & operator= ( LTM_Container2Basic2 const &other);
public:
	std::deque< Point2 > Region; //array of regions, each element corresponds to one design variable
	LTM_Container2Basic2(const int &nVar,const int &nRegion, const Point2 &lower_bound, const Point2 &upper_bound, std::string given_name, std::string filename );
	~LTM_Container2Basic2();

	T RandomNumber(T min, T max);
	//Point2 generate_Random_Point_From_Least_Visited_Region();
	Point2 generate_Random_Point_From_Least_Visited_Region2();
	int leastVisitedRegion();
	void showContents();
	int getSize();
	int getSizePerRegion();
	std::vector<int> getRespetiveRegion(const Point2 &P);
	void save_ltm_container(std::string case_qualifier, int evaluations_counter);
	void load_ltm_container(char const *save_path);
};

typedef LTM_Container2Basic2<double> LTM_Container;

template <typename T>
LTM_Container2Basic2<T>::LTM_Container2Basic2(const int &nVar,const int &nRegion, const Point2 &lower_bound, const Point2 &upper_bound, std::string given_name, std::string filename ) :
local_nVar(nVar),
local_nRegions(nRegion),
local_lower_bound(Point2(lower_bound)),
local_upper_bound(Point2(upper_bound)),
tempVector(Point2(nVar,0.0)),
offset(0.00001),  //DBL_MIN ; //0.0000000000000001;  //DBL_MIN is 1E-37
// nRegions are defined by nRegions+1 points!!!, each regions is identified by its lower_bound point
//i.e. 3rd region should be chosen by range_points[2]
range_points(Point2(local_nRegions + 1, 0.0)),
variable_range_points(std::vector<Point2>(nVar, range_points)), // [variable_index][region_index]
container_name(given_name),
__save_path(filename),
Region (std::deque<Point2>(local_nRegions, lower_bound ))
{


	Point2 local_regions_vector(nVar,local_nRegions);
	Point2 local_step(nVar,local_nRegions);

	for(unsigned int i=0;i<local_step.size(); ++i)
		local_step[i]=(local_upper_bound[i]-local_lower_bound[i]) / local_regions_vector[i];


	//	std::cout << "DBG: received:" << std::endl;
	//	std::cout << "lower bound" << local_lower_bound << std::endl;
	//			std::cout<< "upper bound"  << local_upper_bound << std::endl;
	//			std::cout << "local regions" <<  local_regions_vector << std::endl;
	//			std::cout << "local  step" << local_step << std::endl;
	//			std::cout << " calculation " << (local_upper_bound-local_lower_bound) / local_regions_vector  << std::endl;
	//			std::cout << " calculation2 " << (local_upper_bound-local_lower_bound)   << std::endl;
	//			std::cout << " calculation3 " << Point2(nVar,local_nRegions)   << std::endl;


	//		for(int i=0;i<nVar;++i)
	//			local_range[i]=local_upper_bound[i]-local_lower_bound[i];

	srand(time(NULL));

	std::cout << "so far: " << std::endl;
	std::cout << "number of variables for each design "
			<< variable_range_points.size() << std::endl;
	std::cout << "\t each variable has " << local_nRegions << " regions, defined by "
			<< range_points.size() << " points " << std::endl;

	//create regions
	for (int v = 0; v < nVar; ++v) {
		for (unsigned int r = 0; r < range_points.size(); ++r) {
			range_points[r] = lower_bound[v] + local_step[v] * r;
//			std::cout << "DBG:  v:"<< v << " r:" << r << ": " <<  range_points[r]  << " " << lower_bound[v] << " "<< local_step[v] << " " << r << std::endl;
		}
		variable_range_points[v] = range_points;
		//		variable_range_points.push_back(range_points);
	}

	//report the regions
	for(int v=0; v<nVar; ++v) {
		std::cout<< std::endl << "variable " << v << " range_points: " <<std::endl;
		for(unsigned int r=0; r<range_points.size(); ++r)
			std::cout << " \t" << variable_range_points[v][r];
	}
	std::cout << std::endl;


}

template <typename T>
LTM_Container2Basic2<T>::~LTM_Container2Basic2(){
	std::cout <<"LTM memories deleted!" << std::endl;
}

template <typename T>
T LTM_Container2Basic2<T>::RandomNumber(T min, T max){
	T r = rand() /  ((T)RAND_MAX );
	return (T)(min + r * (max - min));
}





template <typename T>
Point2 LTM_Container2Basic2<T>::generate_Random_Point_From_Least_Visited_Region2(){
	//todo , improve!!! na diavazei apo to HISTORY, OXI apo to LTM!


	//classify the points

	//this counts the occurrences for each region!
	//this is not the same as the range_points above!!!!
	std::vector<int> zero_vector(local_nRegions, 0);
	std::vector<std::vector<int> > region_counters(local_nVar, zero_vector); // [variable_index][region_index]

	//for each entry
	for (const_iterator it = begin(); it != end(); ++it) {
		//for each component of the entry
		for (unsigned int v = 0; v < local_nVar; ++v) {
			for (unsigned int r = 0; r < local_nRegions; ++r)
				if (variable_range_points[v][r] <= (*it)[v]  && (*it)[v] <= variable_range_points[v][r + 1]) {
					++region_counters[v][r];
#ifdef DEBUG
					std::cout << "\t" << (*it)[v] << " belongs to" << variable_range_points[v][r] << " and " << variable_range_points[v][r+1] << std::endl;
#endif
				}
		}
#ifdef DEBUG
		std::cout << std::endl;
#endif
		//choose where does it belong

		//		variable_range_points[v][r];

	}
	//report intermediate classification results
	std::cout << "classification results1" << std::endl;
	for (unsigned int r = 0; r < local_nRegions; ++r) {
		std::cout << " region " << r << " counted: " << std::endl;
		for (unsigned int v = 0; v < local_nVar; ++v) {
			std::cout << " \t" << region_counters[v][r]
			                                         << " Occurrences of variable " << v << std::endl;
		}
	}

	std::cout << "classification results2" << std::endl;
	for (unsigned int v = 0; v < local_nVar; ++v) {
		std::cout << " frequencies for variable " << v << " : " << std::endl;
		for (unsigned int r = 0; r < local_nRegions; ++r) {
			std::cout << " \t" << region_counters[v][r] << " Occurrences in region "
					<< r << std::endl;
		}
	}
	//find minimum region
	std::vector<int> minimum_region(local_nVar, -1);

	for (unsigned int v = 0; v < local_nVar; ++v) {

		int min = region_counters[v][0];
		int temp_min_reg = 0;
		//		std::cout << "v:" << v << "r=0 f:" << region_counters[v][0] << std::endl;
		for (unsigned int r = 1; r < local_nRegions; ++r) {
			//			std::cout << "v:" << v << "r=" << r << " f:" << region_counters[v][r] << std::endl;
			if (region_counters[v][r] < min) {
				min = region_counters[v][r];

				temp_min_reg = r;
			}
		}

		minimum_region[v] = temp_min_reg;
	}
	//report minimum region
	std::cout << "minumum region is" << std::endl;

	for (unsigned int v = 0; v < local_nVar; ++v)
		std::cout << " \t " << minimum_region[v];
	std::cout << std::endl;

	//generate a point from the minimum region

	std::cout << "generated_point" << std::endl;
	std::vector<double> output(local_nVar, 0.0);
	for (unsigned int v = 0; v < local_nVar; ++v) {
		int temp_min_reg = minimum_region[v];
		output[v] = RandomNumber(variable_range_points[v][temp_min_reg],
				variable_range_points[v][temp_min_reg + 1]);
		std::cout << "\t" << output[v];
	}
	std::cout << std::endl;
	return output;
}



template <typename T>
int LTM_Container2Basic2<T>::leastVisitedRegion(){
	int key=0;
	unsigned int min=Region[0].size();
	for (unsigned int i=1 ; i < local_nRegions ; ++i)
		if(Region[i].size()<min)
			key=i;

	return key;
}

//template <typename T>
//void LTM_Container2Basic2<T>::showContents(){
//	int m;
//	int total_size=0;
//	std::cout << "\t LTM contents" << std::endl;
//	for(int i=0; i< local_nRegions; ++i)
//		total_size+=Region[i].size();
//
//	std::cout << "\t ltm size=" << total_size << std::endl;
//
//
//	for(int r=0; r< local_nRegions; ++r){
//		std::cout << std::endl<< "Region =" << r << std::endl;
//		m=1;
//		for( std::deque< Point2 >::iterator it=Region[r].begin() ; it!=Region[r].end() ;++it, ++m)
//			std::cout << m << ".\t" << ( *it) << std::endl;
//
//	}
//
//
//}

template <typename T>
int LTM_Container2Basic2<T>::getSize(){
	int total_size=0;
	for (unsigned int i=0; i< local_nRegions ; ++i)
		total_size+=Region[i].size();

	return total_size;
}

template <typename T>
int LTM_Container2Basic2<T>::getSizePerRegion(){
	int total_size=0;
	for (int i=0; i< local_nRegions ; ++i){
		std::cout << i << ". " << Region[i].size() << std::endl;
		total_size+=Region[i].size();
	}

	return total_size;
}

template <typename T>
std::vector<int> LTM_Container2Basic2<T>::getRespetiveRegion(const Point2 &P){

	std::vector<int> region_point( P.size(), 0);

	for (unsigned int i=0; i< P.size() ; ++i){
		for(int j=0; j<local_nRegions; ++j)
			region_point[i]=j;
	}
	return region_point;
}


template <typename T>
void LTM_Container2Basic2<T>::save_ltm_container(std::string case_qualifier, int evaluations_counter){
	//create the ./memories/LTM_snapXXXXX.txt , where XXXX the current number of evaluations

	std::ostringstream dbg_message, snapshot_name, simple_name;
	std::string temp_s (__save_path+"/" + case_qualifier +  "_" + container_name + "_snap");

	snapshot_name << __save_path+"/" << case_qualifier <<  "_" << container_name << "_snap" << evaluations_counter << ".txt";
	simple_name << __save_path+"/" << case_qualifier <<  "_" << container_name  << ".txt";

	std::ofstream LTMfile(snapshot_name.str().c_str());
	std::ofstream simple_file(simple_name.str().c_str());

	for(  std::set< Point2 >::iterator it=begin() ;  it!=end() ; ++it){ //ola ta shmeia se ka8e deque
		LTMfile << *it  << std::endl;
		simple_file << *it  << std::endl;
	}

	LTMfile.close();
	simple_file.close();

#ifdef DEBUG
	dbg_message<<__save_path<<" saved!";
	std::cout << dbg_message.str() << std::endl;
#endif
}


template <typename T>
void LTM_Container2Basic2<T>::load_ltm_container(char const *save_path){
	std::ifstream LTMfile;
	LTMfile.open ( save_path );
	//T temp_buffer;
	Point2 temp_bufferV(local_nVar, 0.0);
	unsigned int counter=0;

	while( LTMfile>>temp_bufferV[0] ){
		for (unsigned int v=1; v<local_nVar; ++v)
			LTMfile>>temp_bufferV[v];

		insert( temp_bufferV );
		++counter;
	}

	std::cout << " totally " << counter << " designs were loaded into LTM" << std::endl;


	LTMfile.close();
}


typedef LTM_Container2Basic2<double> LTM_Container;

#endif /* LONG_TERM_MEMORY_H_ */

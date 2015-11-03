/*
 * container.cpp
 *
 *  Created on: Jan 19, 2014
 *      Author: ctsotskas
 */
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

#include <cstring>
#include <sstream>

#include "container.h"

Container2::Container2(unsigned int n, unsigned int m, std::string given_name,
		std::string save_path):
		local_nVar(n),
		local_nObj(m),
		container_name(given_name),
		__save_path(save_path){
	if( container_name.compare("MTM")==0   || container_name.compare("IM")==0 || container_name.compare("HISTORY")==0)
		std::cout << "container " << container_name << " has just been created, nVar=" << local_nVar << " and nObj=" << local_nObj << std::endl;

	srand((unsigned)(time(0)));
}

Container2::~Container2(){
	if( container_name.compare("MTM")==0   || container_name.compare("IM")==0 || container_name.compare("HISTORY")==0 )
		std::cout << container_name  << " terminated " << std::endl;
}

//Algorithm 2
int Container2::dominates(const ObjFunction2 &O1, const ObjFunction2 &O2){ //container
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

//Algorithm 3
int Container2::removeDominatedPoints(){
	int number_of_points_prior_removal=size();
	iterator itA, itB;
	std::deque<Point2> toBeRemoved;
	for(itA=begin() ; itA!=end() ; ++itA)
		for(itB=begin() ; itB!=end() ; ++itB)
			if( !(itA->first == itB->first) && dominates(    itA->second   ,     itB->second )==1 )
				toBeRemoved.push_back(itB->first);
	// remove point from container

	//updates in memory cannot be applied on the fly; that's why I log the points to be removed
	for(std::deque<Point2>::iterator itC=toBeRemoved.begin()  ; itC!=toBeRemoved.end() ; ++itC )
		erase(*itC);

	return number_of_points_prior_removal-toBeRemoved.size();
}

//Algorithm 5
int Container2::addIfNotDominated(const Point2 &V1,const ObjFunction2 &O1){
#ifdef DEBUG
	if(V1.empty()){
		std::cout << "variable 1 is empty!" << std::endl;
		exit(-111);
	}
	if(O1.empty()){
		std::cout << "objective 1 is empty!" << std::endl;
		exit(-222);
	}
#endif
	for (iterator it=begin() ; it!= end() ; ++it)
		if( dominates(  it->second  ,  O1 )==1)  //first is the existing point
			return 0;
	insert( tuple(V1, O1) );
	return 1;
}

Point2 Container2::selectRandom(){//container
	iterator it;

	int jump=rand() % size();
#ifdef DEBUG
	std::cout << "[Cont] among" << size() << " points I have choosen the " << jump+1 << "th element" << std::endl;
#endif
	advance(it=begin(), jump);
	return it->first;
}

Point2 Container2::find_extreme(int objective){
	const_iterator it=begin();

	Point2 proposed_vector=it->first;
	double min_objective=it->second[objective];
	++it;

	for(  ; it!=end(); ++it ){
		if(it->second[objective]<min_objective ){
			min_objective=it->second[objective];
			proposed_vector=it->first;

		}
	}
	return proposed_vector;
}

void Container2::displayContainer (){
	int i=1;
	std::cout << "container display" << std::endl;
	for(iterator it=begin() ; it != end(); ++it ){
		std::cout << i << ". " << it->first << "\t" << it->second <<std::endl;
		i++;
	}
}

Point2 Container2::check_memory(){
	Point2 empties(2,0);
	int empty_variables_counter=0, empty_objective_counter=0;
	for (Container2::const_iterator it1=begin(); it1!=end(); ++it1 ){
		if(it1->first.size()==0){
#ifdef DEBUG
			std::cout << "empty variable" << std::endl;
#endif
			++empty_variables_counter;
		}
		if(it1->second.size()==0){
#ifdef DEBUG
			std::cout << "empty objective" << std::endl;
#endif
			++empty_objective_counter;
		}
	}
#ifdef DEBUG
	std::cout << "\tstatus:" << empty_variables_counter << "variables and" << empty_objective_counter << " objectives are empty!" << std::endl;
#endif
	if(empty_variables_counter!=0 ||  empty_objective_counter!=0){
		std::cout <<" empty variables or objectives exist in memories!" <<std::endl;
		exit(-4444);
	}

	empties[0]=empty_variables_counter;
	empties[1]=empty_objective_counter;
	return empties;
}

int Container2::save_container_snapshot(std::string case_qualifier, int evaluations_counter){
	//create the ./memories/MTM_snapXXXXX.txt , where XXXX the current number of evaluations
	//create the ./memories/IM_snapXXXXX.txt , where XXXX the current number of evaluations
	//create the ./memories/HISTORY_snapXXXXX.txt , where XXXX the current number of evaluations

	std::ostringstream dbg_message, snapshot_name, simple_name;
	std::string temp_s (__save_path+"/" + case_qualifier +  "_" + container_name + "_snap");

	snapshot_name << __save_path+"/" << case_qualifier <<  "_" << container_name << "_snap" << evaluations_counter << ".txt";
	simple_name << __save_path+"/" << case_qualifier <<  "_" << container_name  << ".txt";

	std::ofstream SAVEfile(snapshot_name.str().c_str()), simple_file(simple_name.str().c_str());

	int entries_saved=0;
	for(iterator it=begin() ; it != end() ; ++it){
		SAVEfile  << it->first << "\t";
		SAVEfile  << it->second << "\t" << std::endl;

		simple_file << it->first << "\t";
		simple_file  << it->second << "\t" << std::endl;
		++entries_saved;
	}
	SAVEfile.close();
	simple_file.close();
	//this part just copies the same file in a plain format, this will be used from restart (should it is invoked)

#ifdef DEBUG
	dbg_message<<__save_path<<" saved!";
	std::cout << dbg_message.str() << std::endl;
#endif
	return entries_saved;
}


int Container2::report_optima_plot(char const *save_path){
	std::ofstream SAVEfile;
	SAVEfile.open ( save_path ,  std::ios::out  );

	int entries_saved=0;
	for(iterator it=begin() ; it != end() ; ++it){
		SAVEfile  << it->second << "\t" << std::endl;
		++entries_saved;
	}
	SAVEfile.close();
	return entries_saved;
}

int Container2::report_optima_current(char const *save_path){
	std::ofstream SAVEfile;
	SAVEfile.open ( save_path ,  std::ios::out  );

	int entries_saved=0;
	for(iterator it=begin() ; it != end() ; ++it){
		SAVEfile  << " variable vector "<< it->first << "\t ; objectives vector";
		SAVEfile  << it->second << std::endl;
		++entries_saved;
	}
	SAVEfile.close();
	return entries_saved;
}

void Container2::add_point(const entry __entry){
	insert( __entry );
}


int Container2::load_container(char const *load_path){
	unsigned int i;
	int entries_saved=0;
	double buffer;
	std::ifstream LOADfile;
	LOADfile.open ( load_path  );
	Point2 tempPoint(local_nVar, 0.0);
	ObjFunction2 tempObjFun(local_nObj, 0.0);


	std::cout << " check point is reading from" << load_path << " file "<< std::endl;
	std::cout << "\t" << tempPoint.size() << " variables " << tempObjFun.size() << " objectives... " << std::endl;

	while( LOADfile >> buffer){
		//std::cout << "in loop" <<std::endl;
		for(i=0; i<tempPoint.size()  ; ++i){
			//reads one element at a time (ON THE SAME LINE)
			tempPoint[i]=buffer;
			//bufferVector.push_back(buffer); //feeds the above element into the vector, until the line "entry" is replicated into the vector
			LOADfile >> buffer;
		} //that's for the key

		for(i=0; i<tempObjFun.size()-1 ; ++i){
			//reads one element at a time (ON THE SAME LINE)
			tempObjFun[i]=buffer;
			//bufferVector2.push_back(buffer); //feeds the above element into the vector, until the line "entry" is replicated into the vector
			LOADfile >> buffer;

		}
		tempObjFun[i]=buffer;
		insert( entry( tempPoint,tempObjFun ) );
	}

	LOADfile.close();
	std::cout << "\ttotally" << size() << " points loaded" << std::endl;

#ifdef REPORT_CONTAINERS
	displayContainer();
	std::cout << "load finished!" <<std::endl;
#endif
	return entries_saved;
}

double Container2::calculate_hypervolume2D(const ObjFunction2 &reference_point){
	ObjFunction2 tempv1, tempv2;
	std::set<ObjFunction2> counter_table;

	if( size()==0){
		std::cout<<"there are not points to produce HV" << std::endl;
		std::cout<< container_name <<"'s size:" << size() << std::endl;
		exit(-11919);
	}

	for(iterator it=begin(); it != end() ; ++it ){
		counter_table.insert(it->second);
	}
	//uncomment for 2D HV
	double reference_point_x, reference_point_y;
	double newpoint_x, newpoint_y;
	double previous_point_x, previous_point_y;
	ObjFunction2 temp1;
	double hypervolume=0.0;
	reference_point_x= reference_point[0];
	reference_point_y= reference_point[1];

	std::set<ObjFunction2>::iterator it1=counter_table.begin();

	newpoint_x=(*it1)[0];
	newpoint_y=(*it1)[1];
	++it1;

#ifdef DEBUG
	std::cout << "\tnew point values are " << newpoint_x << "," << newpoint_y << std::endl;
#endif

	hypervolume+=fabs(newpoint_y - reference_point_y ) * fabs(newpoint_x - reference_point_x );

	previous_point_x = newpoint_x;
	previous_point_y = newpoint_y;

	while(it1!=counter_table.end()){

		++it1;
		newpoint_x=(*it1)[0]; //update point
		newpoint_y=(*it1)[1];

		temp1[0]=newpoint_x;
		temp1[1]=newpoint_y;

		if(previous_point_x!=newpoint_x)
			hypervolume+=fabs(newpoint_y - previous_point_y ) * fabs(newpoint_x - reference_point_x );


		previous_point_x = newpoint_x;
		previous_point_y = newpoint_y;
	}

	return hypervolume;
}

double Container2::calculate_quality_indicator(const ObjFunction2 &reference_point){
	//container could be const, but I do not change it either way!
#ifdef DEBUG
	std::cout << "hyper volume calculation START" << std::endl;
	std::cout << "\tPoints to evaluate"<< size() << std::endl;
	std::cout << "\t with Reference point set to " << reference_point << std::endl;

#endif
	ObjFunction2 tempv1, tempv2;
	std::set<ObjFunction2> counter_table;

	//copies all the objective vectors to set that uniquely stores them
	//(in other words), stores all the objectives and filters out duplicates
	for(iterator it=begin(); it != end() ; ++it ){
		counter_table.insert(it->second);
	}
	//isws xreiastei gia mellontikh eukairia (stis n-eggrafes na dinw kapoio pososto extra)
	double hv2=0;

	std::set<ObjFunction2>::iterator itHV=counter_table.begin();

	hv2=euclidean_dist( (*itHV), reference_point);

	++itHV;
	while( itHV!= counter_table.end()){
		hv2+=euclidean_dist( (*itHV), reference_point );
		++itHV;
	}

#ifdef DEBUG
	std::cout << "HV = " << hv2 << std::endl <<std::endl;
	std::cout << "hyper volume calculation END" << std::endl;
#endif
	return hv2;
}

int Container2::activate_kick(char const  *save_directory, int const kick_limit, ObjFunction2 failedObjectives){
	/*
	 * checks the current container and counts the frequency of the same objective vector
	 * [the variable vector, might be the same, but the objective is the same - multi-modality]
	 * if any entry exists more than kik_limit times, then the function will return 1.
	 *
	 * combined with the number of consecutive non-improvements, will perform restart move!
	 *
	 */
	//todo kalutera na sunduastei etsi wste AN upervei to kicklimit KAI an den uparxoun sunexeis anakalupseis
	//tote na ginei restart tuxaia apo kapoio allo shmeio ws vase point, an auto den doulepsei gia
	// K fores, tote na kanei reduce move.

	std::ofstream kick_file;
	kick_file.open ( save_directory ,  std::ios::out  );

	int max=0;
	std::vector<int> count(size(),0);


	int counter;

	ObjFunction2 temp;
	std::map<ObjFunction2, int> occurances;
	std::set<ObjFunction2> counter_table;
	ObjFunction2 tempv1, tempv2;

	for(iterator it=begin(); it != end() ; ++it ){
		counter_table.insert(it->second);
	}

	for( std::set<ObjFunction2>::const_iterator it2=counter_table.begin(); it2 != counter_table.end() ; ++it2 ){
		counter=0;
		copy_contents(tempv2, *it2);

		for(iterator it=begin(); it != end() ; ++it ){
			copy_contents(tempv1, it->second);
			if( compare_error (it->second, *it2) && tempv1!=failedObjectives)
				//if( compare_error (tempv1, tempv2) && tempv1!=failedObjectives)
				++counter;
		}
		occurances[*it2]=counter;
	}

	max=0;
	for ( std::map<ObjFunction2, int>::iterator it3=occurances.begin(); it3!=occurances.end(); ++it3){
		if( it3->second > max )
			max=it3->second;

		kick_file << std::fixed <<  it3->first << " appears " << it3->second << "times" << std::endl;
	}

	kick_file.close();
#ifdef DEBUG
	std::cout << "active_kick" << max << std::endl;
#endif
	if(max>kick_limit)
		return 1;

	return 0;
}

int Container2::count_evaluations(ObjFunction2 given_penalty_vector){
	unsigned int number=0;

	for(Container2::const_iterator it=begin(); it!=end(); ++it)
		if(it->second!=given_penalty_vector)
			++number;

	return number;
}

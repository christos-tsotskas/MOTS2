/*
 * short_term_memory.cpp
 *
 *  Created on: Jan 19, 2014
 *      Author: ctsotskas
 */

#include "short_term_memory.h"
#include <sstream>

STM_Container2::STM_Container2(const int &input_STM_size, int const n, std::string given_name, std::string save_path) :
STM_size (input_STM_size),
local_nVar(n),
container_name(given_name),
__save_path(save_path){

	if( container_name.compare("STM")==0 )
			std::cout << "container " << container_name << "can hold up to "<< STM_size << " decision variables, each with (nVar=)" << local_nVar << " elements" << std::endl;


	//TOCHANGE
}; //default constructutor


void STM_Container2::import(Point2 const &P){
	if( !isTabu(P)){
		if(size() < STM_size )
			push_front(P);
		else{
			pop_back();
			push_front(P);
		}
	}

}

void STM_Container2::showTabu(){
	int m=0;
	std::cout << "\t STM contents" << std::endl;
	std::cout << "\t stm size=" << size() << std::endl;
	for (iterator it=begin() ; it!= end() ; ++it){
		std::cout << "\t " << m << "th : " << *it << std::endl;
		++m;
	}
}

bool STM_Container2::isTabu(const Point2 &P){
	//cout << "STM contents" << endl;
	//for(Container::iterator it=STM.begin() ; it!= STM.end() ; it++) //point==tabuPoint, end is thrown if element is not found


	for ( iterator it=begin() ; it!=end() ; ++it)
		if ( P  == *it )
			return true;
	return false;
}

void STM_Container2::save_stm_container(std::string case_qualifier, int evaluations_counter){
	//create the ./memories/STM_snapXXXXX.txt , where XXXX the current number of evaluations

	std::ostringstream dbg_message, snapshot_name, simple_name;
	std::string temp_s (__save_path+"/" + case_qualifier +  "_" + container_name + "_snap");

	snapshot_name << __save_path+"/" << case_qualifier <<  "_" << container_name << "_snap" << evaluations_counter << ".txt";
	simple_name << __save_path+"/" << case_qualifier <<  "_" << container_name  << ".txt";

	std::ofstream STMfile(snapshot_name.str().c_str()); //snapshot file
	std::ofstream simple_file(simple_name.str().c_str()); //general file

	for( iterator it= begin() ; it != end() ; ++it){
		STMfile << *it << std::endl;
		simple_file << *it << std::endl;
	}

	STMfile.close();
	simple_file.close();

#ifdef DEBUG
	dbg_message<<__save_path<<" saved!";
	std::cout << dbg_message.str() << std::endl;
#endif
}

int STM_Container2::load_stm_container(const char* load_path){
	double buffer;
	Point2 tempPoint(local_nVar,0.0);
	unsigned int i;


	std::ifstream LOADfile(load_path);


	std::cout << " STM is reading from file"<< load_path << std::endl;
	std::cout << "\t" << tempPoint.size() << " variables " << std::endl;

	while(LOADfile >> buffer){ //the loop reads through each element
		for(i=0; i<tempPoint.size()-1 ; ++i){
			//reads one element at a time (ON THE SAME LINE)
			tempPoint[i]=buffer; //feeds the above element into the vector, until the line "entry" is replicated into the vector
			LOADfile >> buffer;
		}
		tempPoint[i]=buffer;//for the last element of the current line

		//cout << j++ << ". V1=" << bufferVector << endl;

		push_back(tempPoint);

	}
	LOADfile.close();
	std::cout << "\ttotally " << size() << "STM points loaded" << std::endl;
	std::cout << "----------------------------------------------------" << std::endl;

#ifdef REPORT_CONTAINERS
	showTabu();
	cout << "STM READ" <<endl;
#endif
	return size();
}


//bool STM_Container2::isTabu(const Point2 &P){
//	//cout << "STM contents" << endl;
//	//for(Container::iterator it=STM.begin() ; it!= STM.end() ; it++) //point==tabuPoint, end is thrown if element is not found
//
//
//	for ( iterator it=begin() ; it!=end() ; ++it)
//		if ( P  == *it )
//			return true;
//	return false;
//}





//void STM_Container2::save_stm_container(const char* save_path, int evaluations_counter){
//	//create the ./memories/STM_snapXXXXX.txt , where XXXX the current number of evaluations
//
//	char* snapshot_name = new char[strlen(save_path)+sizeof(evaluations_counter)+50];
//	char* simple_name = new char[strlen(save_path)+sizeof(evaluations_counter)+50];
//
//	sprintf(snapshot_name,"%s%s%d%s",save_path,"_snap",evaluations_counter,".txt");
//	sprintf(simple_name,"%s%s",save_path,".txt");
//
//	std::ofstream STMfile(snapshot_name); //snapshot file
//	std::ofstream simple_file(simple_name); //general file
//
//	for( iterator it= begin() ; it != end() ; ++it){
//		STMfile << *it << std::endl;
//		simple_file << *it << std::endl;
//	}
//
//	STMfile.close();
//	simple_file.close();
//
//	delete[] snapshot_name;
//	delete[] simple_name;
//}

//int STM_Container2::load_stm_container(const char* load_path){
//	double buffer;
//	Point2 tempPoint(local_nVar,0.0);
//	unsigned int i;
//
//	std::ifstream LOADfile(load_path);
//
//
//	std::cout << " STM is reading from file"<< load_path << std::endl;
//	std::cout << "\t" << tempPoint.size() << " variables " << std::endl;
//
//	while(LOADfile >> buffer){ //the loop reads through each element
//		for(i=0; i<tempPoint.size()-1 ; ++i){
//			//reads one element at a time (ON THE SAME LINE)
//			tempPoint[i]=buffer; //feeds the above element into the vector, until the line "entry" is replicated into the vector
//			LOADfile >> buffer;
//		}
//		tempPoint[i]=buffer;//for the last element of the current line
//
//		//cout << j++ << ". V1=" << bufferVector << endl;
//
//		push_back(tempPoint);
//
//	}
//	LOADfile.close();
//	std::cout << "\ttotally " << size() << "STM points loaded" << std::endl;
//	std::cout << "----------------------------------------------------" << std::endl;
//
//#ifdef REPORT_CONTAINERS
//	showTabu();
//	cout << "STM READ" <<endl;
//#endif
//	return size();
//}

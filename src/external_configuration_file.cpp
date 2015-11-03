/*
 * external_configuration_file.cpp
 *
 *  Created on: Jan 20, 2014
 *      Author: ct01
 */


#include <iostream>
#include <cstring>
#include <fstream>
#include <cstdlib>

#include "external_configuration_file.h"

extern int FileExists(char const *strFilename);

unsigned int general_configuration::get_decision_variable_size(){
	return __nVar;
}

unsigned int general_configuration::get_objectives_size(){
	return __nObj;
}

std::string general_configuration::get_filename(){
	return __filename;
}

void general_configuration::list_readings(){

}

const std::string& general_configuration::getAssessment() const {
	return __assessment;
}

unsigned int general_configuration::getDiversify() const {
	return __diversify;
}

unsigned int general_configuration::getEl() const {
	return __EL;
}

const std::string& general_configuration::getFilename() const {
	return __filename;
}

unsigned int general_configuration::getIl() const {
	return __IL;
}

unsigned int general_configuration::getIntensify() const {
	return __intensify;
}

unsigned int general_configuration::getLl() const {
	return __LL;
}

const std::string& general_configuration::getLogtype() const {
	return __logtype;
}

unsigned int general_configuration::getMaximumDuplicates() const {
	return __maximum_duplicates;
}

unsigned int general_configuration::getMaximumImprovements() const {
	return __maximum_improvements;
}

unsigned int general_configuration::getSample() const {
	return __n_sample;
}

unsigned int general_configuration::getObj() const {
	return __nObj;
}

unsigned int general_configuration::getRegions() const {
	return __nRegions;
}

unsigned int general_configuration::getVar() const {
	return __nVar;
}

unsigned int general_configuration::getReduce() const {
	return __reduce;
}

unsigned int general_configuration::getSaveStep() const {
	return __save_step;
}

double general_configuration::getSs() const {
	return __SS;
}

double general_configuration::getSsrf() const {
	return __SSRF;
}

double general_configuration::getStartingPoint() const {
	return __starting_point;
}

unsigned int general_configuration::getStmSize() const {
	return __STM_size;
}

general_configuration::general_configuration(std::string filename,
		unsigned int diversify,
		unsigned int intensify,
		unsigned int reduce,
		double SS,
		double SSRF,
		unsigned int save_step,
		unsigned int n_sample,
		unsigned int nVar,
		unsigned int nObj,
		unsigned int LL, //loop limit
		unsigned int EL, //evaluations limit
		unsigned int IL, //Improvements limit , number of consecutive improvements
		unsigned int nRegions,
		unsigned int STM_size,
		std::string assessment,
		std::string logtype,
		double starting_point,
		unsigned int maximum_improvements,
		unsigned int maximum_duplicates):

						__filename(filename),
						__diversify(diversify),
						__intensify(intensify),
						__reduce(reduce),
						__SS(SS),
						__SSRF(SSRF),
						__save_step(save_step),
						__n_sample(n_sample),
						__nVar(nVar),
						__nObj(nObj),
						__LL(LL), //loop limit
						__EL(EL), //evaluations limit
						__IL(IL), //Improvements limit , number of consecutive improvements
						__nRegions(nRegions),
						__STM_size(STM_size),
						__assessment(assessment),
						__logtype(logtype),
						__starting_point(starting_point),
						__maximum_improvements(maximum_improvements),
						__maximum_duplicates(maximum_duplicates),
						__internal_configuration(1) {

	std::cout << "internal configuration was set" << std::endl;
	show_configurations();
}

void general_configuration::show_configurations() const{
	std::cout << std::endl << "Tabu Search Options" << std::endl;

	std::cout << " 1.diversify " << __diversify << std::endl;
	std::cout << " 2.intensify " <<  __intensify<< std::endl;
	std::cout << " 3.reduce " <<  __reduce<< std::endl;

	std::cout << " 4.SS " <<  __SS<< std::endl;
	if( __SS==0 )
		std::cout << "\t Start step(s) is assigned by user in start_step.txt file!" << std::endl;
	else
		std::cout << "\t Start step(s) (" << __SS <<")  is uniform for all variables" << std::endl;

	std::cout << " 5.SSRF " <<  __SSRF<< std::endl;
	std::cout << " 6.save_step " <<  __save_step<< std::endl;
	std::cout << " 7.n_sample " <<  __n_sample<< std::endl;
	std::cout << " 8.nVar " <<  __nVar<< std::endl;
	std::cout << " 9.nObj " <<  __nObj<< std::endl;
	std::cout << " 10.n_of_loops " <<  __LL<< std::endl;
	std::cout << " 11.n_of_evaluations " <<  __EL<< std::endl;
	std::cout << " 12.n_of_consecutive_improvements " <<  __IL<< std::endl;
	std::cout << " 13.assessment " <<  __assessment<< std::endl;
	std::cout << " 14.nRegions " <<  __nRegions<< std::endl;
	std::cout << " 15.STM_size " <<  __STM_size<< std::endl;
	std::cout << " 16.LogType " <<  __logtype<< std::endl;

	std::cout << " 17.Starting point " <<  __starting_point << std::endl;
	if(__internal_configuration==1){
		std::cout << "\t WARNING!!! Point is set internally!" << std::endl;
	}else{
		if( __starting_point==0)
			std::cout << "\t Starting point is assigned randomly!" << std::endl;
		else
			std::cout << "\t Starting point is by user from the file datum_design_vector.txt" << std::endl;

	}

	if(-1<__starting_point && __starting_point<0){
		std::cout << "\t\tInitial sampling was activated, too!" << std::endl;
	}

	std::cout << " 18.maximum_improvements " <<  __maximum_improvements<< std::endl;
	std::cout << " 19.maximum_duplicates " <<  __maximum_duplicates<< std::endl;
}

general_configuration::general_configuration(std::string filename):
																__filename(filename),
																__diversify(0),
																__intensify(0),
																__reduce(0),
																__SS(0.0),
																__SSRF(0.0),
																__save_step(0),
																__n_sample(0),
																__nVar(0),
																__nObj(0),
																__LL(0),
																__EL(0),
																__IL(0), // available values HV / HB
																__nRegions(0),
																__STM_size(0),
																__assessment(""),
																__logtype(""), // available values full/quick
																__starting_point(0.0),
																__maximum_improvements(0),
																__maximum_duplicates(0),
																__internal_configuration(0)
{

	std::cout << "external configuration file (" << filename << ") ...";
	if( FileExists(filename.c_str()) )
		std::cout << " OK!"<< std::endl;
	else{
		std::cout << "NOT FOUND!!!" << std::endl;
		exit(-1000);
	}

	std::ifstream configuration_file(filename.c_str());

	if(   configuration_file.is_open() ){
		std::cout << "configuration file found! " << std::endl;

		configuration_file >> __diversify; //1
		configuration_file >> __intensify; //2
		configuration_file >> __reduce ; //3
		configuration_file >> __SS ; //4
		configuration_file >> __SSRF ; //5
		configuration_file >> __save_step ; //6
		configuration_file >> __n_sample ; //7
		configuration_file >> __nVar ; //8
		configuration_file >> __nObj ; //9
		configuration_file >> __LL ; //10
		configuration_file >> __EL ; //11
		configuration_file >> __IL ; //12  : available values HV / HB
		configuration_file >> __assessment ; //13
		configuration_file >> __nRegions ; //14
		configuration_file >> __STM_size ; //15
		configuration_file >> __logtype ; //16 -  available values full/quick
		configuration_file >> __starting_point ; //17
		configuration_file >> __maximum_improvements; //18
		configuration_file >> __maximum_duplicates; //19

		configuration_file.close();
		std::cout << "configuration file read! " << std::endl;
		show_configurations();
	}else{
		std::cout << "NO configuration file found!" <<std::endl;
	}
};


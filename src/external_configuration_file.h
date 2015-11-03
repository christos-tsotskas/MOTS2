/*
 * external_configuration_file.h
 *
 *  Created on: Jan 20, 2014
 *      Author: ct01
 */

#ifndef EXTERNAL_CONFIGURATION_FILE_H_
#define EXTERNAL_CONFIGURATION_FILE_H_

#include <cstring>

/**
 * Contains the external settings, as read on the configuration file prescribed
 */

class general_configuration{
	std::string __filename;
	unsigned int __diversify;
	unsigned int __intensify;
	unsigned int __reduce;
	double __SS;
	double __SSRF;
	unsigned int __save_step;
	unsigned int __n_sample;
	unsigned int __nVar;
	unsigned int __nObj;
	unsigned int __LL; //loop limit
	unsigned int __EL; //evaluations limit
	unsigned int __IL; //Improvements limit , number of consecutive improvements
	unsigned int __nRegions;
	unsigned int __STM_size;
	std::string __assessment;
	std::string __logtype;
	double __starting_point;
	unsigned int __maximum_improvements;
	unsigned int __maximum_duplicates;
	unsigned int __internal_configuration;

public:
	general_configuration(std::string filename,
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
	unsigned int maximum_duplicates);

	general_configuration(std::string filename);

	void show_configurations() const;

	unsigned int get_decision_variable_size();

	unsigned int get_objectives_size();

	std::string get_filename();

	void list_readings();

	const std::string& getAssessment() const ;

	unsigned int getDiversify() const ;

	unsigned int getEl() const ;

	const std::string& getFilename() const ;

	unsigned int getIl() const ;

	unsigned int getIntensify() const ;

	unsigned int getLl() const ;

	const std::string& getLogtype() const ;

	unsigned int getMaximumDuplicates() const ;

	unsigned int getMaximumImprovements() const ;

	unsigned int getSample() const ;

	unsigned int getObj() const ;

	unsigned int getRegions() const ;

	unsigned int getVar() const ;

	unsigned int getReduce() const ;

	unsigned int getSaveStep() const ;

	double getSs() const ;

	double getSsrf() const ;

	double getStartingPoint() const ;

	unsigned int getStmSize() const;
};

#endif /* EXTERNAL_CONFIGURATION_FILE_H_ */

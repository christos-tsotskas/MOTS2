/*
 * configuration_settings.h
 *
 *  Created on: Jan 19, 2014
 *      Author: ctsotskas
 */

#ifndef CONFIGURATION_SETTINGS_H_
#define CONFIGURATION_SETTINGS_H_

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
#include "external_configuration_file.h"


class Overall_Optimisation_Configuration_Settings{
	//todo make configurations a map and then load the settings
	std::string __case_name;
	general_configuration __GeneralConfigurationFile;
	ObjFunction2 reference_point;
	ObjFunction2 penalty_point;
	Point2 lower_bound;
	Point2 upper_bound;
	Point2 starting_point;
	Point2 current_step;

	ObjFunction2  read_vector_from_file(char const *filename, char const *file_description, const int elements,const int necessary);

public:
	Overall_Optimisation_Configuration_Settings(std::string case_name,
			general_configuration &external_configuration_file,
			std::string ReferencePointFilename,
			std::string PenaltyPointFilename,
			std::string LowerBoundFilename,
			std::string UpperBoundFilename,
			std::string StartingPointFilename,
			std::string CurrentStepFilename);

	Overall_Optimisation_Configuration_Settings(std::string case_name,
			general_configuration& internal_configuration,
			ObjFunction2 custom_reference_point,
				ObjFunction2 custom_penalty_point,
				Point2 custom_lower_bound,
				Point2 custom_upper_bound,
				Point2 custom_starting_point,
				Point2 custom_current_step);

	~Overall_Optimisation_Configuration_Settings();
	const ObjFunction2& get_reference_point() const;
	const Point2& get_starting_point() const;
	const Point2& get_current_step() const;
	const ObjFunction2& get_penalty_point() const;
	const Point2& get_lower_bound() const;
	const Point2& get_upper_bound() const;
	const general_configuration& getExternalConfigurationFile() const;
	Point2 generateRandomPointWithinBounds() const;
	const std::string& getCaseName() const;
};

#endif /* CONFIGURATION_SETTINGS_H_ */

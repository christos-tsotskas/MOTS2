/*
 * container.h
 *
 *  Created on: Jan 19, 2014
 *      Author: ctsotskas
 */

#ifndef CONTAINER_H_
#define CONTAINER_H_

#include <map>
#include <cstring>
#include "global_defines.h"

class Container2 : public std::map< Point2, ObjFunction2 >{
	unsigned int local_nVar;
	unsigned int local_nObj;
	std::string container_name;
	std::string __save_path;
	typedef std::pair<Point2,ObjFunction2> tuple;
	int dominates(const ObjFunction2 &O1, const ObjFunction2 &O2);
	void displayContainer ();
	int report_optima_plot(char const *save_path);
	int report_optima_current(char const *save_path);
	double calculate_hypervolume2D(const ObjFunction2 &reference_point);
public:
	void add_point(const entry __entry);
	int load_container(char const *load_path);
	double calculate_quality_indicator(const ObjFunction2 &reference_point);
	int activate_kick(char const  *save_directory, int const kick_limit, ObjFunction2 failedObjectives);
	int count_evaluations(ObjFunction2 given_penalty_vector);
	Point2 check_memory();
	int removeDominatedPoints();
	int addIfNotDominated(const Point2 &V1,const ObjFunction2 &O1);
	Point2 selectRandom();
	Point2 find_extreme(int objective);
	int save_container_snapshot(std::string case_qualifier, int evaluations_counter);
	Container2(unsigned int n, unsigned int m, std::string given_name, std::string filename);
	~Container2();
};


#endif /* CONTAINER_H_ */

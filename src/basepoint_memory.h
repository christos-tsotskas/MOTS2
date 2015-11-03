/*
 * basepoint_memory.h
 *
 *  Created on: Jan 25, 2014
 *      Author: christos
 */

#ifndef BASEPOINT_MEMORY_H_
#define BASEPOINT_MEMORY_H_

#include <deque>
#include <cstring>

#include "global_defines.h"

class basepoint_entry{
	Point2 __basepoint;
	std::string __function;
public:
	basepoint_entry(Point2 basepoint,  std::string function): __basepoint(basepoint), __function(function) {}
	const Point2& get__basepoint() const;
	const std::string& getFunction() const;
};


class basepoint_memory: public std::deque<basepoint_entry>{

public:
	basepoint_memory();
	void insert_point(const Point2 &P, std::string const entry_function);
	const Point2& get_current_base_point() const;
	const Point2& get_previous_base_point() const;
	void show_memory();
	void save_memory(std::string __save_path, std::string case_qualifier, std::string container_name,  int evaluations_counter);
	~basepoint_memory();

};



#endif /* BASEPOINT_MEMORY_H_ */

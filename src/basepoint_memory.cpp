/*
 * basepoint_memory.cpp
 *
 *  Created on: Jan 25, 2014
 *      Author: christos
 */

#include <sstream>
#include <fstream>

#include "basepoint_memory.h"


basepoint_memory::basepoint_memory(){
	std::cout << "BasePoint Memory created!" << std:: endl;
}

void basepoint_memory::insert_point(const Point2 &P, std::string const entry_function){
	basepoint_entry new_entry(P,entry_function);
	push_back( new_entry );
}

const Point2& basepoint_memory::get_current_base_point() const{
	return back().get__basepoint();
}

const Point2& basepoint_memory::get_previous_base_point() const{
	return (--(--end()))->get__basepoint();
}

basepoint_memory::~basepoint_memory(){
	std::cout << "BasePoint Memory deleted!" << std:: endl;
}

void basepoint_memory::show_memory(){
	int i=0;
	for( basepoint_memory::const_iterator it=begin() ; it!=end(); ++it){
		std::cout << i << "\t" << it->get__basepoint() << "\t" << it->getFunction() << std::endl;
		++i;
	}

}

void basepoint_memory::save_memory(std::string __save_path, std::string case_qualifier, std::string  container_name,  int evaluations_counter){
	std::ostringstream dbg_message, snapshot_name, simple_name;
	std::string temp_s (__save_path+"/" + case_qualifier +  "_" + container_name + "_snap");

	snapshot_name << __save_path+"/" << case_qualifier <<  "_" << container_name << "_snap" << evaluations_counter << ".txt";
	simple_name << __save_path+"/" << case_qualifier <<  "_" << container_name  << ".txt";

	std::ofstream BASEfile(snapshot_name.str().c_str());
	std::ofstream simple_file(simple_name.str().c_str());


//	BPfile << "#loop" << " "<< "decision variables" << "\t" << "objectives" << " reg=" << "region" << " move=" << "way_it_entered" << endl;
	int i=0;
	for( basepoint_memory::const_iterator it=begin() ; it!=end(); ++it){
		BASEfile << i << "\t" << it->get__basepoint() << "\t" << it->getFunction() << std::endl;
		simple_file << i << "\t" << it->get__basepoint() << "\t" << it->getFunction() << std::endl;
		++i;
	}



	BASEfile.close();
	simple_file.close();

#ifdef DEBUG
	dbg_message<<__save_path<<" saved!";
	std::cout << dbg_message.str() << std::endl;
#endif


}

const Point2& basepoint_entry::get__basepoint() const{
	return __basepoint;
}

const std::string& basepoint_entry::getFunction() const {
	return __function;
}

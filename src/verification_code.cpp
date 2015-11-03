/*
 * verification_code.cpp
 *
 *  Created on: Jan 24, 2014
 *      Author: ct01
 */

#include <fstream>

#include "verification_code.h"

extern int FileExists(char const *strFilename);

ObjFunction2 simple_read_vector_from_file(std::string filename){

	std::cout  << "user requested file ..."<< filename<< std::endl;
	std::ifstream FILE(filename.c_str());
	std::vector<double> buff_vec;

	if( FileExists(filename.c_str()) ){
		std::cout << " OK!"<< std::endl;

		if(FILE.is_open()){
			double buff=0.0;

			while(FILE>>buff)
				buff_vec.push_back(buff);

			FILE.close();
			std::cout << "\t"<< filename<< " read: " << buff_vec << std::endl;
		}else{
			std::cout << " The file " << filename <<"  was not found in root directory" << std::endl;
		}
	}else{
		std::cout << "NOT FOUND!!! (looking for "<< filename << " within the current directory)"<< std::endl;
		std::cout << "(user asked for file)" << std::endl;
		exit(-1002);
	}

	return buff_vec;
}

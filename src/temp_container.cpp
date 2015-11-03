/*
 * temp_container.cpp
 *
 *  Created on: Jan 19, 2014
 *      Author: ctsotskas
 */


#include "temp_container.h"

temp_Container2::temp_Container2 (){
	srand((unsigned)(time(0)+4682));
} //default constructor


//Algorithm 6e
Point2 temp_Container2::selectRandom(){//temp_container

	iterator it;
	int jump=rand() % size();
#ifdef DEBUG
	std::cout << "[tmp2_C] among" << size() << " points I have choosen the " << jump+1 << "th element" << std::endl;
#endif
	advance(it=begin(), jump);
	//advance(it=begin(),  RandomNumber(0, size())  );
	return *it;
}

void temp_Container2::displayContent(){
	std::cout << "\t show temp_container" << std::endl;
	std::cout << "\t its size is=" << size() << std::endl;

	iterator itt=begin();
	int vector_size = itt->size();
	int k=1;
	for( iterator itC=begin()  ; itC != end() ; ++itC){
		std::cout << "\t" << k <<":\t" << *itC << std::endl;
		k++;

	}
	std::cout << "\t its size is=" << size() << std::endl;
	std::cout << "\t entry size is=" << vector_size  << std::endl;
	std::cout << "end of container / bestpoints" << std::endl;
}

int temp_Container2::checkduplicate(){
	int total=0;
	iterator previous;
	iterator itC=begin();
	previous = itC;
	++itC;
	for(    ; itC != end() ; ++itC){
		if(*itC == *previous)
			total++;

		++previous;
	}

	return total;
}

#include <iostream>

void testC11(){
	int array [] = {5, 6, 12};
	for (auto i : array)
		std::cout << "element: " << i << std::endl;

}

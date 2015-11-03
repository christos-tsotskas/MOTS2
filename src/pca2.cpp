

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>




static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)



/**  Error handler  **************************************************/
void erhand(const std::string err_msg){
	/* Error handler */
	std::cout << "Run-time error:\n";
	std::cout << err_msg.c_str();
	std::cout << "Exiting to system.\n";
	exit(1);
}

float pythag(float a,float b){
	float absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}



////template void dumpVector<int>(std::vector<int> v, std::string sep);
//template vector<double>& PCA<double>(vector< vector<double> > &data, const int n, const int m)

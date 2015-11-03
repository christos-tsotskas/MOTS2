#include <iostream>
#include <map>
#include <vector>
#include <cmath>

#include "cutil_inline.h"
#include "container.h"
#include "global_defines.h"

#ifdef GPU_evaluation

//#define DEBUG

__global__ void calc_objectives(const int entries, const double* variables, double* objectives, const int nVar, const int nObj);
__global__ void internal_print_objectives(const int entries, const double* variables, const double* objectives, const int nVar, const int nObj);
__global__ void internal_print_variables(const int entries, const double* variables, const double* objectives, const int nVar, const int nObj);

class CUDA_evaluation{
//	Container2 temp_list;

	double* decision_variables_h;
	double* objectives_h;

	double* dev_variables;
  	double* dev_objective;

	unsigned int entries;
	unsigned int nVar;
	unsigned int nObj;

	int nVar_total_size;
	int nObj_total_size;

	int number_of_threads;
	int number_of_blocks;


	void strip_to_1D(const Container2& in_map, double* out_decision_variables, double* out_objectives);
	void display_data(const double* dataset,const int entries, const int col);
	void display_data(const std::map< std::vector<double>, std::vector<double> > &in_map);
	double calc_obj1(const double &a, const double &b);
	double calc_obj2(const double &a, const double &b);

	void run_thread(unsigned int thread,
		unsigned int block,
		const unsigned int nEntries,
		const double* out_decision_variables,
		const unsigned int nVar,
		double* out_objectives,
		const unsigned int nObj);

	void allocate_memory();
	void send_to_device();

	void print_CUDA();

	void receive_from_GPU();


public:
	CUDA_evaluation(int selected_threads ,const Container2& __temp_list);
	void run();
	void update_container_objectives(Container2& __temp_list);
	~CUDA_evaluation();

	void check_results();
};

#endif //GPU_evaluation

//============================================================================
// Name        : interface_for_simple_data_structures.cpp
// Author      : Christos Tsotskas
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <map>
#include <vector>
#include <cmath>



#include "test.cuh"


#ifdef GPU_evaluation


void CUDA_evaluation::strip_to_1D(const Container2 &in_map, double* out_decision_variables, double* out_objectives){
#define var_index_map(entry_index,var_index) (entry_index*nVar+var_index)
#define obj_index_map(entry_index,obj_index) (entry_index*nObj+obj_index)

	const unsigned int nEntries=in_map.size() ;
	const unsigned int nVar=(in_map.begin())->first.size() ;
	const unsigned int nObj=(in_map.begin())->second.size() ;

	unsigned int entry_index=0;
	for(std::map< std::vector<double>, std::vector<double> >::const_iterator  it=in_map.begin(); it!=in_map.end(); ++it){
		for(unsigned int var_index=0 ; var_index<nVar; ++var_index)
			//			out_decision_variables[entry_index*nVar+var_index]=it->first[var_index];
			out_decision_variables[var_index_map(entry_index,var_index)]=it->first[var_index];

		for(unsigned int obj_index=0 ; obj_index<nObj; ++obj_index)
			//			out_objectives[entry_index*nObj+obj_index]=it->second[obj_index];
			out_objectives[obj_index_map(entry_index,obj_index)]=it->second[obj_index];

		++entry_index;
	}

	//	int a;
	//	std::cout << "type int" << std::endl;
	//	std:: cin >> a;
#undef var_index_map
#undef obj_index_map
}

void CUDA_evaluation::display_data(const double* dataset,const int entries, const int col){
	for(int i=0; i<entries; ++i){
		std::cout << i  << ": ";
		for(int j=0; j<col; ++j){
			std::cout << dataset[i*col+j] << " ";
		}
		std::cout  << std::endl;
	}
}


void CUDA_evaluation::display_data(const std::map< std::vector<double>, std::vector<double> > &in_map){
	int entry=0;
	for(std::map< std::vector<double>, std::vector<double> >::const_iterator  it=in_map.begin(); it!=in_map.end(); ++it){
		std::cout << entry  << ": ";

		for(unsigned int var_index=0 ; var_index<it->first.size(); ++var_index)
			std::cout << it->first[var_index] << " ";

		std::cout << "===>";
		for(unsigned int obj_index=0 ; obj_index<it->second.size(); ++obj_index)
			std::cout << it->second[obj_index] << " ";


		std::cout  << std::endl;
		++entry;
	}
}

double CUDA_evaluation::calc_obj1(const double &a, const double &b){
	return a+b;
}

double CUDA_evaluation::calc_obj2(const double& a, const double& b){
	return a+b;
}
//The design of the optimiser was based on a different basis, compared to the traditional approach.
//More specifically, there are two key differences. First, the algorithm was designed to match the capabilities of the hardware
// and it can scale by adding a more powerful in terms of specifications hardware. In turn, the main principle is to fully utilise
// the available computational resources as much as possible.
//Second, it aims to tackle multi-objective
//problems and with some trivial tweaking single-objective problems can be dealt, too. Consequently, the internal data structures
//and algorithms were selected to implement the Pareto-optimality satisfactorily.

// Since the hardware is not mature enough, it is expected not to be immediately applicable. However, as the technology becomes
//more matures, this approach will tend to be more realistic.


//extra, following the locality of data, there could be re-usability of some components, if the blocks seems similar

//as an extension, the code, could separate the list of samples on blocks, based on similarity patterns (and would also inform the GPU about the similarity).
//This is done to avoid any unexpected searchs and branches at the execution.

//as an extension to the structure of the optimiser, the pattern move should be altered to be more collective, not sequential
//as it is currently.

//as an extentsion, currently all the objectives of a single design are evaluated by a single thread. It worths to investigate how the code performs when each objective is evaluated on a separate thread (which could simplify/alleviate the computational requirements)

//as an extension, the time required to send and receive data from GPU can be improved by using the stream feature from CUDA. Although this was not used here, it is expected to shortent the exchange time when the number of variables grows significantly.

//at conclusions, the presented approch for GPU implementation was selected so as to allow for future implementations to use GPUs' spatial locality more efficiently.

//idea: map an entry to a thread, the total number of threads launched equals the sampling. however, it is not guaranteed that
//this number will always be the same. The GPU can be used in many different ways. Here, the selected approach is to do every evaluation
//on each thread, or simply it serves as an evaluation-factory. Each evaluation is independent from the other, however they are evaluated
// in batches.
//(so it can just scale up as the more resources are added, one requirement is that the problem should ask/produce so many points to be evaluated.).
//The CPU part will do the high level and complex part of the code, the evaluations
void CUDA_evaluation::run_thread(unsigned int thread,
		unsigned int block,
		const unsigned int nEntries,
		const double* out_decision_variables,
		const unsigned int nVar,
		double* out_objectives,
		const unsigned int nObj){

#define var_index_map(entry_index,var_index) (entry_index*nVar+var_index)
#define obj_index_map(entry_index,obj_index) (entry_index*nObj+obj_index)

	const unsigned int Nthreads=2;
	const unsigned int gid=block* Nthreads+thread;

	if(gid<nEntries){

		out_objectives[obj_index_map(gid,0)]=calc_obj1(out_decision_variables[var_index_map(gid,0)], out_decision_variables[var_index_map(gid,1)]);

		out_objectives[obj_index_map(gid,1)]=calc_obj2(out_decision_variables[var_index_map(gid,1)], out_decision_variables[var_index_map(gid,2)]);

	}
	std::cout<<"gid:" << gid << " calculated:" << out_objectives[obj_index_map(gid,0)] << " and " << out_objectives[obj_index_map(gid,1)] << std::endl;

#undef var_index_map
#undef obj_index_map
}

void CUDA_evaluation::run() {

#ifdef DEBUG
	std::cout << "testing_kernel" << std::endl;
#endif //DEBUG	
	calc_objectives<<<number_of_blocks, number_of_threads>>>(entries, dev_variables, dev_objective, nVar, nObj);
	cudaDeviceSynchronize();
	receive_from_GPU();
//	link_objectives_to_map();

}

void CUDA_evaluation::allocate_memory(){
	cudaMalloc( (void**) &dev_variables ,  nVar_total_size);
	cudaMalloc( (void**) &dev_objective ,  nObj_total_size);
}

void CUDA_evaluation::send_to_device(){
	cudaMemcpy( dev_variables, decision_variables_h, nVar_total_size , cudaMemcpyHostToDevice );
	//cudaMemcpy( dev_objective, objectives_h, nObj_total_size , cudaMemcpyHostToDevice );
}




CUDA_evaluation::CUDA_evaluation(int selected_threads, const Container2& temp_list)
				{

	entries= temp_list.size();
	nVar=(temp_list.begin())->first.size();
	nObj=(temp_list.begin())->second.size();

	decision_variables_h= new double[ entries * nVar ];
	objectives_h= new double[ entries * nObj ];

	number_of_threads=selected_threads;
	number_of_blocks= (int) std::ceil (1.0*entries/number_of_threads*1.0);

	std::cout << "received:" << entries << " entries" << std::endl;
	std::cout << "number_of_threads:" << number_of_threads << std::endl;
	std::cout << "number_of_blocks:" << number_of_blocks << std::endl;
#ifdef DEBUG
	temp_list.displayContainer();
#endif //DEBUG
	strip_to_1D(temp_list, decision_variables_h, objectives_h);

	nVar_total_size=entries *nVar*sizeof(double);
	nObj_total_size=entries *nObj*sizeof(double);

#ifdef DEBUG
	std::cout << "variables' size "<< entries << "x" << nVar << std::endl;
	std::cout << "objectives' size "<< entries << "x" << nObj << std::endl;

	std::cout << "received variables" << std::endl;
	display_data(decision_variables_h, entries, nVar);

	std::cout << "received objectives" << std::endl;
	display_data(objectives_h, entries, nObj);
#endif //DEBUG

	allocate_memory();
	send_to_device();
#ifdef DEBUG
	display_data(temp_list);
#endif //DEBUG
}

__global__ void internal_print_variables(const int entries, const double* variables, const double* objectives, const int nVar, const int nObj){

#define var_index_map(entry_index,var_index) (entry_index*nVar+var_index)
#define obj_index_map(entry_index,obj_index) (entry_index*nObj+obj_index)

	int gid=threadIdx.x + blockIdx.x*blockDim.x;


	if(gid<entries){
		printf("cuda_gid %d's variables: %f %f %f\n", gid, variables[var_index_map(gid,0)], variables[var_index_map(gid,1)], variables[var_index_map(gid,2)]);
	}
	__syncthreads();
#undef var_index_map
#undef obj_index_map
}


__global__ void internal_print_objectives(const int entries, const double* variables, const double* objectives, const int nVar, const int nObj){

#define var_index_map(entry_index,var_index) (entry_index*nVar+var_index)
#define obj_index_map(entry_index,obj_index) (entry_index*nObj+obj_index)

	int gid=threadIdx.x + blockIdx.x*blockDim.x;
	if(gid<entries){
		printf("cuda_gid %d's objectives: %f %f\n", gid, objectives[obj_index_map(gid,0)], objectives[obj_index_map(gid,1)] );
	}
	__syncthreads();
#undef var_index_map
#undef obj_index_map
}

__device__ double objective2(const double* variables, const int nVar, const double f1, const int entry) {
#define var_index_map(entry_index,var_index) (entry_index*nVar+var_index)
	double g = 0.0;
	double h = 0.0;

	for (int i = 1; i < nVar; i++){
		g += variables[var_index_map(entry,i)];
	}
	g = 1.0 + 9.0 * g / (nVar-1);
	h = 1.0 - std::pow(f1 / g, 2);

#undef var_index_map
	return g * h;
}

__global__ void calc_objectives(const int entries, const double* variables, double* objectives, const int nVar, const int nObj){

#define var_index_map(entry_index,var_index) (entry_index*nVar+var_index)
#define obj_index_map(entry_index,obj_index) (entry_index*nObj+obj_index)

	int gid=threadIdx.x + blockIdx.x*blockDim.x;
	if(gid<entries){

		const double f1 = variables[var_index_map(gid,0)];

		objectives[obj_index_map(gid,0)]=f1;
		objectives[obj_index_map(gid,1)]= objective2(variables, nVar, f1, gid);
#ifdef DEBUG
		printf("gid:%d, calculated obj1: %f and obj2:%f \n",gid, objectives[obj_index_map(gid,0)], objectives[obj_index_map(gid,1)] );
#endif //DEBUG
	}

	__syncthreads();
#undef var_index_map
#undef obj_index_map	
}


void CUDA_evaluation::print_CUDA(){




	internal_print_variables<<<number_of_blocks, number_of_threads>>>(entries, dev_variables, dev_objective, nVar, nObj);
	internal_print_objectives<<<number_of_blocks, number_of_threads>>>(entries, dev_variables, dev_objective, nVar, nObj);	

}
void CUDA_evaluation::receive_from_GPU(){
	cudaMemcpy( objectives_h, dev_objective, nObj_total_size , cudaMemcpyDeviceToHost );
}


void CUDA_evaluation::check_results(){

	std::vector<double> temp_check_objective(nObj,0.0 );
	/*for(std::map< std::vector<double>, std::vector<double> >::iterator it=temp_list.begin(); it!=temp_list.end(); ++it ){

		for(int obj=0; obj<nObj; ++obj)
			for(int var=0; var<nVar; ++var)
				temp_check_objective[obj]=2;
	}*/
}

void CUDA_evaluation::update_container_objectives(Container2& __temp_list){

	int entry=0;
	ObjFunction2 temp_buffer(nObj,5555.0);
	for(Container2::iterator  it=__temp_list.begin(); it!=__temp_list.end(); ++it){
		for(unsigned int obj_index=0 ; obj_index<nObj; ++obj_index){
			temp_buffer[obj_index]=objectives_h[entry];
			++entry;
		}
		it->second=temp_buffer;
	}



}


//void CUDA_evaluation::link_objectives_to_map(){
//
//	cudaDeviceSynchronize();
//	int entry=0;
//	ObjFunction2 temp_buffer(nObj,5555.0);
//	for(Container2::iterator  it=temp_list.begin(); it!=temp_list.end(); ++it){
//		for(unsigned int obj_index=0 ; obj_index<nObj; ++obj_index){
//			temp_buffer[obj_index]=objectives_h[entry];
//			++entry;
//		}
//		it->second=temp_buffer;
//	}
//	cudaDeviceSynchronize();
//	std::cout << "within the mapping" << std::endl;
//	temp_list.displayContainer();
//}


CUDA_evaluation::~CUDA_evaluation(){
#ifdef DEBUG
	print_CUDA();
#endif //DEBUG

#ifdef DEBUG
	std::cout << "finally, variables" << std::endl;
	display_data(decision_variables_h, entries, nVar);

	std::cout << "finally, objectives" << std::endl;
	display_data(objectives_h, entries, nObj);
#endif //DEBUG
	cutilSafeCall(cudaFree(dev_variables));
	cutilSafeCall(cudaFree(dev_objective));



	free(decision_variables_h);
	free(objectives_h);

	std::cout  << "CUDA evaluations has just finished and evaluated " << entries<<" solutions. Memories were freed." << std::endl;
}



#endif //GPU_evaluation



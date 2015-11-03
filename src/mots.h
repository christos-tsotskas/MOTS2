//============================================================================
// Name        : mots.h
// Author      : Christos Tsotskas
// Version     :
// Copyright   : L-GPL
// Description : MOTS in C++, Ansi-style
//============================================================================

/* This file is part of MOTS_2.

    MOTS_2 is free software: you can redistribute it and/or modify
    it under the terms of the  GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MOTS_2 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MOTS_2.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2011 Christos Tsotskas
 */
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
#include "container.h"
#include "short_term_memory.h"
#include "long_term_memory.h"
#include "objective_function.h"
#include "configuration_settings.h"
#include "basepoint_memory.h"

//define the platform-specific part
#if defined(_WIN32) || defined(_WIN64) || defined(__WIN32) || defined(WINDOWS) && !defined(__CYGWIN__)
#include <Windows.h> //for windows!
void windows_pwd();
#define READ_TIME 1000 //Milliseconds for Windows!!
#endif //WINDOWS_RELEASE

#ifdef _MSC_VER //visual studio
#include "vld.h" //profiling with visual leak detector
#endif

#if defined(__gnu_linux__)   || defined(__linux__)
void linux_pwd();
#define READ_TIME 100000 // 3000000 in microseconds!  for linux
#endif //UNIX_RELEASE


//////////////////////////////////////////////////////////////////
///// library declarations, DO NOT EDIT!!!!!
//////////////////////////////////////////////////////////////////


extern void initialise_vector(Point2 &P, double value);
extern void initialise_objective(ObjFunction2 &P, double value);
extern void correct_accuracy(Point2 &P);

template <typename T> std::ostream& operator<<(std::ostream &os, const std::vector<T>& v);
//template <typename T, size_t S> std::ostream& operator<<(std::ostream &os, const  boost::array<T, S>& A);

extern Point2 operator- (const Point2 &v1, const Point2 &v2);
extern Point2 operator+ (const Point2 &v1, const Point2 &v2);
extern bool operator==(const Point2 &P1, const Point2 &P2);
extern bool operator>(const Point2 &P1, const Point2 &P2);
extern bool operator<(const Point2 &P1, const Point2 &P2);


template<typename T,  size_t S > void copy_contents(std::vector<T> &P, std::vector<T> const &A);
template<typename T,  size_t S > std::vector<T> vector_power(std::vector<T> &P, const double &n);
template<typename T,  size_t S > std::vector<T> vector_power(std::vector<T> const &P, const double &n);
template<typename T,  size_t S > double euclidean_dist(std::vector<T> const &current_point , std::vector<T> const &reference_point);
extern int FileExists(char const *strFilename);
template<typename T> int compare_error ( std::vector<T> const &v1,  std::vector<T> const &v2);

void ID();
double RD();

//template <typename T>
//T RandomNumber2(T min, T max);

class TabuSearch{
	//the order is critical, as the items are initialised
	//according to the order of declaration!!!
	//public:
#ifdef GATAC
	char *name; //represents the points
	char *name2; //represents the objective function values
	ofstream toBeEvaluated;
	deque<Point> temp_deque;
#endif
	std::string __case_name;
	const Overall_Optimisation_Configuration_Settings __DefaultConfigurationSettings;
	int restart_flag;
	unsigned int diversify;
	unsigned int intensify;
	unsigned int reduce;
	double SS;
	double SSRF;
	unsigned int save_step;
	unsigned int n_sample;
	unsigned int nVar;
	unsigned int nObj;
	unsigned int LL; //loop limit
	unsigned int EL; //evaluations limit
	unsigned int IL; //Improvements limit , number of consecutive improvements
	unsigned int nRegions;
	unsigned int STM_size;
	std::string assessment;
	std::string logtype;
	double starting_point;
	unsigned int maximum_improvements;
	unsigned int maximum_duplicates;

	ObjectiveFunction TS_ObjFunc;

	Point2 datumPnt;
	Point2 InitialStep; // vector opou apo8hkeyetai to arxiko step
	Point2 CurrentStep; //vector opou apo8hkeyetai to trexon step!
	ObjFunction2 reference_point;


	//
	// STM - Short Term Memory  - considered tabu - search doesn't revisit these points
	// MTM - Medium Term Memory - optimal or near optimal points are stored - used for intensification - focus on "probably" better solutions
	// LTM - Large Term Memory - logs visited areas  - used for diversification -> pulling the search onto unvisited areas
	// IM - Intensive Memory
	//

	Container2 MTM, IM, HISTORY;
	//HISTORY contains all the evaluated design vectors, if it violates contraints, then the OBJ value is set to a very large value
	//[at least larger than the worst known]

	STM_Container2 STM;
	//	std::deque<Point2> Points; //Base Points
	//	std::deque<Point2>::iterator iIter; //counter gia ton oliko arithmo of the optimisation steps performed, links with Points
	basepoint_memory __BasepointMemory;


	LTM_Container2Basic2<double> LTM;

	unsigned int iLocal; // counter for n_regions
	//This is achieved by dividing each design variable into n_regions regions and counting the number of solutions evaluated in those regions

	unsigned int cur_loops;
	int previous_cur_loops, previous_evaluations ;
	int nextMoveHookeJeeves;
	unsigned int diversification, intensification, reduction;


	typedef std::pair<Point2, const std::string> BP_tuple;
	std::deque<  BP_tuple > BP_tracer;

	double current_trade_off_quality_indicator;
	double previous_trade_off_quality_indicator;
	unsigned int trade_off_consequetive_improvements;

	std::ofstream diversify_file;
	std::ofstream evals_file;
	std::ofstream hypervolume_file;
	std::ofstream i_local_file;
	std::ofstream im_size_file;
	std::ofstream intensify_file;
	std::ofstream step_size_file;
	std::ofstream reduction_file;
	std::ofstream basePoint_file;
	std::ofstream update_memories_file;
	std::ofstream quick_file;


	unsigned int calculate_y_range_for_monitor_plots() const;
	unsigned int calculate_x_range_for_monitor_plots() const;
	const std::vector<int>& calculate_y_range_for_trade_off_plots() const;
	const std::vector<int>& calculate_x_range_for_trade_off_plots() const;
	void check_memories();
	void close_external_monitor_files();
	void create_plot_scripts();
	Point2  decrease(Point2 const &input, unsigned int const index);
	Point2 DiversifyMove2();
	void displayDeque(  std::deque<Point2> D);
	ObjFunction2 evaluate_point(Point2 &P1);
	void evaluate_point_parallel(Container2 &solutions_requested_for_evaluation);
	void first_search_message();
	void final_search_message(time_t start_time);
	std::string give_the_move(const Point2 &P);
	void numerical_check(Point2 input, std::string const location);
	void show_actual_current_step();
	void HookeJeevesMove_parallel(unsigned int loop);
	Point2  increase(Point2 const &input, unsigned int const index);
	Point2 IntensifyMove2();
	void open_external_monitor_files();
	int PatternMove();
	void Push_Base_Point(const Point2 &P,  std::string  const entry_function);
	template<typename T> void readInputFiles(char const *filepath, std::vector<T> &target_vector );
	void read_design_vector_file(double starting_point_setting, char const *filename, char const *file_description, Point2 &Pnt, const Point2 &lower_bound, const Point2 &upper_bound);
	void read_objectives_file(char const *filename, char const *file_description, ObjFunction2 &ObjPnt);
	Point2 ReduceMove2();
	void save_external_monitor_files();
	void save_optimisation_report(std::string name);
	void save_monitor_data();
	void saveParetoFront();
	void save_external_memories(std::string case_qualifier, int loop);
	void save_base_point(std::string case_qualifier, int loop);
	Point2 sample_and_assess(int discovery_evaluations);
	int stoppingCriteriaNotMet(const unsigned int &n_of_loops,const unsigned int &n_of_evaluations, const unsigned int &n_of_consecutive_improvements);
	void UpdateMemories();
	void update_check_point_file();
	void update_optimisation_search_progress();
	void update_trade_off_progress();


	void show_welcome_messages();
	void check_valid_datum_point();
	void evaluate_datum_point();
	void insert_point_to_MTM_IM_HISTORY(const entry &E);
	void set_initial_step();
	void auto_adapt_initial_step();

	void display_short_report_before_the_optimisation_process();
public:

	const Point2& getDatumPnt() const;


	TabuSearch(  Overall_Optimisation_Configuration_Settings &conf_set,  int restart_flag_argument,
			ObjectiveFunction &ObjectiveFunctionModule,
			Container2 &MediumTermMemory,
			Container2 &IntensificationMemory,
			Container2 &HistoryMemory,
			STM_Container2 &ShortTermMemory,
			LTM_Container2Basic2<double> &LongTermMemory
	);



	void reset();
	double search2();
};


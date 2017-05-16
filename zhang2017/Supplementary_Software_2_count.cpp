/*
 Written by Yixin Zhang
 Modified by Hannes Röst, November 2009
 Modified by Fabian Buller, 14.11.2009
 ESACH VERSION 12.2.2010
 Modified by Michael Stravs, 19. Jan. 2012
 Modified by Moreno Wichert & Adrian Rabenseifner, 23. Apr. 2012
 
 Instructions to use it on a Linux box:
 You will need g++ to compile this file. Use the following comamnd to compile
 g++ -o main.o main_prog.cpp 
 To make the program run faster, use compiler optimization:
 g++ -O3 -omain.o main_prog.cpp
 Use then the following command to run the program:
 ./main.o
 
 if you don't have g++, you need to install it:
 sudo apt-get install g++
 
 Instructions for a Mac:
 you might have to include again 
 #include <Carbon/Carbon.h>
 but I'm not sure, so try it out !!!!
 
 ---------------------------------------------------------------------------
 If you have a segmentation fault, its most likely due to limited stacksize.
 Use the following command to unlimit stacksize: 
 ulimit -s unlimited
 
 for convenience, put it in your ./bashrc file so that it is executed each time
 you start a new bash terminal.r encoding region
 Version 3.2 generates output files automatically, for files larger than your ram you have to use a swap file which is big enough
 Version 3.2.1 generates additionally an output file with columns for the three code library, the cutoffvalue is given in the structure file
 Version 3.2.2 generates an output 0, which display codes which were not assigned to the codelists
 Version 3.2.2JPG include sequence cut for 454 sequences
 Version 3.3.1 76bp reads, dynamic matrix, int cut = 0;
 Version 3.4 Esac 
 Version "3.5" Passing by reference; matrix on the heap (MST 2012-01-19)
 */
#include <iomanip>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include <cstdlib>
#include <sys/time.h>
#include <sys/types.h>
#ifndef _WIN32
#include <sys/sysctl.h>
#endif /* not _WIN32 */

using namespace std;


#define CODELIST_MAX 2000

int CFA=0;					// The number of codes passed the constant region check.
int cut=0;					// To decide the length of sequence for analysis
int highest;				// The highest count.
int lowest = 1000;			// The lowest count.

int normFactor = 100;		// Normalization factor for normalized output

int mismatch_limit = 0;

// The input and output file names
string file_name;    		// The base name for all output files.
string outname_eval;	  	// Matrix output of counts
string outname_evalNorm;	// Matrix output of counts, normalized by column (selection)
string outname_summary;

// Define the variables for the structure file, which describes the structure of the sequences for the decoding analysis.
string structure_file;							// The name of the structure file.
string Input_name;								// The name of raw data.
int min_line_length;					// Minimum length of line (ignore others) - added by Moreno
int noc1, noc2, noc3, noc4;			// The number of encoding regions, and the number of codes at each region.
int noc1b, noc2b, noc3b, noc4b;		// The number of codes at each region: when noc1,2,3,4 = 0, noc1,2,3,4b = 1; else, noc1,2,3,4b = noc1,2,3,4.
string Head1;									// The signs of each lines which containing 4 components.
string Head2;									// The code list file name or the constant sequences.
int starting;							// The start number of a code or a constant region.
int ending;							// The end number of a code or a constant region.

// Output type specification
string output_type_text;
enum output_types {
	OT_2BB  = 2,
	OT_3BB = 3,
	OT_4BB = 4,
};
int output_type;


// The variables of the 11 regions + a template for further correction. (12 X 3)
string Codelist1, Codelist2, Codelist3, Codelist4, Con1, Con2, Con3, Con4, Con5, fraction_length;

int codelist1_start = 0, codelist1_end = 0;
int codelist2_start = 0, codelist2_end = 0;
int codelist3_start = 0, codelist3_end = 0;
int codelist4_start = 0, codelist4_end = 0;
int con1_start = 0, con1_end= 0;
int con2_start = 0, con2_end= 0;
int con3_start = 0, con3_end= 0;
int con4_start = 0, con4_end= 0;
int con5_start = 0, con5_end= 0;

// Count reads
int read_number=0;

// To load the code 1 - 4 lists
int codelist_count1 = 0;
string codelist_array_1[CODELIST_MAX];

int codelist_count2 = 0;
string codelist_array_2[CODELIST_MAX];

int codelist_count3 = 0;
string codelist_array_3[CODELIST_MAX];

int codelist_count4 = 0;
string codelist_array_4[CODELIST_MAX];

// This is the temporary variable for reading every coding region.
int read_cc1, read_cc2, read_cc3, read_cc4;


// M. Stravs 2012-01-19 ex main3.6@2.cpp (markus version)
// Multithreading:
// Two mutexes to lock access to file reading (lock_file) and "global variables" (lock_data),
// and a pointer to the thread array (generated at runtime).
pthread_mutex_t lock_file;
pthread_mutex_t lock_data;
pthread_t * threads;

#define MAX_READER_THREADS 8
// -1 = use number of CPU cores

// Define a struct for the return values of a thread.
// A thread returns the number of sequences read by the thread (thread_reads)
// and the number of sequences used for analysis (CFA).
typedef struct {
	int thread_reads;
	int CFA;
} threaded_analyze_ret;


// Define globally the input file pointer and the dynamic matrix pointer
// (they will be generated at runtime)
ifstream * myfile_p;
int * M;
int M1, M2, M3, M4;

// Define a pointer where the selection sum counter will go
long int * atecl_sel;

// Worker function running a thread for sequence analysis
void* threaded_analyze(void *);
// Analyze a line, and write the match with the codelists into the global dynamic matrix M defined with dimensions M1..M4
int analyze_one_line(const string & line);
// Compare two strings and return the number of character differences
int compare_one_by_one (const string & compare1, const string & compare2);
// Output functions
void output_2bb(ofstream &, bool);
void output_3bb(ofstream &, bool);
void output_4bb(ofstream &, bool);

/******************************************************************************************************************************************/
/************************************************************ The main program ************************************************************/
/************************************************************ starts from here! ***********************************************************/
/******************************************************************************************************************************************/

int main (int argc, char * const argv[]){ 
	
	// To give a general head for all output files.
	
	cout << "Please configure the file structure.txt !"<< endl;
	cout << "Please give a name, which will be the head of all output files" << endl;
	cin >> file_name;	
	
	// Start to read the structure file
	
	structure_file = "structure.txt";
	ifstream infile(structure_file.c_str());
	if (!infile) {cout << "Cannot open file " << structure_file << " for reading." << endl; exit(EXIT_FAILURE);}
	if (infile.is_open()){
		infile >> Input_name;
		infile >> min_line_length;
		infile >> noc1 >> noc2 >> noc3 >> noc4; // The number of coding segments, number of code in codelist 1, 2, 3, and 4.
		infile >> output_type_text;
		infile >> mismatch_limit;
		while (! infile.eof())
		{
			infile >> Head1;
			infile >> starting;
			infile >> ending;
			//	cout << Head1 << " this is the input" << endl;
			infile >> Head2;
			if (Head1 == "x"){codelist1_start = starting; codelist1_end = ending; Codelist1 = Head2;}
			if (Head1 == "y"){codelist2_start = starting; codelist2_end = ending; Codelist2 = Head2;}
			if (Head1 == "z"){codelist3_start = starting; codelist3_end = ending; Codelist3 = Head2;}
			if (Head1 == "$"){codelist4_start = starting; codelist4_end = ending; Codelist4 = Head2;}
			if (Head1 == "1"){con1_start = starting; con1_end= ending; Con1 = Head2;}
			if (Head1 == "2"){con2_start = starting; con2_end= ending; Con2 = Head2;}
			if (Head1 == "3"){con3_start = starting; con3_end= ending; Con3 = Head2;}
			if (Head1 == "4"){con4_start = starting; con4_end= ending; Con4 = Head2;}
			if (Head1 == "5"){con5_start = starting; con5_end= ending; Con5 = Head2;}
		}
	}
	infile.close();
	
	// Set output type correctly
	if((output_type_text == "2bb") | (output_type_text == "2BB"))
		output_type = OT_2BB;
	else if((output_type_text == "3bb") | (output_type_text == "3BB"))
		output_type = OT_3BB;
	else if((output_type_text == "4bb") | (output_type_text == "4BB"))
		output_type = OT_4BB;
	else
		output_type = OT_2BB;


	// Produce a 4-dimensional matrix on the heap - one dimension for every code list
	M1 = noc1 + 1;
	M2 = noc2 + 1;
	M3 = noc3 + 1;
	M4 = noc4 + 1;

	M = new int[ M1*M2*M3*M4 ];

	
	// Initialize with zero
	for (int mn1 = 0; mn1 < M1; ++mn1){
		for (int mn2 = 0; mn2 < M2; ++mn2){
			for (int mn3 = 0; mn3 < M3; ++mn3){
				for (int mn4 = 0; mn4 < M4; ++mn4){
					M[M2 * M3 * M4 * mn1 + M3 * M4 * mn2 + M4 * mn3 + mn4] = 0;
				}}}}
	
	// The four "noc1b" variables are used to compute the average (so if a code is not used, they have to be set to 1).
	if (noc1 == 0) {noc1b=1;} else {noc1b=noc1;}
	if (noc2 == 0) {noc2b=1;} else {noc2b=noc2;}
	if (noc3 == 0) {noc3b=1;} else {noc3b=noc3;}
	if (noc4 == 0) {noc4b=1;} else {noc4b=noc4;}
	
	
	// Document the properties of the structure.txt file to the screen
	
	cout << "The code lists are: " << endl;
	cout << codelist1_start << "-" << codelist1_end << "   " << Codelist1 << endl;
	
	if (codelist2_end != 0) {cout << codelist2_start << "-" << codelist2_end << "   " << Codelist2 << endl;}
	if (codelist3_end != 0) {cout << codelist3_start << "-" << codelist3_end << "   " << Codelist3 << endl;}
	if (codelist4_end != 0) {cout << codelist4_start << "-" << codelist4_end << "   " << Codelist4 << endl;}
	
	cout << "The constant regions are: " << endl;
	cout << con1_start << "-" << con1_end << "   " << Con1 << endl; 
	if (con2_end != 0) {cout << con2_start << "-" << con2_end << "   " << Con2 << endl;}	
	if (con3_end != 0) {cout << con3_start << "-" << con3_end << "   " << Con3 << endl;}
	if (con4_end != 0) {cout << con4_start << "-" << con4_end << "   " << Con4 << endl;} 
	if (con5_end != 0) {cout << con5_start << "-" << con5_end << "   " << Con5 << endl;}
	
	// Set the current time for output file names.
	
	time_t t = time(0);
	tm time = *localtime(&t);
	std::stringstream ss;
	std::string strww;
	ss << "_" << time.tm_year + 1900 << "_" << time.tm_mon + 1 << "_" << time.tm_mday << "_" ;
	ss << time.tm_hour << "_" << time.tm_min << "_" << time.tm_sec;
	ss >> strww;
	cout << strww;	
	
	// Start to input the code lists (1 - 4).
	// moved this part up because we need to do this first. hr.
	string line_in;
	
	ifstream codelist_in1(Codelist1.c_str());
	if (codelist_in1.is_open())
	{while (! codelist_in1.eof())
	{getline (codelist_in1, line_in);
		codelist_count1 = codelist_count1 + 1;
		codelist_array_1[codelist_count1] = line_in;
		codelist_array_1[codelist_count1].resize(codelist1_end - codelist1_start +1);
		cout << codelist_count1 << "-----" << codelist_array_1[codelist_count1]  << "\n";}}
	codelist_in1.close();
	
	if (codelist2_end != 0){
		ifstream codelist_in2(Codelist2.c_str());
		if (codelist_in2.is_open())
		{while (! codelist_in2.eof())
		{getline (codelist_in2, line_in);
			codelist_count2 = codelist_count2 + 1;
			codelist_array_2[codelist_count2] = line_in;
			codelist_array_2[codelist_count2].resize(codelist2_end - codelist2_start +1);
			cout << codelist_count2 << "-----" << codelist_array_2[codelist_count2]  << "\n";}}
		codelist_in2.close();
	}

	
	if (codelist3_end != 0){
		ifstream codelist_in3(Codelist3.c_str());
		if (codelist_in3.is_open())
		{while (! codelist_in3.eof())
		{getline (codelist_in3, line_in);
			codelist_count3 = codelist_count3 + 1;
			codelist_array_3[codelist_count3] = line_in;
			codelist_array_3[codelist_count3].resize(codelist3_end - codelist3_start +1);
			cout << codelist_count3 << "-----" << codelist_array_3[codelist_count3]  << "\n";}}
		codelist_in3.close();
	}

	
	if (codelist4_end != 0){
		ifstream codelist_in4(Codelist4.c_str());
		if (codelist_in4.is_open())
		{while (! codelist_in4.eof())
		{getline (codelist_in4, line_in);
			codelist_count4 = codelist_count4 + 1;
			codelist_array_4[codelist_count4] = line_in;
			codelist_array_4[codelist_count4].resize(codelist4_end - codelist4_start +1);
			cout << codelist_count4 << "-----" << codelist_array_4[codelist_count4]  << "\n";}}
		codelist_in4.close();
	}

	cout << "... Pre-loading of code finished successfully." << "\n";

	// To read to raw data file
	
	//if(max_reader_threads < 0) {
		// get number of CPU cores (Mac OS X)
		//int buf;
		//size_t size = sizeof buf;
		//sysctlbyname("machdep.cpu.core_count", &buf, &size, NULL, 0);
		//max_reader_threads = buf;
	//}
	
	
	cout << "Using " << MAX_READER_THREADS << " threads..." << endl;
	
	// Measure processing time: get starting time
	timeval t0;
	gettimeofday(&t0, NULL);
	
	// open the input file
    myfile_p = new ifstream(Input_name.c_str());

    if (myfile_p->is_open())
	{

    	// Create the mutexes to ensure thread safety
		pthread_mutex_init(&lock_file, NULL);
		pthread_mutex_init(&lock_data, NULL);
		
		// Generate threads
		threads = new pthread_t[MAX_READER_THREADS];
		
		// Start concurrent analysis
		for(int n=0; n<MAX_READER_THREADS;n++)
		{
			pthread_create(&threads[n], NULL, threaded_analyze, NULL);
			cout << "started thread " << n << endl;
		}
        
		void * exitstatus;
		threaded_analyze_ret * p_tret;

		// Wait for all threads to finish and count their contributions
		for(int n=0; n<MAX_READER_THREADS;n++)
		{
			pthread_join(threads[n], &exitstatus);
			p_tret = (threaded_analyze_ret *) exitstatus;
			cout << "terminated thread " << n << " ( " << p_tret->thread_reads << " reads, " << p_tret->CFA << " counted )" << endl;
			CFA += p_tret->CFA;
			delete p_tret;
		}
		
		// Destroy the mutexes
		pthread_mutex_destroy(&lock_file);
		pthread_mutex_destroy(&lock_data);
		// all threads are stopped and deleted
		delete[] threads;
        

	}
    myfile_p->close();
    delete(myfile_p);
	
	// Get final time
	timeval t1;
	gettimeofday(&t1, NULL);
	
	cout << "read_number = " << read_number << "\n";	
	
	
	
	// Calculate the highest counts and to assign the counts for statistic evaluations
	for (int mm1 = 0; mm1 < M1; ++mm1){
		for (int mm2 = 0; mm2 < M2; ++mm2){
			for (int mm3 = 0; mm3 < M3; ++mm3){
				for (int mm4 = 0; mm4 < M4; ++mm4){
					// If highest/lowest exceeded, adjust the values
					if (M[M2 * M3 * M4 * mm1 + M3 * M4 * mm2 + M4 * mm3 + mm4] > highest)
					{highest = M[M2 * M3 * M4 * mm1 + M3 * M4 * mm2 + M4 * mm3 + mm4];}
					if (M[M2 * M3 * M4 * mm1 + M3 * M4 * mm2 + M4 * mm3 + mm4] < lowest)
					{lowest = M[M2 * M3 * M4 * mm1 + M3 * M4 * mm2 + M4 * mm3 + mm4];}
				}}}}
	
	cout << endl << "The highest count is: " << highest << endl;
	cout << endl << "The lowest count is: " << lowest << endl;
	float average;
	average = ceil(float(CFA)/float(noc1b)/float(noc2b)/float(noc3b)/float(noc4b));
	cout << endl << "The average count is: " << average << endl;
	
	// Start to report the results
	
	cout << "The decoding has been finished." << endl;
	
	// Set the output filenames
	outname_eval = file_name; outname_eval.append(strww); outname_eval.append("_eval.txt");
	outname_evalNorm = file_name; outname_evalNorm.append(strww); outname_evalNorm.append("_evalNorm.txt");
	outname_summary = file_name; outname_summary.append(strww); outname_summary.append("_Codecounts.txt");


	// Summary report containing:
	// the input data structure
	// an assignment count for each codelist
	ofstream report_summary;
	report_summary.open (outname_summary.c_str());

	// Print input data information to report
	report_summary <<Input_name<< endl;
	report_summary << "The total number of sequences received is: " << read_number << endl;
	report_summary << "Total number of sequences used for Code analysis:	" << CFA << endl;

	// Print data structure to report
	report_summary << "The structure of your inputed data is following: " << endl;
	report_summary << endl;

	report_summary << "The code lists are: " << endl;
	report_summary << codelist1_start << "-" << codelist1_end << "   " << Codelist1 << endl;
	if (codelist2_end != 0) {report_summary << codelist2_start << "-" << codelist2_end << "   " << Codelist2 << endl;}
	if (codelist3_end != 0) {report_summary << codelist3_start << "-" << codelist3_end << "   " << Codelist3 << endl;}
	if (codelist4_end != 0) {report_summary << codelist4_start << "-" << codelist4_end << "   " << Codelist4 << endl;}
	report_summary << "The constant regions are: " << endl;
	report_summary << con1_start << "-" << con1_end << "   " << Con1 << endl;
	report_summary << con2_start << "-" << con2_end << "   " << Con2 << endl;
	if (con3_end != 0) {report_summary << con3_start << "-" << con3_end << "   " << Con3 << endl;}
	if (con4_end != 0) {report_summary << con4_start << "-" << con4_end << "   " << Con4 << endl;}
	if (con5_end != 0) {report_summary << con5_start << "-" << con5_end << "   " << Con5 << endl;}
	report_summary<<endl;

	// Build the codelist assignment count arrays and initialize with zero
	int atecl_1[M1], atecl_2[M2], atecl_3[M3], atecl_4[M4];
	// Build the selection assignment count array and initialize with zero
	int nsel = 0;
	if(output_type == OT_2BB)
		nsel = M1*M2;
	else if(output_type == OT_3BB)
		nsel = M1;
	else
		nsel = 1;

	atecl_sel = new long int[nsel];

	for (int atecl_n = 0; atecl_n < M1; ++atecl_n){atecl_1[atecl_n]=0;}
	for (int atecl_n = 0; atecl_n < M2; ++atecl_n){atecl_2[atecl_n]=0;}
	for (int atecl_n = 0; atecl_n < M3; ++atecl_n){atecl_3[atecl_n]=0;}
	for (int atecl_n = 0; atecl_n < M4; ++atecl_n){atecl_4[atecl_n]=0;}
	for (int atecl_n = 0; atecl_n < nsel; ++atecl_n){atecl_sel[atecl_n]=0;}

	// Count data points assigned to respective codes
	for (int mm1 = 0; mm1 < M1; ++mm1){
		for (int mm2 = 0; mm2 < M2; ++mm2){
			for (int mm3 = 0; mm3 < M3; ++mm3){
				for (int mm4 = 0; mm4 < M4; ++mm4){
					// This counts the code hits for every codelist entry
					atecl_1[mm1]=atecl_1[mm1]+M[M2 * M3 * M4 * mm1 + M3 * M4 * mm2 + M4 * mm3 + mm4];
					atecl_2[mm2]=atecl_2[mm2]+M[M2 * M3 * M4 * mm1 + M3 * M4 * mm2 + M4 * mm3 + mm4];
					atecl_3[mm3]=atecl_3[mm3]+M[M2 * M3 * M4 * mm1 + M3 * M4 * mm2 + M4 * mm3 + mm4];
					atecl_4[mm4]=atecl_4[mm4]+M[M2 * M3 * M4 * mm1 + M3 * M4 * mm2 + M4 * mm3 + mm4];
					// This counts code hits per column; later used in the normalized matrix output
					// (but not in the report)
					// Only columns with all codes assigned are counted (so a column adds up to 1)
					if(mm1!=0 && mm2!=0 && mm3!=0 && mm4!=0)
					{
						// Depending on the selected output style, selections are counted differently!
						int selectionHit = 0;
						if(output_type == OT_2BB)
							selectionHit = M2*mm1 + mm2;
						else if(output_type == OT_3BB)
							selectionHit = mm1;
						else
							selectionHit = 0;
						atecl_sel[selectionHit] += M[M2 * M3 * M4 * mm1 + M3 * M4 * mm2 + M4 * mm3 + mm4];
					}
				}
			}
		}
	}


	// Print assignment counts to report
	report_summary << endl;
	report_summary << "Assigned to codelist 1" << endl;
	for (int mm1 = 0; mm1 < M1; ++mm1){report_summary << mm1 << "\t" << atecl_1[mm1] << endl;}
	
	if (codelist2_end != 0) {report_summary << "Assigned to codelist 2" << endl;
		for (int mm2 = 0; mm2 < M2; ++mm2){report_summary << mm2 << "\t" << atecl_2[mm2] << endl;}}
	
	if (codelist3_end != 0) {report_summary << "Assigned to codelist 3" << endl;
		for (int mm3 = 0; mm3 < M3; ++mm3){report_summary << mm3 << "\t" << atecl_3[mm3] << endl;}}
	
	if (codelist4_end != 0) {report_summary << "Assigned to codelist 4" << endl;
		for (int mm4 = 0; mm4 < M4; ++mm4){report_summary << mm4 << "\t" << atecl_4[mm4] << endl;}}
	
	report_summary<<"The average frequency is: "<< average <<endl;
	report_summary<<"The highest frequency is: "<< highest << endl;
	report_summary<<"The lowest frequency is: "<<lowest << endl;
	report_summary.close();


	// Matrix report with counts for respective code combination:
	// Code1_Code2 combinations are in columns,
	// Code3_Code4 (called CodeA, CodeB in output) combination are in rows

	ofstream report_eval;
	report_eval.open (outname_eval.c_str());
	if(output_type == OT_2BB)
		output_2bb(report_eval, false);
	else if(output_type == OT_3BB)
		output_3bb(report_eval, false);
	else if(output_type == OT_4BB)
		output_4bb(report_eval, false);
	report_eval.close();


	// Matrix report with counts for respective code combination, normalized by column (ie by selection):
	// Code1_Code2 combinations (selections) are in columns,
	// Code3_Code4 (called CodeA, CodeB in output) combination are in rows

	// The normalization is calculated with the previously counted total counts per selection (atecl_sel):
	// normalized count cnorm = codecombis * normFactor * count / totalCounts
	// where normFactor is the baseline factor to normalize to (100),
	// codecombis is the number of possible code combinations of codes 3 and 4,
	// count is the number of hits for the specific code
	// and totalCounts are the total hit counts for the selection (ie code1-code2 combination).

	report_eval.open (outname_evalNorm.c_str());
	if(output_type == OT_2BB)
		output_2bb(report_eval, true);
	else if(output_type == OT_3BB)
		output_3bb(report_eval, true);
	else if(output_type == OT_4BB)
		output_4bb(report_eval, true);
	report_eval.close();

	delete[] M;
	delete[] atecl_sel;
	
	cout << "computation took: " << (t1.tv_sec - t0.tv_sec) << " seconds." << endl;
	
} // end of main()


/***************************************************************************
 * This is the end of the main progam 
 * here are the functions
 ***************************************************************************/

// Output for 2BB library:
// This output interprets codes 1 and 2 as selection codes and codes 3,4 as building block codes A,B.
// It prints a normalized (per selection) or non-normalized table with selections in columns and building block combinations in rows.
// The header contains the combined selection ID in format Code1_Code2.
// The last 4 columns contain the building block IDs and a composed code of format CodeA__CodeB.
void output_2bb(ofstream & out, bool normalize = false)
{
	for (int mm2 = 1; mm2 < M2; ++mm2){
		for (int mm1 = 1; mm1 < M1; ++mm1){
			out << mm1 <<"_"<< mm2 <<",";}}
	out <<"Combined,CodeA,CodeB"<< endl;

	int codecombis = noc3b * noc4b;

	if (codelist3_end != 0) {
		for (int mm4 = 1; mm4 < M4; ++mm4){
			for (int mm3 = 1; mm3 < M3; ++mm3){
				for (int mm2 = 1; mm2 < M2; ++mm2){
					for (int mm1 = 1; mm1 < M1; ++mm1){
						if(!normalize)
							out << M[M2 * M3 * M4 * mm1 + M3 * M4 * mm2 + M4 * mm3 + mm4] << ",";
						else
						{
							double cnorm;
							if(atecl_sel[M2*mm1 + mm2] > 0)
								cnorm = round((double)(codecombis * normFactor) * (double)(M[M2 * M3 * M4 * mm1 + M3 * M4 * mm2 + M4 * mm3 + mm4]) / (double)(atecl_sel[M2*mm1 + mm2]));
							else
								cnorm = 0;
							out << cnorm << ",";
						}

					}
				}
				out << mm3 <<"__"<< mm4 <<"," << mm3 <<"," << mm4 <<  endl;
			}
		}
	}
}

// Output for 3BB library:
// This output interprets code 1 as selection code and code 2,3,4 as building block codes A,B,C.
// It prints a normalized (per selection) or non-normalized table with selections in columns and building block combinations in rows.
// The header contains the selection ID.
// The last 4 columns contain the building block IDs and a composed code of format CodeA__CodeB__CodeC.
void output_3bb(ofstream & out, bool normalize = false)
{
	for (int mm1 = 1; mm1 < M1; ++mm1){
			out << "S" << mm1 <<",";}
	out <<"Combined,CodeA,CodeB,CodeC"<< endl;

	long int codecombis = noc2b * noc3b * noc4b;

	if (codelist3_end != 0) {
		for (int mm4 = 1; mm4 < M4; ++mm4){
			for (int mm3 = 1; mm3 < M3; ++mm3){
				for (int mm2 = 1; mm2 < M2; ++mm2){
					for (int mm1 = 1; mm1 < M1; ++mm1){
						if(!normalize)
							out << M[M2 * M3 * M4 * mm1 + M3 * M4 * mm2 + M4 * mm3 + mm4] << ",";
						else
						{
							double cnorm;
							if(atecl_sel[mm1] > 0)
								cnorm = round((double)(codecombis * normFactor) * (double)(M[M2 * M3 * M4 * mm1 + M3 * M4 * mm2 + M4 * mm3 + mm4]) / (double)(atecl_sel[mm1]));
							else
								cnorm = 0;
							out << cnorm << ",";
						}

					}
					out << mm2 << "__" << mm3 <<"__"<< mm4 <<"," << mm2 << "," << mm3 <<"," << mm4 <<  endl;
				}
			}
		}
	}
}


// Output for 4BB library:
// This output interprets all codes (1,2,3,4) as building block codes A,B,C,D in one selection.
// It prints a normalized (over total counts) or non-normalized table with building block combinations in rows.
// The last 5 columns contain the building block IDs and a composed code of format CodeA__CodeB__CodeC__CodeD.
void output_4bb(ofstream & out, bool normalize = false)
{
	out <<"Counts,Combined,CodeA,CodeB,CodeC,CodeD"<< endl;

	long int codecombis = noc1b * noc2b * noc3b * noc4b;

	if (codelist3_end != 0) {
		for (int mm4 = 1; mm4 < M4; ++mm4){
			for (int mm3 = 1; mm3 < M3; ++mm3){
				for (int mm2 = 1; mm2 < M2; ++mm2){
					for (int mm1 = 1; mm1 < M1; ++mm1){
						if(!normalize)
							out << M[M2 * M3 * M4 * mm1 + M3 * M4 * mm2 + M4 * mm3 + mm4] << ",";
						else
						{
							double cnorm;
							if(atecl_sel[0] > 0)
								cnorm = round((double)(codecombis * normFactor) * (double)(M[M2 * M3 * M4 * mm1 + M3 * M4 * mm2 + M4 * mm3 + mm4]) / (double)(atecl_sel[0]));
							else
								cnorm = 0;
							out << cnorm << ",";
						}
						out << mm1 << "__" << mm2 << "__" << mm3 <<"__"<< mm4 <<"," << mm1 << ","
								<< mm2 << "," << mm3 <<"," << mm4 <<  endl;
					}
				}
			}
		}
	}
}


void * threaded_analyze(void * arg)
{
	string line_pre;	//A line read from raw data.
	string line_in;		//A ligated line from the raw data, before it is terminated by <
	int line_size;		//The variable for line size, it is used to elimate the last blank.
	int line_check;
	int thread_reads = 0;
	int CFA = 0;

	pthread_mutex_lock(&lock_file);
	while (! (myfile_p->eof()))
    {
        line_check = 1;
        line_in.clear();
        
        while (line_check != 0)
        {
            
            getline (*myfile_p, line_pre);
			// change minimal line size from 0 to min_line_length
            if (line_pre[0] != '>' && line_pre.size() > min_line_length) 
            {
                line_in.append(line_pre);
                line_size = line_in.size();
                if (int(line_in[line_size-1]) == 13) { 
                    //0x0A (in decimal 10) is a newline '\n'
                    //0x0D (in decimal 13) is a carriage return '\r'
                    //Windows uses '\r\n'
                    //Mac     uses '\r'
                    //Unix    uses '\n'
                    //thus if its a windows file, the last character could be a 13
                    //and we need to remove it
                    line_in.resize(line_size-1);
                }
                line_check = line_check + 1;}
            if (line_pre[0] == '>')  {line_check = 0;} 
            if (line_pre.size() == 0)  {line_check = 0;}
        }
        pthread_mutex_unlock(&lock_file);
        //analyze the line WITHOUT storing it 
        
        // analyze_one_line() returns 1 if line was processed (ie constant regions were matching), 0 otherwise
        // Add this to the count
        // Note: the analysis results are written directly into the matrix M
        CFA += analyze_one_line(line_in);
        // Increment the thread's read counter (for statistics)
        thread_reads++;

        // This is a potential slowdown but helps us to keep track of processing
        pthread_mutex_lock(&lock_data);
        read_number++;
        pthread_mutex_unlock(&lock_data);
        
        // User notification for progress report
        // Note: we don't threadlock here so the printed read number can already be x000001 (since processing in parallel threads can continue)
        if (read_number % 100000 == 0) cout << "analyzing sequence nr. " << read_number << endl;
        pthread_mutex_lock(&lock_file);
    }
    pthread_mutex_unlock(&lock_file);

    threaded_analyze_ret * p_tret = new threaded_analyze_ret;
    p_tret->thread_reads = thread_reads;
    p_tret->CFA = CFA;
    return p_tret;
}


int analyze_one_line(const string & line_exam)
{
	string str_read;

	// This is the temporary variable for reading every coding region.
	int read_cc1, read_cc2, read_cc3, read_cc4;	
	int CFA = 0;
	int go_on = 0; 
	
	// For each sequence, check the constant regions.
	if (line_exam.size() > cut){

		// Check the constant regions for
		
		if ((con1_end!=0)&(go_on==0)){
			str_read = line_exam.substr((con1_start - 1), (con1_end - con1_start + 1));
			if (compare_one_by_one(str_read, Con1) > mismatch_limit) {go_on = 1;}
		}
		
		if ((con2_end!=0)&(go_on==0)){
			str_read = line_exam.substr((con2_start - 1), (con2_end - con2_start + 1));
			if (compare_one_by_one(str_read, Con2) > mismatch_limit) {go_on = 1;}
		}	
		
		if ((con3_end!=0)&(go_on==0)){
			str_read = line_exam.substr((con3_start - 1), (con3_end - con3_start + 1));
			if (compare_one_by_one(str_read, Con3) > mismatch_limit) {go_on = 1;}
		}	
		
		if ((con4_end!=0)&(go_on==0)){
			str_read = line_exam.substr((con4_start - 1), (con4_end - con4_start + 1));
			if (compare_one_by_one(str_read, Con4) > mismatch_limit) {go_on = 1;}
		}	
		
		if ((con5_end!=0)&(go_on==0)){
			str_read = line_exam.substr((con5_start - 1), (con5_end - con5_start + 1));
			if (compare_one_by_one(str_read, Con5) > mismatch_limit) {go_on = 1;}
		}
		
		// If all the constant regions matched, try to assign the coding regions to their respective code
		if (go_on == 0)
		{
			CFA=1;
			read_cc1=0;
			read_cc2=0;
			read_cc3=0;
			read_cc4=0;
			str_read = line_exam.substr((codelist1_start - 1), (codelist1_end - codelist1_start + 1));

			int min = 0;
			int diff;
			for (int cc1 = 1; cc1 < M1; ++cc1){
				diff = compare_one_by_one(str_read, codelist_array_1[cc1]);
				if (diff  <= min) 
				{read_cc1=cc1; min = diff;}}
			
			if (codelist2_end != 0){
				str_read = line_exam.substr((codelist2_start - 1), (codelist2_end - codelist2_start + 1));
				min = 0;
				for (int cc2 = 1; cc2 < codelist_count2+1; ++cc2){
					diff = compare_one_by_one(str_read, codelist_array_2[cc2]);
					if (diff  <= min) 
					{read_cc2=cc2; min = diff;}}}
			// If there is no codelist2 (also 3, 4 as below), count as a hit to code 1.
			// This is so we can use the same report format.
			else
				read_cc2 = 1;
			
			if (codelist3_end != 0){
				str_read = line_exam.substr((codelist3_start - 1), (codelist3_end - codelist3_start + 1));

				min = 0;
				for (int cc3 = 1; cc3 < codelist_count3+1; ++cc3){
					diff = compare_one_by_one(str_read, codelist_array_3[cc3]);
					if (diff <= min)
					{read_cc3 = cc3; min = diff;}}}
			else
				read_cc3 = 1;
			
			if (codelist4_end != 0){
				str_read = line_exam.substr((codelist4_start - 1), (codelist4_end - codelist4_start + 1));

				min = 0;
				for (int cc4 = 1; cc4 < codelist_count4+1; ++cc4){
					diff = compare_one_by_one(str_read, codelist_array_4[cc4]);
					if (diff <= min) 
					{read_cc4 = cc4; min = diff;}}}
			else
				read_cc4 = 1;
			
            // lock local variables for writing:
            pthread_mutex_lock(&lock_data);
            // Calculate position of hit in the dynamic matrix and increment
            // The hit is for codes read_cc1 .. read_cc4
            M[ M2*M3*M4*read_cc1 + M3*M4*read_cc2 + M4*read_cc3 + read_cc4] = 
            M[ M2*M3*M4*read_cc1 + M3*M4*read_cc2 + M4*read_cc3 + read_cc4] +1;
            pthread_mutex_unlock(&lock_data);
			
		} // if (go_on == 0)
	} // if (line_exam.size() > cut)
	
	// If line was counted for analysis then 1 is returned which can be added to the total counter there
	return(CFA);

} // analyze_one_line



int compare_one_by_one (const string & compare1, const string & compare2){
	int compared_difference = 0;
	int compare_end;
	if (compare1.size() >= compare2.size()){compare_end = compare2.size();}
	if (compare1.size() < compare2.size()){compare_end = compare1.size();}
	for (int compare_counting = 0; compare_counting < compare_end; ++compare_counting){
		if (compare1[compare_counting] != compare2[compare_counting]){
			compared_difference = compared_difference + 1;}
	}
	if (compare1.size() != compare2.size()) {compared_difference = 99;}
	return compared_difference;
}

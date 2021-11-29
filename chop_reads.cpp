/*
 *  chop_reads.cpp
 *  
 *  Chops the desired number of bases from the .fastq / .fa file  
 * 
 *  Created by Aziz Mithani on 20/04/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include <iostream>
#include <fstream>
#include <string>

using namespace std;

//const int NEW_LENGTH = 51;

int main (int argc, char** argv) {

//	string ReadFilename = "/san/001/BGI/Data from BGI/Data from 13th April 2010/Salt-1_2.fq";
//	string outFilename =  "/san/001/BGI/Data from BGI/Data from 13th April 2010/Salt-1_2_x.fq";
	string Read_File(argv[1]);
	int readLength(atoi(argv[2]));
	
	// open the read file stream
	ifstream ifs ( Read_File.c_str() );
	
	// open a file for the output
//	ofstream ofs (OUT_File.c_str());
	
	string str = "";
	while (ifs.good()) {
		// get the first line and write it to the output unmodified - readname
		getline(ifs, str);
		cout << str << endl;
		
		// get the read
		getline(ifs, str);
		// extract the required no. of bases and write it to the output
		str = str.substr(0, readLength);
		cout << str << endl;
		
		// get the third line and write it to the output unmodified
		getline(ifs, str);
		cout << str << endl;
		
		// get the quality score
		getline(ifs, str);
		// chop the required no. of bases and write it to the output
		str = str.substr(0, readLength);
		cout << str << endl;
	}
	
//	ofs.close();
	
	return 0;
}

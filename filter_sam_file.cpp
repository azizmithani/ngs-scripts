/*
 *  filter_sam_file.cpp
 *  
 *	Removes low quality, unpaired and non-unique reads from a sam file
 *
 *  Created by Aziz Mithani on 15/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include "Utilities.h"
#include <stdlib.h>

using namespace std;


bool IGNORE_UNPAIRED_READS = true;
const int MQUAL_THRESHOLD = 20;
const bool IGNORE_NON_UNIQUE_READS = true;

const int FLAG_SEQUENCE_UNMAPPED = 4;
const int FLAG_MATE_UNMAPPED = 8;

const double MINIMUM_CONSENSUS = 0.9;
const int MINIMUM_COVERAGE = 20;


void FilterSAMFile(string SAM_File, int Max_Insert_Size) {
/**
	Removes unwanted reads
 */
	
	// open the sam file stream
	ifstream ifs_sam ( SAM_File.c_str() );

	string strRead1="";
	string strRead2="";
	while (ifs_sam.good()) {
		getline(ifs_sam, strRead1);
		
		if (strRead1.empty()) { // ignore emtry lines
			continue;
		} else if (strRead1[0] == '@') { // ignore the header
			cout << strRead1 << endl;
			continue;
		}

		// get the 2nd read in the pair
		getline(ifs_sam, strRead2);
		
		if (strRead2.empty()) { // ignore emtry lines
			continue;
		} 
		
		vector<string> entries1;
		strsplit(strRead1, entries1, "\t");

		vector<string> entries2;
		strsplit(strRead2, entries2, "\t");

		// read 1
		bool PrintRead1 = true;
		if (atoi(entries1[4].c_str()) < MQUAL_THRESHOLD) {
			PrintRead1 = false;
		} else if (atoi(entries1[1].c_str()) & FLAG_SEQUENCE_UNMAPPED) {
			PrintRead1 = false;
		} else if (IGNORE_NON_UNIQUE_READS && entries1[11].find("XT:A:U") == string::npos) {
			PrintRead1 = false;
		} else if (abs(atoi(entries1[8].c_str())) > Max_Insert_Size) {
			PrintRead1 = false;
		}

		// read 2
		bool PrintRead2 = true;
		if (atoi(entries2[4].c_str()) < MQUAL_THRESHOLD) {
			PrintRead2 = false;
		} else if (atoi(entries2[1].c_str()) & FLAG_SEQUENCE_UNMAPPED) {
			PrintRead2 = false;
		} else if (IGNORE_NON_UNIQUE_READS && entries2[11].find("XT:A:U") == string::npos) {
			PrintRead2 = false;
		} else if (abs(atoi(entries2[8].c_str())) > Max_Insert_Size) {
			PrintRead2 = false;
		}
		
		if (IGNORE_UNPAIRED_READS) {
			if (PrintRead1 && PrintRead2) {
				cout << strRead1 << endl;
				cout << strRead2 << endl;
			}
		} else {
			if (PrintRead1) {
				cout << strRead1 << endl;
			}
			if (PrintRead2) {
				cout << strRead2 << endl;
			}
		}

	} // end while
	
	// close the file
	ifs_sam.close();

}

int main (int argc, char** argv) {
	
	if (argc < 2) {
		cout << "SAM File not specified" << endl;
		exit(0);
	}
	string SAM_File (argv[1]);
	int Max_Insert_Size = atoi(argv[2]);
	if ( argc >= 4 && atoi(argv[3]) == 0 ) {
		IGNORE_UNPAIRED_READS = false;
	}
	 

	FilterSAMFile(SAM_File, Max_Insert_Size);
			
	return 0;
}


/* Changed on 1st August 2010. The following code keeps a read, if it's mapped in pair but the mate is either low quality or non-uniquely mapped
 void FilterSAMFile(string SAM_File) {
 
 // open the sam file stream
 ifstream ifs_sam ( SAM_File.c_str() );
 
 string strRead="";
 while (ifs_sam.good()) {
 getline(ifs_sam, strRead);
 
 if (strRead.empty()) { // ignore emtry lines
 continue;
 } else if (strRead[0] == '@') { // ignore the header
 cout << strRead << endl;
 continue;
 }
 
 vector<string> entries;
 strsplit(strRead, entries, "\t");
 
 if (atoi(entries[4].c_str()) <= MQUAL_THRESHOLD) {
 continue;
 } else {
 
 if (IGNORE_UNPAIRED_READS && (atoi(entries[1].c_str()) & FLAG_SEQUENCE_UNMAPPED || atoi(entries[1].c_str()) & FLAG_MATE_UNMAPPED) ) {
 continue;
 } else if (IGNORE_NON_UNIQUE_READS && entries[11].find("XT:A:U") == string::npos) {
 continue;
 }
 }
 
 
 
 // This read has passed the filter ... output the read
 cout << strRead << endl;
 
 } // end while
 
 // close the file
 ifs_sam.close();
 
 }
 */

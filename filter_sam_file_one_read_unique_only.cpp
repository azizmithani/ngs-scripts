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


const int MQUAL_THRESHOLD = 0;
const bool IGNORE_UNPAIRED_READS = true;

const int FLAG_SEQUENCE_UNMAPPED = 4;
const int FLAG_MATE_UNMAPPED = 8;

const double MINIMUM_CONSENSUS = 0.9;
const int MINIMUM_COVERAGE = 20;

/*
 void strsplit(const string& str, vector<string>& tokens, const string& delimiters = " ") {
 // Skip delimiters at beginning.
 string::size_type lastPos = str.find_first_not_of(delimiters, 0);
 // Find first "non-delimiter".
 string::size_type pos = str.find_first_of(delimiters, lastPos);
 
 while (string::npos != pos || string::npos != lastPos) {
 // Found a token, add it to the vector.
 tokens.push_back(str.substr(lastPos, pos - lastPos));
 // Skip delimiters.  Note the "not_of"
 lastPos = str.find_first_not_of(delimiters, pos);
 // Find next "non-delimiter"
 pos = str.find_first_of(delimiters, lastPos);
 }
 }
 */

void FilterSAMFile(string SAM_File) {
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
		
		if (atoi(entries1[4].c_str()) < MQUAL_THRESHOLD || atoi(entries2[4].c_str()) < MQUAL_THRESHOLD) { // ignore if either of the reads is low quality
			continue;
		} else {
			
			if (IGNORE_UNPAIRED_READS && (atoi(entries1[1].c_str()) & FLAG_SEQUENCE_UNMAPPED || atoi(entries1[1].c_str()) & FLAG_MATE_UNMAPPED) ) {
				continue;
			} else if ((entries1[11].find("XT:A:U") == string::npos && entries2[11].find("XT:A:U") != string::npos) || (entries1[11].find("XT:A:U") != string::npos && entries2[11].find("XT:A:U") == string::npos)) {
				// do nothing
			} else {
				continue;
			}
			
		}
		
		
		
		// This read has passed the filter ... output the read
		cout << strRead1 << endl;
		cout << strRead2 << endl;
		
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
	
	FilterSAMFile(SAM_File);
	
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

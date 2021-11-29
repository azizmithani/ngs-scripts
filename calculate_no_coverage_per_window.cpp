/*
 *
 *  Created by Aziz Mithani on 07/05/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include "Utilities.h"

using namespace std;


void CalculateNoCoverage(string Header_File, string No_Coverage_File, string OUT_File, int WindowSize) {
	
	map < string, int* > PositionCount;
	map < string, int > RefLength;

	ifstream ifs_header( Header_File.c_str() );
	string strHeader = "";
	while (ifs_header.good()) {
		getline(ifs_header, strHeader);
		
		if (strHeader.empty()) { // ignore emtry lines
			continue;
		} else if (strHeader.find("@SQ") != string::npos) { // get sequence length from the header
			vector<string> entries;
			strsplit(strHeader, entries, "\t");
			
			string reference = entries[1].substr(3, entries[1].length());
			string strLength = entries[2].substr(3, entries[2].length());
			int length = atoi(strLength.c_str());
			
			int nWindows = length / WindowSize + 1;
			// memory initialisation
			int *RefPositionCount = new int[nWindows];
			for(int i = 0; i < nWindows; i++) {
				*(RefPositionCount+i) = 0;
			}
			
			RefLength[reference] = length;
			PositionCount[reference] = RefPositionCount;
		}
	} // end while
	
	ifs_header.close();
	
	ifstream ifs_sam( No_Coverage_File.c_str() );
	string strNo_Coverage = "";
	string reference = "";
	string previous_reference = "";
	int* RefPositionCount;
	while (ifs_sam.good()) {
		getline(ifs_sam, strNo_Coverage);
		
		if (strNo_Coverage.empty()) { // ignore emtry lines
			continue;
		} else {
			vector<string> entries;
			strsplit(strNo_Coverage, entries, "\t");
			
			string reference = entries[0];
			int pos = atoi(entries[1].c_str());
			
			if (reference.compare(previous_reference) != 0) {
				if (previous_reference.length() > 0) {
					//break;
					PositionCount[previous_reference] = RefPositionCount;
				}
				
				RefPositionCount = PositionCount[reference];
				previous_reference = reference;
			}
			
			int window = pos / WindowSize;
			*(RefPositionCount + window) = *(RefPositionCount + window) + 1;
			
			//cout << reference << "\t" << pos << "\t" << window << "\t" << atoi(entries[3].c_str()) << "\t" << *(RefCoverage + window) << endl;
		}
	} // end while
	PositionCount[previous_reference] = RefPositionCount;
	
	ifs_sam.close();
	
	ofstream ofs(OUT_File.c_str());
	map < string, int* >::iterator itReadCount = PositionCount.begin();
	for (; itReadCount != PositionCount.end(); itReadCount++) {
		reference = itReadCount->first;
		RefPositionCount = itReadCount->second;
		
		int length = RefLength[reference];
		int nWindows = length / WindowSize + 1;
		for (int i = 0; i < nWindows; i++) {
			ofs << reference << "\t" << WindowSize * i + 1 << "\t" << WindowSize * (i + 1) << "\t" << *(RefPositionCount + i) << endl;
		}

	}
	ofs.close();
}


int main (int argc, char** argv) {
	
	if (argc < 3) {
		cout << "Files not specified" << endl;
		exit(0);
	}
	
	string Header_File (argv[1]);
	string No_Coverage_File (argv[2]);
	string OUT_File (argv[3]);
	int WindowSize (atoi(argv[4]));

	CalculateNoCoverage(Header_File, No_Coverage_File, OUT_File, WindowSize);

	return 0;
}

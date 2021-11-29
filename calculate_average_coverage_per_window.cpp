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
#include <stdlib.h>

using namespace std;

int WINDOW_SIZE = 10000;

void CalculateCoverage(string Header_File, string Pileup_File, string OUT_File) {
	
	map < string, int* > Coverage;
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
			
			int nWindows = length / WINDOW_SIZE + 1;
			// memory initialisation
			int *RefCoverage = new int[nWindows];
			for(int i = 0; i < nWindows; i++) {
				*(RefCoverage+i) = 0;
			}
			
			RefLength[reference] = length;
			Coverage[reference] = RefCoverage;
		}
	} // end while
	
	ifs_header.close();
	
	ifstream ifs_pileup( Pileup_File.c_str() );
	string strPileup = "";
	string reference = "";
	string previous_reference = "";
	int* RefCoverage;
	while (ifs_pileup.good()) {
		getline(ifs_pileup, strPileup);
		
		if (strPileup.empty()) { // ignore emtry lines
			continue;
		} else {
			vector<string> entries;
			strsplit(strPileup, entries, "\t");
			
			string reference = entries[0];
			int pos = atoi(entries[1].c_str());
			
			if (reference.compare(previous_reference) != 0) {
				if (previous_reference.length() > 0) {
					//break;
					Coverage[previous_reference] = RefCoverage;
				}
				
				RefCoverage = Coverage[reference];
				previous_reference = reference;
			}
			
			int window = pos / WINDOW_SIZE;
			*(RefCoverage + window) = *(RefCoverage + window) + atoi(entries[3].c_str());
			
			//cout << reference << "\t" << pos << "\t" << window << "\t" << atoi(entries[3].c_str()) << "\t" << *(RefCoverage + window) << endl;
		}
	} // end while
	Coverage[previous_reference] = RefCoverage;
	
	ifs_pileup.close();
	
	ofstream ofs(OUT_File.c_str());
	map < string, int* >::iterator itCoverage = Coverage.begin();
	for (; itCoverage != Coverage.end(); itCoverage++) {
		reference = itCoverage->first;
		RefCoverage = itCoverage->second;
		
		int length = RefLength[reference];
		int nWindows = length / WINDOW_SIZE + 1;
		for (int i = 0; i < nWindows; i++) {
			
			int current_window_length = WINDOW_SIZE; 
			if ( i == nWindows - 1 ) { // last window is likely to be smaller
				current_window_length = length % (WINDOW_SIZE * i);
			}
			
			double avg_coverage = (double)*(RefCoverage + i)/(double)current_window_length;
			ofs << reference << "\t" << WINDOW_SIZE * i + 1 << "\t" << WINDOW_SIZE * (i + 1) << "\t" << avg_coverage << endl;
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
	string Pileup_File (argv[2]);
	string OUT_File (argv[3]);
	WINDOW_SIZE = atoi(argv[4]);

	CalculateCoverage(Header_File, Pileup_File, OUT_File);

	return 0;
}

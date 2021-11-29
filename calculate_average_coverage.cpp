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

void CalculateCoverage(string Header_File, string Pileup_File, string Coordinate_File) {
	
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
			
			// memory initialisation
			int *RefCoverage = new int[length];
			for(int i = 0; i < length; i++) {
				*(RefCoverage+i) = 0;
			}
			
			RefLength[reference] = length;
			Coverage[reference] = RefCoverage;
		}
	} // end while
	
	ifs_header.close();
	
	ifstream ifs_pileup( Pileup_File.c_str() );
	string strPileup = "";
	while (ifs_pileup.good()) {
		getline(ifs_pileup, strPileup);
		
		if (strPileup.empty()) { // ignore emtry lines
			continue;
		} else {
			vector<string> entries;
			strsplit(strPileup, entries, "\t");
			
			string reference = entries[0];
			int pos = atoi(entries[1].c_str());
			
			int* RefCoverage = Coverage[reference];
			
			*(RefCoverage+pos) = atoi(entries[3].c_str());
		}
	} // end while
	
	ifs_pileup.close();
	
	ifstream ifs_coordinate( Coordinate_File.c_str() );
	string strCoordinate = "";
	while (ifs_coordinate.good()) {
		getline(ifs_coordinate, strCoordinate);
		
		if (strCoordinate.empty()) { // ignore emtry lines
			continue;
		} else {
			strCoordinate.erase(remove(strCoordinate.begin(), strCoordinate.end(), '\r'), strCoordinate.end());
			strCoordinate.erase(remove(strCoordinate.begin(), strCoordinate.end(), '\n'), strCoordinate.end());

			vector<string> entries;
			strsplit(strCoordinate, entries, "\t");
			
			
			string GeneName = entries[0];
			string reference = entries[1];
			int start = atoi(entries[2].c_str()) - 1;
			int end = atoi(entries[3].c_str()) - 1;
			
			int* RefCoverage = Coverage[reference];

			double sum = 0.0;
			for (int i = start; i <= end; i++) {
				sum += *(RefCoverage+i);
			}

			cout << strCoordinate << "\t" << sum / (double)(end - start + 1) << endl;
		}
	} // end while
	
	ifs_coordinate.close();
	
}


int main (int argc, char** argv) {
	
	if (argc < 3) {
		cout << "Files not specified" << endl;
		exit(0);
	}
	
	string Header_File (argv[1]);
	string Pileup_File (argv[2]);
	string Coordinate_File (argv[3]);

	CalculateCoverage(Header_File, Pileup_File, Coordinate_File);

	return 0;
}

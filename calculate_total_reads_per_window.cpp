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

void CalculateReadCount(string Header_File, string SAM_File, string OUT_File, int Window_Size, int Mapping_Quality) {
	
	map < string, int* > ReadCount;
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
			
			int nWindows = length / Window_Size + 1;
			// memory initialisation
			int *RefReadCount = new int[nWindows];
			for(int i = 0; i < nWindows; i++) {
				*(RefReadCount+i) = 0;
			}
			
			RefLength[reference] = length;
			ReadCount[reference] = RefReadCount;
		}
	} // end while
	
	ifs_header.close();
	
	ifstream ifs_sam( SAM_File.c_str() );
	string strSAM = "";
	string reference = "";
	string previous_reference = "";
	int* RefReadCount;
	while (ifs_sam.good()) {
		getline(ifs_sam, strSAM);
		
		if (strSAM.empty()) { // ignore emtry lines
			continue;
		} else if (strSAM[0] =='@') { // ignore header
			continue;
		} else {
			vector<string> entries;
			strsplit(strSAM, entries, "\t");
			
			if (atoi(entries[4].c_str()) < Mapping_Quality) {
				continue;
			}
			string reference = entries[2];
			int pos = atoi(entries[3].c_str());
			
			if (reference.compare("*") == 0 || reference.compare("chloroplast") == 0 || reference.compare("mitochondria") == 0 || reference.compare("ChrC") == 0 || reference.compare("ChrM") == 0) {
				continue;
			}			

			if (reference.compare(previous_reference) != 0) {
				if (previous_reference.length() > 0) {
					//break;
					ReadCount[previous_reference] = RefReadCount;
				}
				
				RefReadCount = ReadCount[reference];
				previous_reference = reference;
			}
			
			int window = pos / Window_Size;
			*(RefReadCount + window) = *(RefReadCount + window) + 1;
			
			//cout << reference << "\t" << pos << "\t" << window << "\t" << atoi(entries[3].c_str()) << "\t" << *(RefCoverage + window) << endl;
		}
	} // end while
	ReadCount[previous_reference] = RefReadCount;
	
	ifs_sam.close();
	
	ofstream ofs(OUT_File.c_str());
	map < string, int* >::iterator itReadCount = ReadCount.begin();
	for (; itReadCount != ReadCount.end(); itReadCount++) {
		reference = itReadCount->first;
		RefReadCount = itReadCount->second;
		
		int length = RefLength[reference];
		int nWindows = length / Window_Size + 1;
		for (int i = 0; i < nWindows; i++) {
			ofs << reference << "\t" << Window_Size * i + 1 << "\t" << Window_Size * (i + 1) << "\t" << *(RefReadCount + i) << endl;
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
	string SAM_File (argv[2]);
	string OUT_File (argv[3]);
	int Window_Size (atoi(argv[4]));

	int Mapping_Quality = 0;
	if (argc == 6) {
		Mapping_Quality = atoi(argv[5]);
	}
	CalculateReadCount(Header_File, SAM_File, OUT_File, Window_Size, Mapping_Quality);

	return 0;
}

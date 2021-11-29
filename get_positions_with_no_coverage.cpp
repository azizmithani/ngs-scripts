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

void GetPositionsWithNoCoverage(string Pileup_File, string OUT_File) {
	
	ofstream ofs(OUT_File.c_str());
	ifstream ifs_pileup( Pileup_File.c_str() );
	string strPileup = "";
	string reference = "";
	string previous_reference = "";
	int previous_pos=0;
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
				previous_pos = 0;
				previous_reference = reference;
			}
			
			if (pos - previous_pos > 1) {
				for (int i = previous_pos + 1; i < pos; i++) {
					ofs << reference << "\t" << i << endl;
				}
			}
			previous_pos = pos;

			//cout << reference << "\t" << pos << "\t" << window << "\t" << atoi(entries[3].c_str()) << "\t" << *(RefCoverage + window) << endl;
		}
	} // end while
	
	ofs.close();
}


int main (int argc, char** argv) {
	
	if (argc < 3) {
		cout << "Files not specified" << endl;
		exit(0);
	}
	
	string Pileup_File (argv[1]);
	string OUT_File (argv[2]);

	GetPositionsWithNoCoverage(Pileup_File, OUT_File);

	return 0;
}

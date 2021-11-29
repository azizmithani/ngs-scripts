/*
 *
 *  Created by Aziz Mithani on 29/07/2010.
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

using namespace std;

typedef pair<long, string> RefINS;
typedef pair<string, RefINS > INS;


void ReadInsertionList (string Ins_File, set<INS> &InsList) {
	
	// open the insertion file stream
	ifstream ifs_ins ( Ins_File.c_str() );
	
	int offset = -1;
	string strIns = "";
	while (ifs_ins.good()) {
		getline(ifs_ins, strIns);
		
		if (strIns.empty()) { // ignore emtry lines
			continue;
		}
		
		vector<string> entries;
		strsplit(strIns, entries, "\t");
		
		if (offset < 0) {
			if (entries[2].compare(entries[0] + ":" + entries[1]) == 0) {
				offset = 1;
			} else {
				offset = 0;
			}
		}
		
		string reference = entries[0];
		long pos = atol(entries[1].c_str());
		string bases = entries[2+offset];

		RefINS refins (pos, bases);
		INS ins (reference, refins);
		
		
		InsList.insert(ins);
		
	}
	// close the SNP file
	ifs_ins.close();
		
}

int main (int argc, char** argv) {
	
	if (argc < 3) {
		cout << "Files not specified" << endl;
		exit(0);
	}
	
	string Ins_File_1 (argv[1]);
	set<INS> InsList1;
	ReadInsertionList(Ins_File_1, InsList1);

	for (int v = 2; v < argc; v++) {
		string Ins_File_2 (argv[v]);
		
		set<INS> InsList2;
		ReadInsertionList(Ins_File_2, InsList2);

		// allocate a vector for the differences
		vector<INS> InsList(InsList1.size());

		// find the common Ins
		vector<INS>::iterator it = set_intersection(InsList1.begin(), InsList1.end(), InsList2.begin(), InsList2.end(), InsList.begin());
		
		InsList1.clear();
		InsList1.insert(InsList.begin(), it);
		
	}
	
	set<INS>::iterator itIndel = InsList1.begin();
	for (; itIndel != InsList1.end(); itIndel++) {
		cout << itIndel->first << "\t" << itIndel->second.first << "\t" << itIndel->second.second << endl;
		//cout << itIndel->first << "\t" << itIndel->second.first << "\t" << itIndel->second.second << endl;
	}
	
	return 0;
}

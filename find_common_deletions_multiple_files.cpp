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

typedef pair<long, string> RefINDEL;
typedef pair<string, RefINDEL > INDEL;


void ReadIndelList (string Indel_File, set<INDEL> &IndelList) {
	
	// open the indel file stream
	ifstream ifs_indel ( Indel_File.c_str() );
	
	int offset = -1;
	string strIndel = "";
	while (ifs_indel.good()) {
		getline(ifs_indel, strIndel);
		
		if (strIndel.empty()) { // ignore emtry lines
			continue;
		}
		
		vector<string> entries;
		strsplit(strIndel, entries, "\t");
		
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

		RefINDEL refindel (pos, bases);
		INDEL indel (reference, refindel);
		
		
		IndelList.insert(indel);
		
	}
	// close the SNP file
	ifs_indel.close();
		
}

int main (int argc, char** argv) {
	
	if (argc < 3) {
		cout << "Files not specified" << endl;
		exit(0);
	}
	
	string Indel_File_1 (argv[1]);
	set<INDEL> IndelList1;
	ReadIndelList(Indel_File_1, IndelList1);

	for (int v = 2; v < argc; v++) {
		string Indel_File_2 (argv[v]);
		
		set<INDEL> IndelList2;
		ReadIndelList(Indel_File_2, IndelList2);

		// allocate a vector for the differences
		vector<INDEL> IndelList(IndelList1.size());

		// find the common Indel
		vector<INDEL>::iterator it = set_intersection(IndelList1.begin(), IndelList1.end(), IndelList2.begin(), IndelList2.end(), IndelList.begin());
		
		IndelList1.clear();
		IndelList1.insert(IndelList.begin(), it);
		
	}
	
	set<INDEL>::iterator itIndel = IndelList1.begin();
	for (; itIndel != IndelList1.end(); itIndel++) {
		cout << itIndel->first << "\t" << itIndel->second.first << "\t" << itIndel->second.second << endl;
		//cout << itIndel->first << "\t" << itIndel->second.first << "\t" << itIndel->second.second << endl;
	}
	
	return 0;
}

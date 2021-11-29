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

typedef pair<long, char> RefSNP;
typedef pair<string, RefSNP > SNP;


void ReadSNPList (string SNP_File, set<SNP> &SNPList) {
	
	// open the snp file stream
	ifstream ifs_snp ( SNP_File.c_str() );
	
	int offset = -1;
	string strSNP = "";
	while (ifs_snp.good()) {
		getline(ifs_snp, strSNP);
		
		if (strSNP.empty()) { // ignore emtry lines
			continue;
		}
		
		vector<string> entries;
		strsplit(strSNP, entries, "\t");
		
		if (offset < 0) {
			if (entries[2].compare(entries[0] + ":" + entries[1]) == 0) {
				offset = 1;
			} else {
				offset = 0;
			}
		}
		
		string reference = entries[0];
		long pos = atol(entries[1].c_str());
		string ConsensusBase = entries[3+offset];

		RefSNP refsnp (pos, ConsensusBase[0]);
		SNP snp (reference, refsnp);
		
		
		SNPList.insert(snp);
		
	}
	// close the SNP file
	ifs_snp.close();
		
}

int main (int argc, char** argv) {
	
	if (argc < 3) {
		cout << "Files not specified" << endl;
		exit(0);
	}
	
	string SNP_File_1 (argv[1]);
	set<SNP> SNPList1;
	ReadSNPList(SNP_File_1, SNPList1);

	for (int v = 2; v < argc; v++) {
		string SNP_File_2 (argv[v]);
		
		set<SNP> SNPList2;
		ReadSNPList(SNP_File_2, SNPList2);

		// allocate a vector for the differences
		vector<SNP> SNPList(SNPList1.size());

		// find the common snps
		vector<SNP>::iterator it = set_intersection(SNPList1.begin(), SNPList1.end(), SNPList2.begin(), SNPList2.end(), SNPList.begin());
		
		SNPList1.clear();
		SNPList1.insert(SNPList.begin(), it);
		
	}
	

	set<SNP>::iterator itSNP = SNPList1.begin();
	for (; itSNP != SNPList1.end(); itSNP++) {
		cout << itSNP->first << "\t" << itSNP->second.first << "\t.\t" << itSNP->second.second << endl;
		//cout << itSNP->first << "\t" << itSNP->second.first << "\t" << itSNP->second.second << endl;
	}
	
	return 0;
}

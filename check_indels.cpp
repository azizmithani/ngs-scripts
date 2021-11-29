/*
 *  check_snps.cpp
 *  
 *	Accept or reject SNPs based on the set criteria
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

using namespace std;

const int nCols = 5;

void ProcessIndelFile(string INDEL_File, string OUT_ACCEPT_File, string OUT_REJECT_File, int MinimumCoverage, double MinimumConsensus) {
	
	// open the indel file stream
	ifstream ifs_indel ( INDEL_File.c_str() );

	// also open the output file
	ofstream ofs_accept (OUT_ACCEPT_File.c_str());
	ofstream ofs_reject (OUT_REJECT_File.c_str());	
	
	int offset = 0;
	string strINDEL = "";
	while (ifs_indel.good()) {
		getline(ifs_indel, strINDEL);
		
		if (strINDEL.empty()) { // ignore emtry lines
			continue;
		}
		
		
		vector<string> entries;
		strsplit(strINDEL, entries, "\t");

		if (entries[2].compare(entries[0] + ":" + entries[1]) == 0) {
			offset = 1;
		} else {
			offset = 0;
		}

		cout << "Checking Indel at position " << entries[0] << ":" << entries[1] << endl;

		int coverage = atoi(entries[4+offset].c_str());

		//cout << "Checking Indel at position " << entries[0] << ":" << entries[1] << "\t" << coverage ;

		if ( coverage < MinimumCoverage ) {
			ofs_reject << strINDEL << endl;
			continue;
		}
		
		string reference = entries[0];
		int pos = atoi(entries[1].c_str());
		
		double consensus = atof(entries[5+offset].c_str());

		//cout << consensus << endl;
		
		if ( consensus < MinimumConsensus ) {
			ofs_reject << strINDEL << endl;
			continue;
		}
				
		// This snp has passed the filter ... accept the snp
		ofs_accept << strINDEL << endl;
		
	} // end while
	
	// close the files
	ifs_indel.close();
	ofs_accept.close();
	ofs_reject.close();
	
}



int main (int argc, char** argv) {
	
	if (argc < 2) {
		cout << "Parameters: INDEL_File<string>, Minimum_Coverage<int>, Minimum_Consensus_Percentage<double>" << endl;
		exit(0);
	}
	string INDEL_File (argv[1]);

	int MinimumCoverage;
	double MinimumConsensus;
	
	if (argc < 3) {
		MinimumCoverage = 20;
		MinimumConsensus = 90;
	} else if (argc < 4) {
		MinimumCoverage = atoi(argv[2]);
		MinimumConsensus = 90;
	} else {
		MinimumCoverage = atoi(argv[2]);
		MinimumConsensus = atof(argv[3]);
	}
	
	string OUT_ACCEPT_File = INDEL_File + ".acc";
	string OUT_REJECT_File = INDEL_File + ".rej";
	
	ProcessIndelFile(INDEL_File, OUT_ACCEPT_File, OUT_REJECT_File, MinimumCoverage, MinimumConsensus);
			
	return 0;
}

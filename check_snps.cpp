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

void ProcessSNPFile(string SNP_File, string BASE_DIST_File, string OUT_ACCEPT_File, string OUT_REJECT_File, int MinimumCoverage, double MinimumConsensus) {
	
	map < string, int** > Frequency;
	map < string, int > RefLength;
	
	ifstream ifs_dist ( BASE_DIST_File.c_str() );
	string str = "";
	while (ifs_dist.good()) {
		getline(ifs_dist, str);
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} else if (str.find("@SQ") != string::npos) { // get sequence length from the header
			vector<string> entries;
			strsplit(str, entries, "\t");
			
			string reference = entries[1].substr(3, entries[1].length());
			string strLength = entries[2].substr(3, entries[2].length());
			int length = atoi(strLength.c_str());
			
			// memory initialisation
			int **RefFrequency;
			RefFrequency = new int* [length];
			for(int i = 0; i < length; i++) {
				*(RefFrequency+i) = new int[nCols];		
				for (int j = 0; j < nCols; j++) {
					(*(RefFrequency+i))[j] = 0;
				}
			}
			
			RefLength[reference] = length;
			Frequency[reference] = RefFrequency;
		} else {
			vector<string> entries;
			strsplit(str, entries, "\t");
			
			string reference = entries[0];
			int pos = atoi(entries[1].c_str());
			
			int** RefFrequency = Frequency[reference];
			int* PosFreqency = *(RefFrequency + pos - 1);
			
			for ( int i = 3; i < entries.size(); i++ ) {
				PosFreqency[i-3] = atoi(entries[i].c_str());
			}
			
		}
	} // end while
	
	ifs_dist.close();
	
	// open the snp file stream
	ifstream ifs_snp ( SNP_File.c_str() );

	// also open the output file
	ofstream ofs_accept (OUT_ACCEPT_File.c_str());
	ofstream ofs_reject (OUT_REJECT_File.c_str());	
	
	int offset = 0;
	string strSNP = "";
	while (ifs_snp.good()) {
		getline(ifs_snp, strSNP);
		
		if (strSNP.empty()) { // ignore emtry lines
			continue;
		}
		
		
		vector<string> entries;
		strsplit(strSNP, entries, "\t");

		if (entries[2].compare(entries[0] + ":" + entries[1]) == 0) {
			offset = 1;
		} else {
			offset = 0;
		}

		cout << "Checking SNP at position " << entries[0] << ":" << entries[1] << " ... ";

		int coverage = atoi(entries[7+offset].c_str());

		if ( coverage < MinimumCoverage ) {
			ofs_reject << strSNP << endl;
			cout << "Rejected (Low Coverage) " << endl;
			continue;
		}
		
		string reference = entries[0];
		int pos = atoi(entries[1].c_str());
		
		if (pos >= RefLength[reference]) { // outside the reference .. reject it
			ofs_reject << strSNP << endl;
			cout << "Rejected (Outside Reference) " << endl;
			continue;
		}
		
		string ConsensusBase = entries[3+offset];
		int base;
		switch (ConsensusBase[0]) {
			case 'A':
				base = 0;
				break;
			case 'C':
				base = 1;
				break;
			case 'G':
				base = 2;
				break;
			case 'T':
				base = 3;
				break;
			default:
				base = 4;
				break;
		} // end switch
		
		if (base == 4) { //don't accept ambiguous characters
			ofs_reject << strSNP << endl;
			cout << "Rejected (Ambiguous Character) " << endl;
			continue;
		}
		
		int** RefFrequency = Frequency[reference];
		int* PosFreqency = *(RefFrequency + pos - 1);
	
		double consensus = (coverage == 0 ? 0.0 : (double) PosFreqency[base]/(double) coverage);
		
//		if (pos == 2380) {
//			cout << "Debug:\t" << consensus << "\t" << PosFreqency[base] << "\t" << coverage << endl;
//		}
		
		if ( consensus < MinimumConsensus ) {
			ofs_reject << strSNP << endl;
			cout << "Rejected (Low Consensus) " << endl;
			continue;
		}
				
		// This snp has passed the filter ... accept the snp
		ofs_accept << strSNP << endl;
		cout << "Accepted " << endl;
		
	} // end while
	
	// close the files
	ifs_snp.close();
	ofs_accept.close();
	ofs_reject.close();
	
}



int main (int argc, char** argv) {
	
	if (argc < 2) {
		cout << "Parameters: SNP_File<string>, BASE_Dist_File<string>, Minimum_Coverage<int>, Minimum_Consensus_Percentage<double>" << endl;
		exit(0);
	}
	if (argc < 2) {
		cout << "SNP File not specified" << endl;
		exit(0);
	} else if (argc < 3) {
		cout << "Base Distribution File not specified" << endl;
		exit(0);
	}
	string SNP_File (argv[1]);
	string BASE_DIST_File (argv[2]);

	int MinimumCoverage;
	double MinimumConsensus;
	
	if (argc < 4) {
		MinimumCoverage = 10;
		MinimumConsensus = 0.9;
	} else if (argc < 5) {
		MinimumCoverage = atoi(argv[3]);
		MinimumConsensus = 0.9;
	} else {
		MinimumCoverage = atoi(argv[3]);
		MinimumConsensus = atof(argv[4]);
	}
	
	string OUT_ACCEPT_File = SNP_File + ".acc";
	string OUT_REJECT_File = SNP_File + ".rej";
	
	ProcessSNPFile(SNP_File, BASE_DIST_File, OUT_ACCEPT_File, OUT_REJECT_File, MinimumCoverage, MinimumConsensus);
			
	return 0;
}

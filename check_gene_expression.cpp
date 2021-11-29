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

const double AVG_COVERAGE_THRESHOLD = 5.0;
const double PERCENTAGE_COVERAGE_THRESHOLD = 50.0;

void ReadCoverageFile (string Header_File, string Coverage_File, map < string, int* >& CoverageList) {
	
	map < string, long > RefLength;
	
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
			
			RefLength[reference] = length;
		}
	} // end while
	
	ifs_header.close();
	
	// open the coverage file stream
	ifstream ifs_Coverage ( Coverage_File.c_str() );
	
	string reference;
	string PreviousReference = "";
	string str="";
	int* RefCoverageList;
	while (ifs_Coverage.good()) {
		getline(ifs_Coverage, str);
		
		
		if (str.empty()) { // ignore emtry lines
			continue;
		} 
		
		vector<string> entries;
		strsplit(str, entries, "\t");
		
		reference = entries[0];
		long pos = atol(entries[1].c_str());
		int coverage = atoi(entries[2].c_str());
		
		if (PreviousReference.compare(reference) != 0) {
			if (PreviousReference.length() > 0) {
				CoverageList[PreviousReference] = RefCoverageList;
			}
			
			long length = RefLength[reference];
			// memory initialisation
			RefCoverageList = new int[length];
			for(long i = 0; i < length; i++) {
				*(RefCoverageList+i) = 0;		
			}
			CoverageList[reference] = RefCoverageList;
			
			PreviousReference = reference;
		}

		*(RefCoverageList + pos - 1) = coverage;
		
	}
	CoverageList[PreviousReference] = RefCoverageList;
	// close the SNP file
	ifs_Coverage.close();
		
}

void CheckGeneExpression(string GFF_File, map < string, int* >& CoverageList, string OUT_File) {
	// open the stream for output
	ofstream ofs (OUT_File.c_str());

	// gff file
	ifstream ifs_gff (GFF_File.c_str());
	string strGFF = "";
	string reference = "";
	while (ifs_gff.good()) {
		// read the next line from the gff file
		getline(ifs_gff, strGFF);
		
		if (strGFF.empty()) { // ignore emtry lines
			continue;
		} else if (strGFF[0] == '#') {
			continue;
		}
		
		vector<string> entries;
		strsplit(strGFF, entries, "\t");
		
		reference = entries[0];
		// get start and end coordinates
		long start = atol(entries[3].c_str());
		long end = atol(entries[4].c_str());
		// calculate the gene length
		double GeneLength = end - start + 1;

		int* RefCoverageList = CoverageList[reference];
		
		int TotalCoverage = 0;
		int PositionsWithCoverage = 0;
		for (int pos = start - 1; pos < end; pos++) {
			TotalCoverage += RefCoverageList[pos];
			if (RefCoverageList[pos] > 0) {
				PositionsWithCoverage++;
			}
		}
		
		double PercentageGeneCovered = (double)PositionsWithCoverage / (double)GeneLength * 100.00;
		
		double AvgCoverage = (double)TotalCoverage / (double)PositionsWithCoverage * 100.00;
				
		string ID = entries[8].substr(entries[8].find("ID=",0) + 3, entries[8].find_first_of(';') - 3);
				
		if (PercentageGeneCovered < PERCENTAGE_COVERAGE_THRESHOLD || AvgCoverage < AVG_COVERAGE_THRESHOLD) {
			cout << reference << "\t" << ID << "\t" << start << "\t" << end << "\t" << GeneLength << "\t" << PositionsWithCoverage << "\t" << PercentageGeneCovered << "\t" << AvgCoverage << "\tUNEXPRESSED" << endl;
		} else {
			cout << reference << "\t" << ID << "\t" << start << "\t" << end << "\t" << GeneLength << "\t" << PositionsWithCoverage << "\t" << PercentageGeneCovered << "\t" << AvgCoverage << "\tEXPRESSED" << endl;
		}
	} // end while ifs_gff is good

	ifs_gff.close();
	ofs.close();

	
}

int main (int argc, char** argv) {
	
	if (argc < 2) {
		cout << "SAM File not specified" << endl;
		exit(0);
	}
	string Header_File (argv[1]);
	string GFF_File (argv[2]);
	string Coverage_File (argv[3]); 
	string OUT_File (argv[4]);
	
	// variable to hold related lists 
	map < string, long > RefLength;
	map < string, int* > CoverageList;
	
	// Read Coverage
	ReadCoverageFile(Header_File, Coverage_File, CoverageList);
	
	// check gene expression
	CheckGeneExpression(GFF_File, CoverageList, OUT_File);
	
	return 0;
}


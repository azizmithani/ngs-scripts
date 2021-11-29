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


typedef map < pair<int, int>, pair< string, int > > REF_GFF_DATA;
typedef map < string, REF_GFF_DATA > GFF_DATA;

int GetReadEndCoordinates (int ReadStart, string cigar) {
	
	long sum = 0;
	string length = "";
	for (int i = 0; i < cigar.length(); i++) {
		if (isdigit(cigar[i])) {
			length += cigar[i];
		} else if (cigar[i]=='S') {
			// ignore
			length="";
		} else if (cigar[i] =='D') {
			sum+= atol(length.c_str());
			length="";
		} else if (cigar[i] =='I') {
			length = "";
		} else if (cigar[i]=='M' || cigar[i]=='N') {
			sum += atol(length.c_str());
			length="";
		}
	} // end for
	
	return ReadStart + sum - 1;
}

void ReadGFFFile(string GFF_File, GFF_DATA& GFFData) {
	
	ifstream ifs_gff ( GFF_File.c_str() );
	string strGFF = "";
	
	string reference = "";
	string previous_reference = "";
	REF_GFF_DATA RefGFFList;
	while (ifs_gff.good()) {
		// read the next line from the gff file
		getline(ifs_gff, strGFF);
		
		if (strGFF.empty()) { // ignore emtry lines
			continue;
		}
		if (strGFF[0] == '#') { // ignore header
			continue;
		}
		
		vector<string> entries;
		strsplit(strGFF, entries, "\t");
		
		reference = entries[0];
		int start = atoi(entries[3].c_str());
		int end = atoi(entries[4].c_str());
		string ID = entries[8].substr(entries[8].find("ID=",0) + 3, entries[8].find_first_of(';') - 3);
		
		if (reference.compare(previous_reference) != 0) {
			if (previous_reference.length() > 0) {
				GFFData[previous_reference] = RefGFFList;
			}
			
			// get the GFF list for this reference
			RefGFFList = GFFData[reference];
			
			// save the current reference as previous reference
			previous_reference = reference;
		}
		
		RefGFFList[pair<int, int>(start,end)] = pair<string, int>(ID, 0);
	} // end while
	// update the last reference
	GFFData[reference] = RefGFFList;
	
}

int ProcessSAMFile(string SAM_File, GFF_DATA& GFFData){//, map<string, int>& RefLength) {
	
	// open the sam file stream
	ifstream ifs_sam ( SAM_File.c_str() );
		
	string strSAM = "";
	string Read = "";
	string reference = "";
	string previous_reference = "";
	int TotalReads = 0;
//	bool GetNextRead = true;

	REF_GFF_DATA RefGFFData;
	REF_GFF_DATA::iterator itRefGFFData;

	while (ifs_sam.good()) {
//		if (GetNextRead) {
			getline(ifs_sam, strSAM);
//		}
		
//		GetNextRead = true;
		if (strSAM.empty()) { // ignore emtry lines
			continue;
		} else if (strSAM[0] == '@') { // ignore the header
/*AMT
			vector<string> entries;
			strsplit(str, entries, "\t ");
			
			reference = entries[1].substr(3, entries[1].length());
			string strLength = entries[2].substr(3, entries[2].length());
			int length = atoi(strLength.c_str());
			
			RefLength[reference] = length;
AMT*/
			continue;
		}

		bool invalidGene = false;
		// get individual fields
		vector<string> entries;
		strsplit(strSAM, entries, "\t");
				
		// get the reference names
		reference = entries[2];
		
		// read coordinates
		int ReadStart = atoi(entries[3].c_str());
		int ReadEnd = GetReadEndCoordinates(ReadStart, entries[5]);
		
		// get the GFF data for this reference
		if (reference.compare(previous_reference) != 0) {
			
			if (previous_reference.length() > 0) {
				// save the GFF data for the previous reference
				GFFData[previous_reference] = RefGFFData;
			}			
			
			// get the new data
			RefGFFData = GFFData[reference];
			itRefGFFData = RefGFFData.begin();
			
			while (itRefGFFData != RefGFFData.end() && ReadStart > itRefGFFData->first.second) {
				itRefGFFData++;
			}
			
			if (itRefGFFData == RefGFFData.end() || ReadEnd < itRefGFFData->first.first) {
				invalidGene = true;
			}
			
			previous_reference = reference;
		}		
		
		if (ReadStart > itRefGFFData->first.second) {
			// read lies after the current gene			
						
//			GetNextRead = false;
			while (itRefGFFData != RefGFFData.end() && ReadStart > itRefGFFData->first.second) {
				itRefGFFData++;
			}
			
			if (itRefGFFData == RefGFFData.end() || ReadEnd < itRefGFFData->first.first) {
				invalidGene = true;
			}
						
//			continue;
		}
		
		if (!invalidGene) 
			itRefGFFData->second.second++;
		
		TotalReads++;
		
	} // end if sam file is good
		
	// save the GFF data for the last reference
	GFFData[previous_reference] = RefGFFData;
	
	ifs_sam.close();
	
	return TotalReads;
}


int main (int argc, char** argv) {
	
	if (argc < 3) {
		cout << "SAM File not specified" << endl;
		exit(0);
	}
	string GFF_File (argv[1]);
	string SAM_File (argv[2]); 
//	string OUT_File (argv[4]);
	
	// variable to hold related lists 
	map < string, long > RefLength;

	// read the GFF file
	GFF_DATA GFFData;
	ReadGFFFile(GFF_File, GFFData);
	
	// Process the SAM file
	int TotalReads = ProcessSAMFile(SAM_File, GFFData);
	
	//cout << "Total Reads:\t" << TotalReads << endl;
	GFF_DATA::iterator itGFFData = GFFData.begin();
	for (; itGFFData != GFFData.end(); itGFFData++) {
		REF_GFF_DATA RefGFFData = itGFFData->second;
		
		REF_GFF_DATA::iterator itRefGFFData = RefGFFData.begin();
		for (; itRefGFFData != RefGFFData.end(); itRefGFFData++) {
			int GeneLength = itRefGFFData->first.second - itRefGFFData->first.first + 1;
			
			double RPKM = (1000000000.00 / (double)(TotalReads)) * ((double)itRefGFFData->second.second / (double)(GeneLength));
			cout << itRefGFFData->second.first << "\t" << GeneLength << "\t" << itRefGFFData->second.second << "\t" << RPKM << endl;
		}
	}
		
	return 0;
}

